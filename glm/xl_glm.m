function xl_glm(xlsname, numDV, contrast, workSheetName)
%conduct statistics on Excel worksheet 
% xlsname : excel format file to analyze
% numDV : number of dependent variables: e.g. brain regions
% contrast : statistical contrast: if not provided user will be asked
% workSheetName : name of Excel worksheet to analyze (if empty "GLM" [case sensitive])
%Example Worksheet with 2 DVs (ROI1 and ROI2 are two brain regions
% ID ROI1 ROI2 PATIENT CONTROL
% S1  123  234 0       1
% S2  111  222 0       1
% S3  119  212 0       1
% S4  106  228 1       0
% S5  107  242 1       0
% S6  108  202 1       0
%Examples
% xl_glm; %use GUI
% xl_glm('kk.xlsx', 6, [1 -1 0]); %t-contrast 
% xl_glm('kk.xlsx', 6, [1 -1], 'GLM_NoNuisance');

fprintf('Version 30 June 2016 of %s %s %s\n', mfilename, computer, version);
ignoreColumnOne = true; %e.g. 1st column is subject ID
if ~exist('xlsname','var')
   [file,pth] = uigetfile({'*.xls;*.xlsx;*.txt;*.tab','Excel/Text file';'*.txt;*.tab','Tab-delimited text (*.tab, *.txt)';'*.val','VLSM/NPM text (*.val)'},'Select the design file');
   if isequal(file,0), return; end;
   xlsname=[pth file];
end
if exist(xlsname,'file') ~= 2
    error('%s Unable to find Excel file named %s\n', mfilename, xlsname);
end
if exist('xl_xls2mat','file') ~= 2
    error('%s requires xl_xls2mat to be installed', mfilename);
end
if exist('spm_Tcdf','file') ~= 2
    error('%s requires SPM to be installed', mfilename);
end;
if ~exist('workSheetName','var')
    workSheetName = 'GLM';
end;
%read file
designMat = xl_xls2mat(xlsname, workSheetName,[],true);
if numel(designMat) < 2
   error('%s unable to find multiple subjects in %s', mfilename, xlsname);
end
SNames = fieldnames(designMat);
if ignoreColumnOne
    designMat = rmfield(designMat,SNames{1}); %remove first column - mat name
    SNames = fieldnames(designMat);
end
if numel(SNames) < 2
   error('%s not enough columns %s ', mfilename, xlsname);
end
%check number of DVs
if ~exist('numDV','var')
   %show columnes so user can determine number of DVs
   fprintf('Columns:\n');
   for i = 1 : numel(SNames)
    fprintf(' %d\t%s\n', i, SNames{i});
   end
   answer = inputdlg(sprintf('Number of Dependent Variables (1..%d)', numel(SNames)-1), ...
       'Specify number of DVs (e.g. brain regions)', 1,  {'1'});
   numDV = str2num(answer{1});
end
if (numDV < 1) || (numDV >= numel(SNames))
   error('Number of Dependent variables must be 1..%d', numel(SNames)-1);
end
%extract values
numObs = numel(designMat); %number of observations (participants)
dat = zeros(numObs, numel(SNames));
for i = 1:numObs
    for j = 1:numel(SNames)
        v = designMat(i).(SNames{j});
        if isempty(v) || ~isnumeric(v) || ~isfinite(v)
            error('All cells must be numeric: unable to understand "%s"', v);
        end
        dat(i,j) = v;
    end %for j,  each column
end %for i, each row/observation/participant
DV = dat(:, 1:numDV);
IV = dat(:, numDV+1:end);
numIV = size(IV,2);

fprintf('Independent Variables:\n');
for i = numDV+1 : numel(SNames)
    fprintf(' %d\t%s\n', i - numDV, SNames{i});
end

if numIV == 1
    contrast = [1];
else
    if ~exist('contrast','var')
        contrast = [1 -1]; %assume t-test
        contrast = contrast(:); %ensure column vector
        contrast = [contrast; zeros(numIV - numel(contrast),1)]; %zero-pad
        answer = inputdlg('Statistical Contrast', 'Specify contrast', 1,  {sprintf('%d ', contrast)});
        contrast = str2num(answer{1});
    end
    contrast = contrast(:); %ensure column vector
    contrast = [contrast; zeros(numIV - numel(contrast),1)]; %zero-pad
end
fprintf('Contrast: [%s]\n', strtrim(sprintf('%d ', contrast)));
nPerms = 10000;
fprintf('Permutations: %d\n', nPerms);
%compute/display one-tailed results
fprintf('Results: (one-tailed)\n');
[pUnc pFWE peak volu mass t df] = glm_perm_fl(DV, IV, contrast, nPerms, 0.001);
for i = 1:numDV
    p =  1.0 - spm_Tcdf( t(i), df); %classic distribution
    fprintf(' %s\tt(%d)=\t%g\tp<\t%0.8f\tpPerm<\t%0.8f\tpFWE<\t%0.8f\n', ...
        SNames{i}, df, t(i), p, pUnc(i), pFWE(i));
end
%compute/display two-tailed results
fprintf('Results: (two-tailed)\n');
[pUnc pFWE peak volu mass t df] = glm_perm_fl(DV, IV, contrast, nPerms, 0.001, true);
for i = 1:numDV
    p = 2 * (1.0 - spm_Tcdf( abs(t(i)), df)); %classic distribution: n.b. ABS for two tail predicted!
    fprintf(' %s\tt(%d)=\t%g\tp<\t%0.8f\tpPerm<\t%0.8f\tpFWE<\t%0.8f\n', ...
        SNames{i}, df, t(i), p, pUnc(i), pFWE(i));
end


%end xl_glm()

%function [uncZ, threshMin, threshMax, t, df] = glm_perm_flSub(Y, X, c, nPerms, kPcrit, good_idx)
function [pUnc pFWE peak volu mass t df] = glm_perm_fl(Y, X, c, nPerms, pClus, is2Tail)
%glm_perm_fl: A simple Freedman-Lane permutation test for a t-contrast
% Usage: [pUnc pFWE peak volu mass t] = glm_perm_fl(Y, X, c, nPerms, pClus)
%
% peak, volu and mass are all set-level statistics (thus not requiring
% computation of contiguous clusters)

% Copyright 2010 Ged Ridgway
% http://www.mathworks.com/matlabcentral/fileexchange/authors/27434

% Defaults
if ~exist('is2Tail', 'var')
    is2Tail = false;
end
if nargin < 5
    pClus = 0.001;
end
if nargin < 4
    nPerms = 10000;
end

% Basics and reusable components
[n p] = size(X);
[nY v] = size(Y);
if nY ~= n, error('Size of data and design are inconsistent'); end

c0 = eye(p) - c * pinv(c);
X0 = X * c0;
R0 = eye(n) - X0 * pinv(X0);
Y  = R0 * Y; % pre-regress, but keep whole design below, for Freedman-Lane
df = n - rank(X);
tClus = spm_invTcdf(1 - pClus, df);
pXX = pinv(X)*pinv(X)'; % = pinv(X'*X), which is reusable, because
pX  = pXX * X';         % pinv(P*X) = pinv(X'*P'*P*X)*X'*P' = pXX * (P*X)'
% Things to track over permutations
pUnc = ones(1, v); % count times each voxel exceeded or equalled
peak = nan(1, nPerms);
volu = nan(1, nPerms);
mass = nan(1, nPerms);
% Original design (identity permutation)
t = glm_quick_t(Y, X, pXX, pX, df, c);
if is2Tail, t = abs(t); end;
peak(1) = max(t);
volu(1) = sum(t >= tClus);
mass(1) = sum(t(t >= tClus));
% Permutations
%h = waitbar(0, 'Permuting');
for p = 2:nPerms
    Xp  = X(randperm(n), :);
    pXp = pXX * Xp'; % = pinv(Xp)
    tp  = glm_quick_t(Y, Xp, pXX, pXp, df, c);
    if is2Tail, tp = abs(tp); end;
    pUnc = pUnc + (tp >= t);
    peak(p) = max(tp);
    volu(p) = sum(tp >= tClus);
    mass(p) = sum(tp(tp >= tClus));
    %waitbar(p/nPerms, h);
end
%close(h)
% Comput pUnc and pFWE
pUnc = pUnc / nPerms;
pFWE = nan(1, v);
for i = 1:v
    pFWE(i) = sum(peak >= t(i)) / nPerms;
end
%glm_perm_fl

function s = glm_quick_t(y, X, pXX, pX, df, c)
b = pX * y;                     % parameters
y = y - X * b;                  % residuals
s = sum(y .* y);                % sum squared error
s = sqrt(s .* (c'*pXX*c) / df); % standard error
s = c'*b ./ s;                  % t-statistic
%end glm_quick_t()