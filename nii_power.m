function nii_power(fnm, zCrit, pPower, flipR)
%Compute number of people required to replicate observed effects
%(optional) inputs:
% fnm : filename of power.mat file created by NiiStat
% zCrit : threshold, e.g. 1.64 for 1-tailed p<0.05 uncorrected for FWE
% pPower : 0..1 Probability of rejecting the null hypothesis (1-Beta)
%Note
% To create a power.mat file type the following in the Matlab command line
% BEFORE running NiiStat, type in command line:
%   global global_powerMap; global_powerMap = true;
% Then run 
%   NiiStat
% The analysis will create a file named "power.mat" you can analyze with
% this script
%
% http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-1-Sided
% http://www.sample-size.net/correlation-sample-size/
%Note:
% You can determine zCrit with spm_invNcdf(), for example for 0.05%
%  zCrit = spm_invNcdf(1-0.05); %e.g. for alpha = 0.05 this would be 1.6449, as zCrit = spm_invNcdf(1-alpha)
%Example
% nii_power; %use GUI
% nii_power('power.mat', 3.5, 0.6, 1)

if (nargin < 4)  
    prompt = {'Enter threshold Z-score:','Power (0..1, e.g. 0.7=70% of studies will replicate)'};
    dlg_title = 'Power analysis preferences';
    num_lines = 1;
    def = {'3.5','0.6'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer), return; end;
    zCrit = str2double(answer{1});
    pPower = str2double(answer{2});
end;
zCrit = abs(zCrit);
if ~exist('fnm','var'), fnm = 'power.mat'; end;
if ~exist(fnm,'file') 
   [A,Apth] = uigetfile({'*.mat;';'*.*'},'Select first image');
   fnm = [Apth, A];
end;
m = load(fnm);
%check inputs
nObs = size(m.les,1);
nVox = size(m.les,2);
nBeh = size(m.beh,2);
fprintf('%d participants, %d voxels/regions tested, %d behaviors tested\n', nObs, nVox, nBeh);
if size(m.beh,1) ~= nObs, error('behavioral data must have one entry per participant.\n'); end;
if size(m.beh,2) ~= 1, error('%s only currently supports one behavior\n', mfilename); end;
fprintf('Threshold z = %g, Power = %g, beta (type2 errors) = %g\n', zCrit, pPower, 1-pPower);
if (sum(m.les(:,1)==0)+sum(m.les(:,1)==1) ) == numel(m.les(:,1)) 
    powerTTestSub(m, zCrit, pPower); %lesions are binomial: each voxel either spared or destroyed
else
    powerRegressionSub(m, zCrit, pPower); %lesions are continuous, e.g. ROI or fMRI
end
%end nii_power()

function powerTTestSub(m, zCrit, pPower)
%http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-1-Sided
%the commented code below replicates means of 1 and 0, std of 1 and 1 and sampling ratio of 1
% rnd = [0.8 1.2 -1.1 0.2 0 -1.4 0.4 -0.9 1.5 -0.7]; %mean=0, std=1
% m.beh = [rnd (rnd + 1)]'; %effect size is 1
% m.les = [zeros(numel(rnd),1); ones(numel(rnd),1)];
% zCrit = abs(spm_invNcdf(0.05)); %  p < 0.05
% pPower = 0.8;
% powerTTestSub(m, zCrit, pPower)
fprintf('Note: with t-test effect at every voxel/region assumed in predicted direction (predicted direction ignored).\n');
nObs = size(m.les,1);
nVox = size(m.les,2);
nBeh = size(m.beh,2);
%compute cohen's D
d = zeros(nVox,1);
kappa = zeros(nVox,1);
startTime = tic;
for i = 1: nVox
   [d(i), kappa(i)] = cohensD(m.les(:,i),m.beh); 
end
%compute power
pow = zeros(nVox,1);
invBeta = spm_invNcdf(pPower);
for i = 1: nVox
   pow(i)= powerFX(d(i), zCrit, invBeta, kappa(i));
end
if (~isfield(m, 'hdr'))
    for i = 1: nVox
        fprintf('Test %d Sample size (nA+nB) = %d\n', i, pow(i));
    end
    return;
end;
pow(isnan(pow)) = 0; %e.g. kappaR is 1/0
saveStatMap(m.hdr, m.roi_names, m.good_idx, d, 'CohensD');
saveStatMap(m.hdr, m.roi_names, m.good_idx, kappa, 'kappa');
saveStatMap(m.hdr, m.roi_names, m.good_idx, pow, 'Power');
fprintf('required %g seconds\n', toc(startTime));
%end powerTTestSub()

function [d, kappa] = cohensD(Y,X)
sdpooled = sqrt((std(X((Y==0),1))^2+std(X((Y==1),1))^2)/2);
effectSize = abs(mean(X((Y==0),1))-mean(X((Y==1),1)));
d = effectSize/sdpooled;
n1 = sum(Y(:));
n0 = numel(Y) - n1;
if n1 > n0 %avoid divide by zero, kappa will range from 1 (equal group size) to 0 (no people in smaller grouper)
    kappa = n0/n1;
else
   kappa = n1/n0; 
end
%end cohensD()

function powerRegressionSub(m, zCrit, pPower)
%http://www.sample-size.net/correlation-sample-size/
%compute correlation coefficient r
% m.les = [-0.9 -1.1 -1.4 0 0.8 -0.7 1.2 0.2 0.4 1.5]'; 
% m.beh = [1:numel(m.les)]';
% corrcoef(m.les, m.beh) %correlation = 0.7743
% zCrit = abs(spm_invNcdf(0.025)); %  p < 0.05
% pPower = 0.8;
% powerRegressionSub(m, zCrit, pPower,0);
nObs = size(m.les,1);
nVox = size(m.les,2);
nBeh = size(m.beh,2);
r = zeros(nVox,1);
startTime = tic;
for i = 1: nVox 
   cc = corrcoef(m.les(:,i),m.beh);
   r(i) = cc(2);
end
%if flipR 
%    r = -r;
%    fprintf('Inverting correlations: predicting negative correlations\n');
%end
%fprintf('Range of r (correlation coefficient) is %g..%g\n', min(r(:)), max(r(:)) );
r = abs(r);
fprintf('Max r = %g, absolute values (assuming each voxel/region in predicted direction)', max(r(:)));
pow = zeros(nVox,1);
invBeta = spm_invNcdf(pPower);
for i = 1: nVox
    c = 0.5 * log((1+r(i))/(1-r(i))); 
   pow(i)= round(((zCrit+invBeta)/c)^2 + 3);
end
pow(isnan(pow)) = 0; %e.g. kappaR is 1/0
if (~isfield(m, 'hdr'))
    for i = 1: nVox
        fprintf('Test %d r %g zAlpha %g zBeta %g Sample size (nA+nB) = %g\n', i, r(i), zCrit, invBeta, pow(i));
    end
    return;
end;
saveStatMap(m.hdr, m.roi_names, m.good_idx, pow, 'Power');
saveStatMap(m.hdr, m.roi_names, m.good_idx, r, 'R');
fprintf('required %g seconds\n', toc(startTime));
%end powerRegressionSub()

function saveStatMap(hdr, roi_names, good_idx, masked_img, statname)
if (isempty(hdr)) || (numel(hdr.dim) ~= 3), error('Header corrupt'); end;
if isempty(roi_names) %voxelwise analysis
    img = zeros(hdr.dim(1), hdr.dim(2), hdr.dim(3));
    img(good_idx) = masked_img(:);
else
    valIndex = zeros(numel(roi_names));
    valIndex(good_idx) = masked_img(:);
    if ~exist(hdr.fname,'file'), error('Unable to find %s',hdr.fname); end;
    %hdr = spm_vol(hdr.fname);
    imgIndex = spm_read_vols(hdr);
    img = zeros(size(imgIndex));
    for i = 1: numel(roi_names)
        img(imgIndex == i) = valIndex(i);
    end

end
%if numel(img) ~= prod(m.hdr.dim), 
%    error('unable to save image: expected %d voxels not %d', prod(m.hdr.dim), numel(img)); 
%end;
%img = reshape(img,m.hdr.dim);
hdr.fname = [deblank(statname) '.nii'];
hdr.pinfo = [1;0;0];
hdr.private.dat.scl_slope = 1;
hdr.private.dat.scl_inter = 0;
hdr.private.dat.dtype = 'FLOAT32-LE';%'INT16-LE', 'FLOAT32-LE';
hdr.dt    =[16,0]; %4= 16-bit integer; 16 =32-bit real datatype
spm_write_vol(hdr,img);
%end saveSumMapROI

function nTotal = powerFX(effectSize, zCrit, invBeta, kappa)
nA = ceil((1+1/kappa)*((zCrit+invBeta)/(effectSize))^2);
nB = ceil(kappa * nA);
nTotal = nA + nB;
%end powerFX()

function [nA, nB] = powerX(effectSize, alpha, beta, kappa)
%reports number needed for group A and B to achieve desired power
% effectSize: Cohen's effect size 
% alpha : p-value threshold
% beta : power requested
% kappa : sampling ratio of nB/nA, 1=equal numbers, 0.5 = nA twice as common
%http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-1-Sided
%nA=(sdA^2+sdB^2/kappa)*((qnorm(1-alpha)+qnorm(1-beta))/(muA-muB))^2)
nA = ceil((1+1/kappa)*((spm_invNcdf(1-alpha)+spm_invNcdf(beta))/(effectSize))^2);
nB = ceil(kappa * nA);