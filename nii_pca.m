function nii_pca (xlsname, numFactor)
%compute Principle Components Analysis on NiiStat Behavioral Data
%This can allow you to conduct an analysis on a reduced number of factors
% xlsname : name of design file to read
% numFactor : number of PCA factors to compute (optional, if not specified all)
%Examples
% nii_pca('pca43.xlsx');
% nii_pca('pca43.xlsx', 2); %only compute 2 factors

%xlsname = 'pca430.xlsx';
if ~exist('xlsname','var')  
   [file,pth] = uigetfile({'*.xls;*.xlsx;*.txt;*.tab','Excel/Text file';'*.txt;*.tab','Tab-delimited text (*.tab, *.txt)';'*.val','VLSM/NPM text (*.val)'},'Select the design file'); 
   if isequal(file,0), return; end;
   xlsname=[pth file];
end
designMat = nii_read_design(xlsname);
[p,n] = fileparts(xlsname);
statname = fullfile(p,[n,'_',timeStampSub(),'_pca.tab']);
diary (statname);
[numMat, designMat] = onlyNumericSub(designMat);


numCond = size(numMat,2);
fprintf('PCA for %d subjects and %d conditions\n', size(numMat,1), numCond);
numMat = normColSub(numMat);
%compute PCA
%http://stats.stackexchange.com/questions/104274/factor-analysis-for-given-data-with-help-of-matlab
means=mean(numMat);
centered=numMat-repmat(means,size(numMat,1),1);
covariance=(centered'*centered)/(9);
[V,D]=eig(covariance);
[e,i]=sort(diag(D),'descend');
sorted=V(:,i);
if ~exist('numFactor','var')
    numFactor = numCond;
    fprintf('Estimating ALL factors as conditions! (reduce manually)\n'); 
end
if (numFactor <1) or (numFactor > numCond)
    error('Please specify between 1 and %d factors', numCond);
end
factorLoadings=sorted(:,1:numFactor);
PCA=numMat*factorLoadings;
%report results
%NEXT: report variability predicted by each factor
pct = (e./sum(e))*100;
for c = 1: numFactor
    fprintf('Factor\t%d\tPct\t%g\tPctTotal\t%g\n',c,pct(c), sum(pct(1:c)));
end
%NEXT: report FACTOR LOADINGS
SNames = fieldnames(designMat);
fprintf('\nFactorLoadings\n');
fprintf('Task\t');
for f = 1: numFactor
   fprintf('Factor%d\t', f); 
end
fprintf('\n');
for c = 1: numCond
    fprintf('%s\t', SNames{c+1} );
    for f = 1: numFactor
        fprintf('%g\t',factorLoadings(c,f));
    end;
    fprintf('\n');
end
%Next: report PCA
fprintf('\nSubjectID\t');
for f = 1: numFactor
   fprintf('Factor%d\t', f); 
end
fprintf('\n');
numSubj = size(numMat,1);
for s = 1: numSubj
    fprintf('%s\t', deblank( designMat(s).(SNames{1})));
    for f = 1: numFactor
       fprintf('%g\t', PCA(s,f)); 
    end
    fprintf('\n');
end
diary off
fprintf('Saved as: %s\n', statname);
%end nii_pca


function [numMat, designMat] = onlyNumericSub(designMat)
SNames = fieldnames(designMat);
numCond = numel(SNames)-1;
numSubj = numel(designMat);
if numCond < 2 || numSubj < 4
    error('PCA requires at least two conditions and 4 subjects.');
end
%remove rows with nonsense data
numMat = zeros(numSubj,numCond);
numMat(:) = nan;
for s=1:numSubj
    sname = deblank( designMat(s).(SNames{1}));
    for c=1:numCond
        v =   designMat(s).(SNames{c+1});
        if ~isnumeric(v)
            fprintf('Excluding %s: "%s" is not a number\n',sname, v);
            break
        else
            numMat(s,c) = v;
        end
    end
end
rejectRows =any(isnan(numMat),2);
numMat(rejectRows,:)=[];
designMat(rejectRows)=[];
if size(numMat,1) < 4
    error('PCA requires at least 4 subjects.');
end
%end onlyNumericSub()

function datetime = timeStampSub
datetime=datestr(now);
datetime=strrep(datetime,':',''); %Replace colon with underscore
datetime=strrep(datetime,'-','');%Replace minus sign with underscore
datetime=strrep(datetime,' ','_');%Replace space with underscore
%end timeStampSub()

function x = normColSub(x)
%normalize each column for range 0..1
% x = [1 4 3 0; 2 6 2 5; 3 10 2 2] -> x = [0 0 1 0; 0.5 0.333 0 1; 1 1 0 0.4]
if size(x,1) < 2 %must have at least 2 rows
    fprintf('Error: normalizing columns requires multiple rows\n');
    return
end
if min(max(x,[],1)-min(x,[],1)) == 0
    fprintf('Error: unable to normalize columns: some have no variability\n');
    return;
end
x = bsxfun(@minus,x,min(x,[],1)); %translate so minimum = 0
x = bsxfun(@rdivide,x,max(x,[],1)); %scale so range is 1
%end normColSub()