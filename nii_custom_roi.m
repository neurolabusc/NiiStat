function nii_custom_roi (xlsname)
%compute Principle Components Analysis on NiiStat Behavioral Data
%This can allow you to conduct an analysis on a reduced number of factors
% xlsname : name of design file to read
% numFactor : number of PCA factors to compute (optional, if not specified all)
%Examples
% nii_pca('pca43.xlsx');
% nii_pca('pca43.xlsx', 2); %only compute 2 factors

clc
xlsname = 'P50x.xlsx';
modality = 'lesion';
region = 6; %e.g. 1=aal 2=aalcat 3=bro 4=cat 5=fox 6=jhu 
%regionIndices = [11 13 15]; rName = 'BrocaFrac';%Brocas
regionIndices = [35]; rName = 'WernickeFrac';%Wernickes
%regionIndices = [11 15 184]; rName = 'DSdFrac';  %DSd
%regionIndices = [186 188 41]; rName = 'DSvFrac'; %DSv




[kROIs, kROInumbers] = nii_roi_list();
roi = [deblank(kROIs(6,:)), '.nii'];
%next check inputs
if ~exist('xlsname','var')  
   [file,pth] = uigetfile({'*.xls;*.xlsx;*.txt;*.tab','Excel/Text file';'*.txt;*.tab','Tab-delimited text (*.tab, *.txt)';'*.val','VLSM/NPM text (*.val)'},'Select the design file'); 
   if isequal(file,0), return; end;
   xlsname=[pth file];
end
if ~exist(xlsname,'file')  
    error('Unable to find %s', xlsname);
end
if ~exist(roi,'file')  
    error('Unable to find %s', roi);
end
hdrROI = spm_vol(roi);
imgROI = spm_read_vols(hdrROI);
imgROI = ismember(imgROI,regionIndices);
nROI = sum(imgROI(:));
if nROI < 1
   error('No ROIs survive %s', roi); 
end
fprintf('%g voxels survive in [ ', nROI);
fprintf('%g ',regionIndices(:))
fprintf('] of %s\n', roi);


designMat = nii_read_design(xlsname);
numSubj = numel(designMat);
SNames = fieldnames(designMat);
SNames = SNames{1};
%[outhdr, outimg] = nii_reslice_target(inhdr, inimg, tarhdr, interp) 
fprintf('%s\n',rName);
for s = 1 : numSubj
    sname = deblank( designMat(s).(SNames));
    m = load(sname);
    img = m.(modality).dat;
    hdr = m.(modality).hdr;
    [hdrX, imgX] = nii_reslice_target(hdr, img, hdrROI, 1);
    nTotal = sum(imgX(:));
    imgX(imgROI == 0) = 0;
    nMask = sum(imgX(:));
    
    fprintf('%s\ttotal\t%g\tmask\t%g\tfrac\t%g\n', sname, nTotal, nMask, nMask/nTotal);
    %fprintf('%g\n', nMask/nTotal);
    
end


