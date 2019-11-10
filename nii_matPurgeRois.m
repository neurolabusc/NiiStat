function nii_matPurgeRois(matnames)
%NiiStat mat files cache ROI means for fast access. This script removes this cache
%Use this if you want update a ROI, for example using ROI in MNI vs SPM space
% matname = (optional) filename(s) to purge
%Example
% nii_mat2nii; %use GUI
% nii_mat2nii('vx.mat');
% nii_mat2nii(strvcat('LM1001.mat','P051.mat'))

if ~exist('matnames','var') %no files specified
   [A,Apth] = uigetfile({'*.mat;';'*.*'},'Select image(s) to purge','MultiSelect', 'on'); 
   matnames = strcat(Apth,char(A));
end

if length(matnames) < 1, return; end
for f = 1 : size(matnames,1)
    matname = deblank(matnames(f,:));
    if ~exist(matname, 'file'), error('Unable to find %s', matname); end
    mat = load(matname);
    %auto-detect modalities
    fld=fieldnames(mat);
    nROI = 0;
    for i = 1: numel(fld)
        if ~isfield( mat.(fld{i}),'label'), continue; end;
        mat = rmfield(mat,fld{i});
        nROI = nROI + 1;
    end
    if (nROI < 1)
       fprintf('No ROIs to purge from %s\n', matname)
       continue; 
    end
    fprintf('Purged %d modalities*ROIs from %s\n', nROI,matname);
    save(matname,'-struct', 'mat');
end
