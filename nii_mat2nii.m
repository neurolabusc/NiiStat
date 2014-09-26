function nii_mat2nii (matname)
%Convert NiiStat Mat files to NIfTI images for viewing
% matname = (optional) filename to convert
%
%Example
% nii_mat2nii; %use GUI
% nii_mat2nii('vx.mat');

if ~exist('matname','var') %no files specified
   [A,Apth] = uigetfile({'*.mat;';'*.*'},'Select first image');
   matname = [Apth, A];
end;
if length(matname) < 1, return; end;
[kModalities, ~] = nii_modality_list();
mat = load(matname);
for m = 1 : size(kModalities,1)
    modality = deblank(kModalities(m,:));
    
    if isfield(mat, modality) && isfield(mat.(modality),'hdr')
        hdr = mat.(modality).hdr;
        img = mat.(modality).dat;
        [~,nam] = fileparts(matname);
        hdr.fname = [nam '_' modality '.nii'];
        spm_write_vol(hdr,img);
    end;
end