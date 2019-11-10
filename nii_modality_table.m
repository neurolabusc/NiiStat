function nii_modality_table(pth)
%Create a table that shows modalities available for NiiStat
% pth: (optional) folder to search for Mat files, else current directory
%Examples
% nii_modality_table; %current workng directory
% nii_modality_table('~/MatFiles');

if ~exist('pth','var')
    pth = pwd;
end
if ~isdir(pth), error('Invalid folder %s', pth); end;
%1st row: table labels
modality = nii_modality_list;
outtext = fullfile(pth,'NiiStat_modalities.txt');
if exist(outtext,'file'), delete(outtext); end;
diary(outtext);
fprintf('ID\t');
for i = 1: size(modality,1)
    fprintf('%s\t', deblank(modality(i,:)))
end
fprintf('\n');
%subsequent rows: individuals
fnames = dir(fullfile(pth,'*.mat'));
if isempty(fnames), return; end;
for f = 1 : numel(fnames)
    nam = deblank(fnames(f).name);
    mat = load(fullfile(pth,nam));
    if ~isfield(mat,'T1') || ~isfield(mat.T1,'hdr') continue; end;
    fprintf('%s\t', nam)
    for i = 1: size(modality,1)
        fld = isfield(mat, deblank(modality(i,:)));
        if ~fld
            fname = [deblank(modality(i,:)), '_jhu'];
            fld = isfield(mat, fname);
        end
        fprintf('%d\t', fld)
    end %for each modality
    fprintf('\n');
end %for each mat file
diary OFF