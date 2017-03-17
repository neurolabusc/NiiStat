function nii_powermap (fname, nSubj)
%Generate power map image
% fname : name of sum map (32-bit real)
% nSubj : number of participants used to make sum image
%Chris Rorden 1/2015
% for details see http://www.ncbi.nlm.nih.gov/pubmed/19285465/
% http://opensource.org/licenses/BSD-2-Clause
%Examples
%   nii_powermap; %GUI
%   nii_powermap('41sum.nii', 41);

warning('Obsolete: please use nii_power');
if ~exist('fname','var') %no image
    fname = spm_select(1,'image','Select sum image');
    [pth,nam,ext] = spm_fileparts( fname);
    fname = fullfile(pth,[ nam, ext]); %'img.nii,1' -> 'img.nii'
end
if ~exist('nSubj','var')
	nSubj = cell2mat(inputdlg('Enter sample size (theoretical maximum for sum map):', 'Details', 1,{'41'}));
    nSubj = str2double(nSubj);
end
fprintf('Warning: assumes continuous behavioral data (binomial data less extreme\n');
[pth,nam,ext] = spm_fileparts( fname);
hdr = spm_vol (deblank (fname));
if (hdr.dt ~= 16), error('Currently only supports 32-bit real input images'); end;
img = spm_read_vols (hdr);
fprintf('Image voxel range %g..%g\n', min(img(:)), max(img(:)));
mx = max(img(:));
if (mx > nSubj) || (min(img(:)) < 0) , error('Values in sum image must be in range 0..%d\n', nSubj); end;
if ~(abs(round(mx)-mx)) <= eps('double')
    error('Image should be integers (number of people with injury at location)');
end
for v=1 : numel(img)
    if (img(v) > 0)
        img(v) = spm_invNcdf( 1/nchoosek(nSubj, img(v)) );
    end
end
img = abs(img);
hdr.fname = fullfile(pth, ['powermap_' nam ext]);
spm_write_vol(hdr,img);