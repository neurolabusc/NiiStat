function nii_array2roi (array, roifname, outname)
%Saves an array to an index region of interest
% array : e.g. [2 3 3 5]; means region1=2, region2=3, region3=3, etc
% roifname : name of region of interext file
%Examples
% nii_array2roi([11 12 13 14 15 16 17 18 19 20], 'fox.nii')

if ~exist('array','var') || ~exist('roifname','var')
    return
end
if exist(roifname,'file') ~= 2
    fprintf('No statistcal map saved: unable to find image named %s\n',roifname);
    return
end
if nargin < 3
    [pth, nam, ext] = fileparts(roifname);    
    outname = fullfile(pth, ['array2', nam, ext]);
end
hdr = spm_vol (roifname);
img = spm_read_vols(hdr); %load ROI
if (max(img(:)) ~= numel(array))
    fprintf('ROI %s has %d regions, whereas the array has %d\n',roifname, max(img(:)), numel(array));
    return;
end
%create new headers with same dimensions as ROI
hdrC = hdr;
[pth, nam] = fileparts(outname);
hdrC.fname = fullfile(pth, [nam, '.nii']);
hdrC.dt(1) = 16; %make 32 bit real
hdrC.private.dat.dtype = 'FLOAT32-LE';
hdrC.private.dat.scl_slope = 1;
hdrC.private.dat.scl_inter = 0;
imgC = zeros(size(img));
for i = 1:numel(array)
    imgC(img == i) = array(i);
end
spm_write_vol(hdrC,imgC);
%end nii_array2roi()