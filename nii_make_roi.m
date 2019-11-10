function nii_make_roi (roiNames)
%Given multiple binary images create single 'roi' image and text file
%  roiNames : region of interest to use (BA1.nii,BA2.nii...etc)
%Examples
% nii_make_roi %use GUI
% nii_make_roi(strvcat('ba6voidorsalcut.nii','BA45.nii','BA44.nii'))

if ~exist('roiNames','var') %roiNames not provided
    roiNames = spm_select(inf,'image','Select regions of interest');
end
if size(roiNames,1) > 255
    error('Please edit this script to create ROIs with more than 255 regions');
end
fid = fopen('roi.txt','w');
for i=1:size(roiNames,1)  
    roiName = deblank(roiNames(i,:));
    %read image data
    hdr = spm_vol (deblank (roiName));
    img = spm_read_vols (hdr);
    if i == 1
        hdr1 = hdr;
        imgROI = zeros(size(img));
    end
    if min(hdr1.dim ~= hdr.dim)
        fprintf('Warning: dimensions do not match: reslicing %s to match %s\n',hdr.fname,hdr1.fname);
        [hdr, img] = nii_reslice_target(hdr, img, hdr1, false);
    end
    img(img ~= max(img(:))) = nan;
    imgROI(~isnan(img)) = i; 
    %create text file
    [~,nam,~] = spm_fileparts(roiName);
    nam = regexprep(nam,' ','_');
    fprintf(fid,'%d|%s|%s\n', i, nam, nam);
end
fclose(fid); 
hdr1.fname = 'roi.nii';
hdr1.dt    =[2,0];%2=uint8
hdr1.pinfo = [1;0;0];%scale=1, intercept = 0
spm_write_vol(hdr1,imgROI);

%end nii_make_roi()
