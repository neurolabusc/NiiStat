function nii_roi2mm (ROIIndex)
%report coordinates for each parcel in a region-of-interest map.
% ROIindex : (optional) number for region of interest
%Examples
% nii_roi2mm; %use GUI
% nii_roi2mm(2);

[kROI, kROINumbers] = nii_roi_list() ;
if ~exist('ROIIndex','var') %region of interest not specified
    ROIIndex = str2double(cell2mat(inputdlg(['RoiIndex (' sprintf('%s',kROINumbers) ')'], 'Choose image', 1,{'2'})));
end;
ROIname = deblank(kROI(ROIIndex,:));
label = labelSub ([ROIname '.txt']);
hdr = spm_vol ([ROIname '.nii']);
img = spm_read_vols (hdr);
fprintf('Region\tNumber\tX\tY\tZ\n');
for i = 1: max(img(:))
    %identify voxels in region
    img1 = zeros(size(img));
    img1(img == i) = 1;
    sumTotal = sum(img1(:));
    %compute center of mass
    vox(1) = sum(sum(sum(img1,3),2)'.*(1:size(img1,1)))/sumTotal; %dimension 1
    vox(2) = sum(sum(sum(img1,3),1).*(1:size(img1,2)))/sumTotal; %dimension 2
    vox(3) = sum(squeeze(sum(sum(img1,2),1))'.*(1:size(img1,3)))/sumTotal; %dimension 3
    %convert from voxels to mm
    mm=vox*hdr.mat(1:3, 1:3)'+hdr.mat(1:3, 4)';
    %report coordinates
    fprintf('%s\t%d\t%g\t%g\t%g\n',deblank(label(i,:)), i, mm(1), mm(2), mm(3)); 
end
%hdr.mat
%end nii_roi2mm()

function label = labelSub (ROIname)
if ~exist(ROIname,'file'), error('Unable to find %s',ROIname); end;
fid = fopen(ROIname);  % Open file
label=[];
tline = fgetl(fid);
while ischar(tline)
    %disp(tline)
    label=strvcat(label,tline); %#ok<REMFF1>
    tline = fgetl(fid);
end
fclose(fid); 
%end labelSub()