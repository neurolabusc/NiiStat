function nii_nii2matReport (niinames, ROIIndex, ThreshType)
%Convert NIfTI images to mat files for analysis with nii_stat_xls
% niinames = (optional) filenames to convert
% ROIIndex = number of region of interest
% ThreshType = -1=only negative, 0=non-zero, 1=only positive
%Example
% nii_nii2matReport
if ~exist('niinames','var') %no files specified
 niinames = spm_select(inf,'image','Select images to convert to mat');
end;
if length(niinames) < 1, return; end;
[kROI, kROINumbers] = nii_roi_list() ;
if ~exist('ThreshType','var') %get preferences
    prompt = {['RoiIndex (' sprintf('%s',kROINumbers) ')'],['Threhsold [-1=negativeVoxels,0=nonZeroVoxels,1=positiveVoxels']};
    dlg_title = 'Select region of interest';
    num_lines = 1;
    def = {'1','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer), return; end;
    ROIIndex = str2double(answer{1});
    ThreshType = str2double(answer{2});
end
%convert images
for i=1:size(niinames,1)
    niiname = deblank(niinames(i,:));
    niiname = unGzSub (niiname);
    [p,n,~] =spm_fileparts(niiname);
	hdr = spm_vol(niiname);
    img = spm_read_vols(hdr);
    if (ThreshType < 0) 
        img(img > 0) = 0;
        img(img < 0) = 1;        
    elseif (ThreshType > 0)
        img(img > 0) = 1;
        img(img < 0) = 0;
    else
        img(img ~= 0) = 1;
    end
    stat = [];
    ROIname = deblank(kROI(ROIIndex,:));
    [pth, proi] = fileparts(ROIname);
    new = nii_roi2stats(ROIname,hdr,img,'');
    vxls = nii_stats2roi (ROIname, hdr, img);
    fprintf('ThresholdedImage:\t%s\tROI:\t%s\tThresholdType\t%d\tNonZeroVoxelsImage\t%d\n',niiname,proi,ThreshType,sum(img(:) ~= 0));
    fprintf('Region\t%%ThresholdedImage\t%%ROI\tLabel\n');
    
    for i = 1:numel(new.(proi).mean)
        fprintf('%d\t%.3f\t%.3f\t%s\n',i,vxls(i)*100,new.(proi).mean(i)*100, new.(proi).label(i,:));
    end
end;
%end main function nii2matReport

function vxls = nii_stats2roi (roiName, inhdr, inimg);
rhdr = spm_vol ( [roiName '.nii']);
rimg = spm_read_vols (rhdr);
if ~isequal(inhdr(1).mat, rhdr.mat) || ~isequal(inhdr(1).dim(1:3), rhdr.dim(1:3))
    %we could warp the ROI to match the image, however since ROIs are indexed we would need to use nearest neighbor interpolation
    % therefore, we warp all images to match the ROI
    inimgOrig = inimg; %load input image 
    fprintf('reslicing %s to match %s\n',inhdr(1).fname,rhdr.fname);
    inimg = zeros([rhdr.dim(1:3),size(inimg,4)]);
    %inhdr1 = inhdr(1); %first volume
    imgdim = rhdr.dim;
    for v = 1: size(inimg,4) 
        inimgOrig1 = inimgOrig(:,:,:,v); 
        for i = 1:imgdim(3)
            M = inv(spm_matrix([0 0 -i])*inv(rhdr.mat)*inhdr(1).mat); %#ok<MINV>
            inimg(:,:,i,v) = spm_slice_vol(inimgOrig1, M, imgdim(1:2), 1); % (linear interp)
        end        
    end; %for each volume
end %if image must be resliced to match ROI
tot=sum(inimg(:));
mx = max(rimg(:));
vxls = zeros(mx,1);
for i=1:mx
    vxls(i) = sum(inimg(rimg==i))/tot;
    %fprintf('%d\t%d\n',i,sum(inimg(rimg==i)));
end %for i each roi
%end

function fnm = unGzSub (fnm)
[pth,nam,ext] = spm_fileparts(fnm);
if strcmpi(ext,'.gz') %.nii.gz
    fnm = char(gunzip(fnm));    
elseif strcmpi(ext,'.voi') %.voi -> 
    onam = char(gunzip(fnm));
    fnm = fullfile(pth, [nam '.nii']);
    movefile(onam,fnm);
end;  
%end unGzSub()