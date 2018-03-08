function [stat] = nii_roi2stats (roiName, inhdr, inimg, statName,diskName)
%compute mean intensity of input image for every region of interest
% NB: replaces nii_roi_stats: now able to reslice data if required
%  roiName : region of interest to use ('jhu', 'bro', 'fox', 'tpm','aal')
%  inhdr : either name of NIfTI image or NIfTI header strucutre
%  inimg : (only used if inhdr is a structure): input image data
%  statName : (optional) prefix for struct, e.g. 'rest', 'cbf'
%  diskName : if not empty, data appended to mat file on disk
%Examples
% s = nii_roi2stats('jhu','LS_LM1001.nii','','lesion_') %3d volume
% s = nii_roi2stats('jhu1mm','REST_LM1001.nii,'','lesion_') %4d volumes

if ~exist('roiName','var')
    roiName = 'fox';
end
if ~exist('inhdr','var') %no image
 inhdr = spm_select(1,'image','Select voxelwise NIfTI image');
end;
if ~exist('diskName','var')
    diskName = '';
end;
if ~exist('statName','var')
    statName = '';
end;
if (length(inhdr) < 1), return; end; %escape if no selection
if ~isstruct(inhdr)
    if (exist(inhdr,'file') == 0)
        fprintf('%s error: unable to find image named %s\n',mfilename,inhdr);
        return;
    end;
    inhdr = spm_vol(inhdr); %load input header
	inimg = spm_read_vols(inhdr); %load input image
end
if (length(inimg) < 1),
    fprintf('No image data to reslice');
    return;
end; %escape if no selection
%niiName = [roiName '.nii'];
niiName = [fileparts(which(mfilename))  filesep 'roi' filesep roiName '.nii']; %precedence expected NiiStat location
if exist(niiName, 'file') ~= 2 %unable to find image
    %niiName = [fileparts(which(mfilename))  filesep 'roi' filesep roiName '.nii'];
    niiName = [roiName '.nii']; %fallback to somewhere in path
    if exist(niiName, 'file') ~= 2
        fprintf('%s did not find region of interest named %s\n',mfilename, niiName);
        return;
    end
end
%niiName = which(niiName); %provide full path

[pth, roiShortName] = fileparts(niiName);
if isempty(pth)
    niiName = which(niiName);
    [pth, roiShortName] = fileparts(niiName);
end

txtName = fullfile(pth, [roiShortName '.txt']);
if exist(txtName, 'file') ~= 2 %unable to find image
    fprintf('%s did not find region of interest named %s\n',mfilename, txtName);
    return;
end
proi = [statName  roiShortName];
% 666 todo make sure .label has same number of elements as .mean

stat.(proi).label = labelSub (txtName);
rhdr = spm_vol (deblank (niiName));
rimg = spm_read_vols (rhdr);
if ~isequal(inhdr(1).mat, rhdr.mat) || ~isequal(inhdr(1).dim(1:3), rhdr.dim(1:3))
    if ndims(inimg) < 4
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

%         h =rhdr(1);
%         h.fname = 'test.nii';
%         h.dt(1) = 16; %float
%         spm_write_vol(h,inimg(:,:,:,1));
    else
        %for 4D images, we reslice the ROI to match the imag
        rimgOrig = rimg; %load input image
        fprintf('reslicing %s to match %s\n',rhdr.fname, inhdr(1).fname);
        rimg = zeros(inhdr(1).dim(1:3));
        %inhdr1 = inhdr(1); %first volume
        imgdim = inhdr(1).dim(1:3);
        for i = 1:imgdim(3)
            M = inv(spm_matrix([0 0 -i])*inv(inhdr(1).mat)*rhdr.mat); %#ok<MINV>
            rimg(:,:,i) = spm_slice_vol(rimgOrig, M, imgdim(1:2), 0); % (nearest neighbor interp)
        end
%         h =inhdr(1);
%         h.fname = 'test.nii';
%         h.dim = inhdr(1).dim(1:3);
%         spm_write_vol(h,rimg);
    end
end %if image must be resliced to match ROI
if size(inimg,4) == 1 %3D data
    stat.(proi).mean = roiMeanSub (inimg, rimg);
    if (size(stat.(proi).label,1) ~= size(stat.(proi).mean,1))
        error('Error: somthing is wrong with the ROI file: %d labels but %d regions', size(stat.(proi).label,1), size(stat.(proi).mean,1));
    end;
    
else
    [stat.(proi).r, stat.(proi).p] = roiCorrelSub (inimg, rimg);
    if (size(stat.(proi).label,1) ~= size(stat.(proi).r,1))
        error('Error: somthing is wrong with the ROI file: %d labels but %d regions', size(stat.(proi).label,1), size(stat.(proi).r,1));
    end;
end
if length(diskName) < 1, return; end
if exist(diskName,'file')
    old = load(diskName);
    stat = nii_mergestruct(stat,old);
end
save(diskName, '-struct', 'stat');
%end nii_roi2stats()

function label = labelSub (txtName)
fid = fopen(txtName);  % Open file
label=[];
tline = fgetl(fid);
while ischar(tline)
    %disp(tline)
    label=strvcat(label,tline); %#ok<REMFF1>
    tline = fgetl(fid);
end
fclose(fid);
%end labelSub()

function [mn] = roiMeanSub (img, rimg)
%find mean intensity for each region of interest in a parcellated image
if ~isequal(size(img(1:3)),size(rimg(1:3)))
   fprintf('roiMeanSub error: images must have same dimensions\n');
   return;
end
nroi = max(rimg(:));
%save('img.mat','-struct','img');
%save('rimg.mat','-struct','rimg');
%rimg(isnan(img)) = 0; %<- a simpler way to handle NaNs?
mn = zeros(nroi,1);
for r = 1:nroi
    %idx = (rimg(:) == r);
    idx = (rimg(:) == r) & (~isnan(img(:)) );
    mn(r) = mean(img(idx));
    %fprintf('Region %d has %d voxels with a mean of %f\n',r,sum(idx), mn(r));
end

%end roiMeanSub

function [r,p] = roiCorrelSub (img, rimg, showGraph)
%extract correlation coefficients from a 4D NIfTI image using a parcellated 3D region of interest image
if ~exist('showGraph','var') %no image
    showGraph = false;
end
if ~isequal(size(img(1:3)),size(rimg(1:3)))
   fprintf('roiCorrelSub error: image and ROI must have same dimensions\n');
   return;
end
[~, ~, ~, nV] = size(img);
if (nV < 2)
   fprintf('%s requires 4D images, unable to use %s\n',mfilename,imgName);
   return;
end
nroi = max(rimg(:));
dat = zeros(nV,nroi); %pre-allocate
tmp = reshape(img, [], nV);
img1 = img(:,:,:,1);
rimg(isnan(img1)) = 0; %remove voxels outside brain mask
for r = 1:nroi
    idx = (rimg(:) == r);
    dat(:,r) = mean(tmp(idx,:));
end
% for r = 1:nroi
%     for v = 1: nV
%         img1 =  img(:,:,:,v);
%         idx = (rimg(:) == r) & (~isnan(img1(:)) );
%         if isempty(idx)
%             fprintf(' ROI does not have index %d\n',r);
%         end
%         dat(v,r) = mean(img1(idx));
%     end
% end

[r,p] = corrcoef(dat);
if sum(~isfinite(p(:))) > 0
    fprintf('%s generated some not-a-number correlations. Perhaps numerical gaps in region of interest file\n',mfilename);
    r(~isfinite(r))=0;
    p(~isfinite(p))=1;
end
if ~showGraph, return; end;
imagesc(r);
colormap('jet');colorbar;
%end roiCorrelSub
