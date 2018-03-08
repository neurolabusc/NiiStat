function [stat, roi_names] = nii_nii2roi (roiNames, inhdr, inimg)
%compute mean intensity of input image for every region of interest
%  roiName : region of interest to use ('jhu', 'bro', 'fox') 
%  inhdr : either name of NIfTI image or NIfTI header strucutre
%  inimg : (only used if inhdr is a structure): input image data
%Examples
% [s, r] = nii_nii2roi('b1.nii','alfREST_LM1035.nii') %single ROI
% [s, r] = nii_nii2roi(strvcat('b1.nii','b2.nii'),'alfREST_LM1035.nii','') %2 ROIs

normalizeRoi = true; %sets ROI range from 0..1, disable for probability maps where peak should be less than 1.0
if ~exist('roiNames','var')
    roiNames = spm_select(inf,'image','Select regions of interest');
end
if ~exist('inhdr','var') %no image
 inhdr = spm_select(1,'image','Select voxelwise NIfTI image');
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
if size(inimg,4) > 1
        fprintf('%s not designed for 4D images!\n',mfilename);
end
if ~isequal(inhdr.dim,size(inimg))
    inimg = reshape(inimg,inhdr.dim);
end
%stat.label = [];
roi_names =[];
for j=1:size(roiNames,1)  
    roiName = deblank(roiNames(j,:));
    [pth,nam,ext] = spm_fileparts(roiName);
    if ~exist(roiName,'file') 
        roiName = fullfile(pth,[nam ext]); %remove volume, e.g. 'img.nii,1' -> 'img.nii'
    end
    if ~exist(roiName,'file') 
        roiName = [nam ext]; %remove volume and path (check cwd)
    end
    if ~exist(roiName,'file')
        [mpth,~,~] = spm_fileparts( deblank (which(mfilename)));
        roiName = fullfile(mpth,[nam '.nii']);
        if ~exist(roiName,'file')
            error('Unable to find file %s',roiName);
        end
    end
    rhdr = spm_vol (deblank (roiName)); %load ROI header
    rimg = spm_read_vols (rhdr); %load ROI image
    if normalizeRoi
        if (min(rimg(:)) ~= 0) || (max(rimg(:)) ~= 0) %only if required
            fprintf('Scaling %s to be in the range 0..1\n', roiName);
            rimg = rimg - min(rimg(:)); %translate so range 0..max
            if ~isfinite(max(rimg(:))) || (max(rimg(:)) == 0)
                error('Not a valid region of interest %s',roiName);
            end
            rimg = rimg ./ max(rimg(:)); %scale so range now 0..1
        end
    end
    if sum(~isfinite(rimg(:))) > 0 
        fprintf('WARNING: not all voxels are finite - these are set to zero %s\n', rhdr.fname);
        rimg( ~isfinite(rimg)) = 0;
    end
    if ~isequal(inhdr.mat, rhdr.mat) || ~isequal(inhdr.dim(1:3), rhdr.dim(1:3))
        %we could warp the ROI to match the image, however since ROIs are indexed we would need to use nearest neighbor interpolation
        % therefore, we warp all images to match the ROI
        fprintf('Reslicing image to match %s\n',rhdr.fname);
        imgdim = rhdr.dim;
        inimgR = rimg; %preallocate
        for i = 1:imgdim(3)
            M = inv(spm_matrix([0 0 -i])*inv(rhdr.mat)*inhdr(1).mat); %#ok<MINV>
            inimgR(:,:,i) = spm_slice_vol(inimg, M, imgdim(1:2), 1); % (linear interp)         
        end; %for each slice
    else %if reslice else image matches dimension of ROI
        inimgR = inimg;
    end %if image must be resliced to match ROI
    if sum(~isfinite(inimgR(:))) > 0 
        fprintf('WARNING: not all voxels are finite - these are set to zero %s\n', inhdr.fname);
        inimgR( ~isfinite(inimgR)) = 0;
    end
    [~, roiNam, ~] = fileparts(roiName); 
    %stat.label = strvcat(stat.(statName).label,roiName);
    roi_names = strvcat(roi_names, roiNam); %#ok<REMFF1>
    %roi_names{j} = roiNam;
    idx = (~isnan(rimg(:))) & (~isnan(inimgR(:)) ); %only include numeric voxels
    %stat.mean(j) = mean(inimgR(idx).*rimg(idx)); %modulate inimg with ROI <- mean for volume
    stat(j) = sum(inimgR(idx).*rimg(idx)) / (sum(rimg(idx))); %#ok<AGROW> %modulate inimg with ROI
    %stat(j) = mean(inimgR(idx).*rimg(idx)); %#ok<AGROW> %modulate inimg with ROI
    %Test routines save output to disk    
    %rhdr.fname = fullfile(pth,['r' roiNam '.nii']);
    %spm_write_vol(rhdr,inimgR); %save resized image
    % inhdr.fname = fullfile(pth,['i' roiNam '.nii']);
    % inimgR = reshape(inimg,inhdr.dim(1:3));
    % spm_write_vol(inhdr,inimgR); %save resized image
end
%end nii_nii2roi()
