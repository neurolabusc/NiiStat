function nii_array2node (nodeArray, edgeArray, roifname, outname)
%Saves an array as BrainNet Viewer node file
% array : e.g. [2 3 3 5]; means region1=2, region2=3, region3=3, etc
% roifname : name of region of interext file
%Examples
% nii_array2node([1 2 3 4 nan 6 7 8 9 10],50*rand(10,10),'/Users/rorden/NiiStat/roi/fox.nii','theNodes.node')
% nii_array2node([1 2 3 4 nan 6 7 8 9 10],[],'/Users/rorden/NiiStat/roi/fox.nii')
%Example that displays results
%  nam = which('nii_stat');
%  pth = fileparts(nam);
%  nam = fullfile(pth, [filesep 'roi' filesep 'fox.nii']);
%  nii_array2node([1 2 3 4 5 6 7 8 9 10],50*rand(10,10),nam,'theNodes.node')
%  MRIcroS('closeAllItems'); 
%  MRIcroS('addLayer','BrainMesh_ICBM152.nv'); %open image
%  MRIcroS('addNodes','theNodes.node','theNodes.edge'); %add hot spots  

if ~exist('roifname','var') %region of interest not specified
    [kROI, kROINumbers] = nii_roi_list() ;
    ROIIndex = str2double(cell2mat(inputdlg(['RoiIndex (' sprintf('%s',kROINumbers) ')'], 'Choose image', 1,{'2'})));
    roifname = [deblank(kROI(ROIIndex,:)) '.nii'];
end;
if ~exist('nodeArray','var') || ~exist('roifname','var')
    fprintf('Invalid inputs\n');
    return
end
if exist(roifname,'file') ~= 2
    fprintf('No statistcal map saved: unable to find image named %s\n',roifname);
    return
end
[pth, nam, ext] = fileparts(roifname);
label = labelSub (fullfile(pth, [nam, '.txt']));
if  ~exist('outname','var') || isempty(outname)
    outname = fullfile(['array2', nam, '.node']);
end
hdr = spm_vol (roifname);
img = spm_read_vols(hdr); %load ROI
if (max(img(:)) ~= numel(nodeArray))
    fprintf('ROI %s has %d regions, whereas the array has %d\n',roifname, max(img(:)), numel(nodeArray));
    return;
end
%create new headers with same dimensions as ROI
fid = fopen(outname,'w');
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
    if isfinite(nodeArray(i))
        l = deblank(label(i,:));
        l = regexprep(l,' ','_');
        fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%s\n', mm(1), mm(2), mm(3),1,nodeArray(i),l);
    end
end
fclose(fid);
if  ~exist('edgeArray','var') || ( numel(edgeArray) ~= numel(nodeArray)^2) 
    exit;
end;
[pth, nam, ~] = fileparts(outname);
outname = fullfile(pth, [ nam, '.edge']);
save(outname,'edgeArray','-ascii')
%end nii_array2node()

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
