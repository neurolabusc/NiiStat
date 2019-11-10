function nii_mat2ortho(fnm, psnm)
%Display all images stored in a NiiStat matfile
% fnm : name of matfile to view
% psnm : (optional) append data to post-script file
%Examples
% nii_mat2ortho %use GUI
% nii_mat2ortho('M2020.mat')
% nii_mat2ortho('M2020.mat', 'my.ps')
% v = dir( '*.mat');
% v ={v.name}';
% for i = 1: size(v,1)
%   nii_mat2ortho(char(deblank(v(i,:))), 'out.ps');
% end
%Chris Rorden, 3/2016

if exist('fnm', 'var')
   [Apth, A, ext] = fileparts(fnm);
   A = {[A,ext]};
else
   [A,Apth] = uigetfile({'*.mat;';'*.*'},'MultiSelect', 'on','Select NiiStat-format Mat file');
else
    [A,Apth] = fileparts(file);
end
if ~exist('psnm', 'var') || isempty(psnm)
    psnm = fullfile(Apth,[mfilename date '.ps']);
end
if ~iscell(A)
    A = {A};
end
nfiles = numel(A);
for n = 1:nfiles
    fnm = fullfile(Apth, A{n});
    m = load(fnm);
    %count number of images
    f=fieldnames(m);
    f=sort(f);
    nImgs = 0;
    for i = 1: numel(f)
        if isfield( m.(f{i}),'dat') && isfield( m.(f{i}),'hdr')
            nImgs = nImgs + 1;
        end
    end
    if nImgs < 1
        fprintf('No images found in %s\n', fnm);
        return;
    end
    %set crosshairs lesion to center of mass
    XYZmm = [0;0;0];
    if isfield(m,'lesion') && isfield(m.lesion,'dat')
        XYZmm = getCenterOfIntensitySub(m.lesion.hdr, m.lesion.dat);
    end
    %create a new figure
    [~,fnm] = fileparts(fnm); % /home/cr/m2020.mat -> 'm2020'
    clf;
    nImg = 0;
    for i = 1: numel(f)
        if isfield( m.(f{i}),'dat') %&& sfield( m.(f{i}),'hdr')
            nImg = nImg + 1;
            str = sprintf('%s\n%s', f{i}, fnm);
            plotOrthoSub(m.(f{i}).hdr, m.(f{i}).dat, XYZmm, nImg, nImgs, str);
        end
    end
    if ~exist('psnm', 'var') || isempty(psnm), return; end;
    print('-dpsc', '-append', psnm);
    clear m;
end
%end nii_mat2ortho








function plotOrthoSub(hdr, img, XYZmm, Slot, NumSlots, Caption)
xhair = mm2voxSub(hdr, XYZmm);
set(gcf,'color','w');
%img = img - min(img(:)); %set minimum to zero
img(img < 0) = 0;
thresh = max(img(:));
if  numel(img) > 1000 
    imgS = sort(img(:));
    pct = round(numel(imgS) * 0.995);
    pct = imgS(pct);
    if pct > 0
       thresh = pct; 
    end
end

%img(img > thresh) = thresh;
img = 63 * (img / thresh); %matlab color scheme have 64 indices "size(bone)"
sz = size(img);
colormap(gray); %colormap(bone)
ax = img(:,:,xhair(3));
ax(xhair(1),:) = 128;
ax(:,xhair(2)) = 128;
cor = squeeze(img(:,xhair(2),:));
cor(xhair(1),:) = 128;
cor(:,xhair(3)) = 128;
sag = squeeze(img(xhair(1),:,:));
sag(xhair(2),:) = 128;
sag(:,xhair(3)) = 128;
im = zeros( sz(3)+sz(2), sz(1)+sz(2));
im(1:sz(1),1:sz(2)) = ax;
im(1:sz(1),(1:sz(3))+sz(2)) = cor;
im((1:sz(2))+sz(1),(1:sz(3))+sz(2)) = sag;
%scale output
rows = ceil(sqrt(NumSlots)); %e.g. 2..4 items shown in 2x2 mosaic, 5..9 in 3x3
scale = 1/rows;
col = mod(Slot-1,rows) * scale;
row = floor((Slot-1)/rows) * scale;
plotImgSub ( flipud(im'), Caption, col, row, scale, scale);
%plotOrthoSub

function xhair = mm2voxSub(hdr, XYZmm)
mInv = inv(hdr.mat);
xhair = mInv * [XYZmm; 1]; %convert from voxels to millimeters
xhair = round(xhair(1:3))';
xhair(xhair < 1) = 1;
xhair = min(xhair, hdr.dim);
%end mm2voxSub

function plotImgSub ( Img, Caption, X, Y, wid, ht)
subplot('Position',[X Y wid ht]); %width, height
image((Img));
set(gca,'XTickLabel', [],'XTick',[],'YTick',[]);
axis image
w = wid/2;
h = ht/2;
annotation('textbox',[X+w Y w h],'String',Caption,'FontSize',18,'fontn','Arial', 'color','red', 'BackgroundColor', 'white', 'LineStyle','none')
%end plotImgSub()

function XYZmm = getCenterOfIntensitySub(hdr, img)
XYZmm = ones(3,1);
img = img - min(img(:));
img(isnan(img)) = 0;
%find center of mass in each dimension (total mass divided by weighted location of mass
% img = [1 2 1; 3 4 3];
sumTotal = sum(img(:));
coivox = ones(4,1);
coivox(1) = sum(sum(sum(img,3),2)'.*(1:size(img,1)))/sumTotal; %dimension 1
coivox(2) = sum(sum(sum(img,3),1).*(1:size(img,2)))/sumTotal; %dimension 2
coivox(3) = sum(squeeze(sum(sum(img,2),1))'.*(1:size(img,3)))/sumTotal; %dimension 3
XYZmm = hdr.mat * coivox; %convert from voxels to millimeters
XYZmm = XYZmm(1:3);
%end setCenterOfIntensitySub()
