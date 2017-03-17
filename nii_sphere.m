function nii_sphere(XYZmm, Radius, fnm)
%Create sphere with "Radius" voxel radius at location XYZmm using image fnm
% XYZmm : center of sphere, e.g. [0,0,0];
% Radius : size in voxels
% fnm : template image to use


if ~exist('XYZmm', 'var') || ~exist('Radius','var')
    prompt = {'Xmm','Ymm','Zmm','Radius (voxels)'};
    dlg_title = 'Sphere options';
    num_lines = 1;
    def = {'0','0','0','4'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    XYZmm = [0,0,0];
    XYZmm(1) = str2double(answer{1});
    XYZmm(2) = str2double(answer{2});
    XYZmm(3) = str2double(answer{3});
    Radius = str2double(answer{4});
end;
if ~exist('fnm', 'var') && exist('nii_roi_list','file')
    [files] = nii_roi_list;
    fnm = [deblank(files(1,:)), '.nii'];
elseif ~exist('fnm', 'var') 
    fnm = spm_select(1,'image','Select image to use as template'); 
end;
%create sphere
hdr = spm_vol(fnm);
hdr = hdr(1); %in case 4D
img = zeros(hdr.dim);
XYZvox = mm2voxSub(hdr, XYZmm);
img = addSphereSub (img,  Radius, XYZvox);
%save data
hdr.fname = sprintf('z%dy%dz%dd%d.nii', XYZvox(1), XYZvox(2), XYZvox(3), Radius+1+Radius);  
spm_write_vol(hdr,img);
fprintf('Created volume %s with %d voxels in %d voxel diameter sphere (same dimensions as %s)\n', hdr.fname, sum(img(:) > 0), Radius+1+Radius, fnm);
%end nii_sphere()

function img = addSphereSub (img,  Radius, XYZvox)
dx = zeros(size(img));
dy = zeros(size(img));
dz = zeros(size(img));

for x = 1 : size(img, 1)
   dx(x,:,:) = (x - XYZvox(1)).^2; 
end
for y = 1 : size(img, 2)
   dy(:,y,:) = (y - XYZvox(2)).^2; 
end
for z = 1 : size(img, 3)
   dz(:,:,z) = (z - XYZvox(3)).^2; 
end
dx = sqrt(dx + dy + dz);
img(dx <= Radius) = 1;
%end addSphereSub()

function xhair = mm2voxSub(hdr, XYZmm)
XYZmm = XYZmm(:);
mInv = inv(hdr.mat);
%xhair = mInv * [XYZmm; 1]; %convert from voxels to millimeters
%mInv = inv(hdr.mat);
xhair = hdr.mat \ [XYZmm; 1]; %convert from voxels to millimeters
xhair = round(xhair(1:3))';
xhair(xhair < 1) = 1;
xhair = min(xhair, hdr.dim);
%end mm2voxSub()




function findpeak (img,  Radius, XYZvox)
fprintf(' Sphere with radius of %f voxels\n',Radius); 
%draw a sphere centered at the peak
d = size(img);
lo = p-Radius;
hi = p+Radius;
lo((lo<1)) = 1;
if hi(1)>d(1) hi(1)=d(1); end;
if hi(2)>d(2) hi(2)=d(2); end;
if hi(3)>d(3) hi(3)=d(3); end;
img(:) = 0;    
    for z = lo(3) : hi(3) 
        for y = lo(2) : hi(2) 
                for x = lo(1) : hi(1)
                    %compute distance from origin using Pythagorean theorem
            dx = sqrt((p(1)-x)^2+(p(2)-y)^2+(p(3)-z)^2 );
                    if (dx<Radius) img(x,y,z)=1; end;
                end;%X
        end;%Y
end;%Z
%save the image
img(isnan(img)) = 0;
spm_write_vol(VO,img);