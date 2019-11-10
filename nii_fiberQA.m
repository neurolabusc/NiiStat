function nii_fiberQA(baseDir)
%Given folder with many NiiStat Mat files, open all with DTI data and report which
% have unusual signal in the right hemisphere

if exist('baseDir','var')
    cd(baseDir);
end
isRightHemisphereSpared = true; %set to false for RHD and true for LHD
nii_correlSub('fa_jhu', isRightHemisphereSpared, true, true); %WHITE MATTER!!!
nii_correlSub('fmri_AICHA', isRightHemisphereSpared, true);
nii_correlSub('cbf_AICHA', isRightHemisphereSpared, true);
nii_correlSub('md_AICHA', isRightHemisphereSpared,  true);
nii_correlSub('mk_AICHA', isRightHemisphereSpared, true);
nii_correlSub('palf_AICHA', isRightHemisphereSpared);
nii_correlSub('alf_AICHA', isRightHemisphereSpared);
nii_correlSub('i3mT1_AICHA', isRightHemisphereSpared, true);
nii_correlSub('dtimn_jhu', isRightHemisphereSpared);
nii_correlSub('dtimn_AICHA', isRightHemisphereSpared);
nii_correlSub('rest_AICHA', isRightHemisphereSpared);
%these next to are slow, comment the following line for faster speed
nii_correlSub('RestAve_jhu',isRightHemisphereSpared, true);
nii_correlSub('fMRIave_jhu',isRightHemisphereSpared, true);
%end nii_modality_table_qa()

function nii_correlSub(fieldname, isRightHemisphereSpared, isNorm, isWhiteMatter)
if ~exist('isNorm','var')
    isNorm = false; %assume use does not want to scale intensities from 0..1
end
if ~exist('isWhiteMatter','var')
    isWhiteMatter = false; %assume use does not want to scale intensities from 0..1
end
m = dir('*.mat');
if isempty(m), error('Unable to find mat files'); end;
fprintf('Found %d subjects in %s\n',numel(m), pwd);
label = labelSub (fieldname);
isGM_Hemi = isHemiSub(fieldname, label, isRightHemisphereSpared, isWhiteMatter);
nROI = sum(isGM_Hemi(:));
msk = triu(ones(nROI,nROI),1);
mat = [];
nam = [];
matOK = 0;
tic
for s = 1:numel(m)
    snam = m(s).name;
    mt = getMatSub(snam, fieldname, msk, isGM_Hemi);
    if isempty(mt)
        %fprintf('Skipping %s\n', snam); %optional: report individuals with missing modality
        continue;
    end
    if isNorm
       mt = normSub(mt); %normalize intensity for range of 0..1 
    end
    matOK = matOK + 1;
    mat(matOK,:) = mt;
    nam = [nam, {snam}]; %#ok<*AGROW>
end %for each subject - generate population statistics
if matOK == 0, return; end;
%fprintf('Found %d files with field %s\n',  matOK, fieldname);
if matOK < 2, fprintf('Not enough participants to compute "%s" QA statistics\n', fieldname); return; end; %stdev requires multiple samples
mn = mean(mat);
rOK = zeros(matOK,1);
fid = fopen([fieldname '.txt'], 'w');
for s = 1:matOK
    [r] = corrcoef (mn, mat(s,:));
    rOK(s) = r(2,1);
end
fclose(fid);
%report individual values relative to Group
mnG = mean(rOK);
stdG = std(rOK);
[~, idx] = sort(rOK); %list individuals sorts by QA
for s = 1:matOK
    %subj = s; %for input order
    subj = idx(s); %for sorted order
    z = (rOK(subj) - mnG) / stdG;
    fprintf('%s\tr=\t%g\tz=\t%g\n',  nam{subj}, rOK(subj), z );
end
fprintf('For %d files with %s, r min=\t%g\tmean=\t%g\tmax=\t%g\tstdev=\t%g\n',  matOK, fieldname, min(rOK), mnG, max(rOK), stdG);
toc
%end fiberQA()

function img = normSub(img);
finite = isfinite(img);
imgFinite = img(finite);
if isempty(imgFinite), return; end;
mn = min(imgFinite);
rng = max(imgFinite) - mn;
if rng == 0, return; end; %no range
img = (img - mn)/rng; 
img(~finite) = 0;
%end normSub()

function mt = getMatSub(matname, fieldname, msk, isGM_R)
mt = [];
fld = fieldSub(matname, fieldname);
if  isempty(fld) %compute ROI e.g. fieldname=RestAve_jhu with m.RestAve.hdr
    [~,nam,roi] = fileparts(strrep(fieldname, '_', '.')); %dtimn_jhu -> .jhu
    roi = strrep(roi, '.', ''); %.jhu -> jhu
    fld = fieldSub(matname, nam);
    if  isempty(fld), return; end; %e.g. RestAve does not exist
    if ~isfield(fld,'hdr'), return; end; %e.g. RestAve.hdr does not exist
    s = nii_roi2stats (roi, fld.hdr, fld.dat);
    fld = s.(roi);
    if ~isfield(fld,'mean')
        error('Unable to load "%s" from the file %s', fieldname, matname);
    end
    mt =  fld.mean; %189 row vector
    mt = mt(isGM_R);
    %mt = mt(msk(1,:) ~=0);
    mt(~isfinite(mt)) = 0;
    return; 
end;

if isfield(fld,'r')
    mt =  fld.r; %189x189 matrix
    mt = mt(isGM_R,isGM_R); %53x53 matrix
    mt = mt(msk ~=0);
else
    if ~isfield(fld,'mean')
        error('Unable to load "%s" from the file %s', fieldname, matname);
    end
    mt =  fld.mean; %189 row vector
    mt = mt(isGM_R);
    %mt = mt(msk(1,:) ~=0);
    mt(~isfinite(mt)) = 0;
end
%end


function label = labelSub(fieldname)
[~,~,ext] = fileparts(strrep(fieldname, '_', '.')); %dtimn_jhu -> .jhu
ext = strrep(ext, '.', ''); %.jhu -> jhu
pth = which('NiiStat');
[pth] = fileparts (pth);
pth = fullfile(pth , 'roi', [ext '.txt']);
if ~exist(pth,'file'), error('Unable to find %s\n',pth); end;
fid = fopen(pth);  % Open file
label=[];
tline = fgetl(fid);
while ischar(tline)
    label=strvcat(label,tline); %#ok<DSTRVCT,REMFF1>
    tline = fgetl(fid);
end
fclose(fid);
%end labelSub()

function fld = fieldSub(matname, fieldname)
m = load(matname);
if isfield(m, fieldname)
    fld = m.(fieldname);
else
    fld = [];
end
%end fieldSub()

function isGM_Hemi = isHemiSub(fieldname, label, isRightHemi, isWhiteMatter)
n = size(label,1);
isGM_Hemi = zeros(n,1);
if isRightHemi
    HemiKey = 'R|';
else
    HemiKey = 'L|';
end
for i = 1: n
    if isWhiteMatter
        if ~isempty(strfind(label(i,:), HemiKey)) && ~isempty(strfind(label(i,:), '|2')) %right hemisphere,  white matter
            isGM_Hemi(i) = 1;
            %fprintf('%s\n',label(i,:) );
        end        
    else
        if ~isempty(strfind(label(i,:), HemiKey)) && isempty(strfind(label(i,:), '|2')) && isempty(strfind(label(i,:), '|3')) %right hemisphere, not CSF or white matter
            isGM_Hemi(i) = 1;
            %fprintf('%s\n',label(i,:) );
        end
    end
end
fprintf('%s has %d regions of type "%s"\n', fieldname, sum(isGM_Hemi(:)), HemiKey);
isGM_Hemi = logical(isGM_Hemi);
%end isGMr_Sub()
