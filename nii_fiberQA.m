function nii_fiberQA(baseDir)
%Given folder with many NiiStat Mat files, open all with DTI data and report which
% have unusual signal in the right hemisphere

if exist('baseDir','var')
    cd(baseDir);
end
nii_correlSub('dtimn_jhu');
nii_correlSub('dtimn_AICHA')
%nii_modality_table_qa



function nii_correlSub(fieldname)
m = dir('*.mat');
if isempty(m), error('Unable to find mat files'); end;
fprintf('Found %d subjects in %s\n',numel(m), pwd);
label = labelSub (fieldname);
isGM_Hemi = isHemiSub(fieldname, label, true);
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
        fprintf('Skipping %s\n', snam);
        continue;
    end
    matOK = matOK + 1;
    mat(matOK,:) = mt;
    nam = [nam, {snam}]; %#ok<*AGROW>
end %for each subject - generate population statistics
if matOK == 0, return; end;
%fprintf('Found %d files with field %s\n',  matOK, fieldname);
mn = mean(mat);
rOK = zeros(matOK,1);
fid = fopen([fieldname '.txt'], 'w');
for s = 1:matOK
    [r] = corrcoef (mn, mat(s,:));
    fprintf('%s\tr=\t%g\n',  nam{s}, r(2,1) );
    fprintf(fid, '%s\tr=\t%g\n',  nam{s}, r(2,1) );
    rOK(s) = r(2,1);
end
fclose(fid);
fprintf('For %d files with %s, r min=\t%g\tmean=\t%g\tmax=\t%g\n',  matOK, fieldname, min(rOK), mean(rOK), max(rOK));
toc
%end fiberQA()


function mt = getMatSub(matname, fieldname, msk, isGM_R)
fld = fieldSub(matname, fieldname);
mt = [];
if  isempty(fld), return; end;
%fprintf('Processing %s\n',  matName);
mt =  fld.r; %189x189 matrix
mt = mt(isGM_R,isGM_R); %53x53 matrix
mt = mt(msk ~=0);
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

function isGM_Hemi = isHemiSub(fieldname, label, isRightHemi)
n = size(label,1);
isGM_Hemi = zeros(n,1);
if isRightHemi
    HemiKey = 'R|';
else
    HemiKey = 'L|';
end
for i = 1: n
    if ~isempty(strfind(label(i,:), HemiKey)) && isempty(strfind(label(i,:), '|2')) && isempty(strfind(label(i,:), '|3')) %right hemisphere, not gray or white matter
        isGM_Hemi(i) = 1;
        %fprintf('%s\n',label(i,:) );
    end
end
fprintf('%s has %d regions of type "%s"\n', fieldname, sum(isGM_Hemi(:)), HemiKey);
isGM_Hemi = logical(isGM_Hemi);
%end isGMr_Sub()