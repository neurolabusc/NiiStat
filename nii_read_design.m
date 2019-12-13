function [designMat, designUsesNiiImages, CritN, nuisanceMat] = nii_read_design (xlsname, worksheetname)
designUsesNiiImages = false;
CritN = [];
if exist(xlsname,'file') ~= 2
    fprintf('%s Unable to find file named "%s"\n',mfilename, xlsname);
    return
end
if ~exist('worksheetname','var')
    worksheetname = 'NiiStat';
end
[~,~,x] = fileparts(xlsname);
if ~exist(xlsname,'file')
   error('Unable to find file %s\n', xlsname);
end
if strcmpi(x,'.tab') || strcmpi(x,'.txt')  || strcmpi(x,'.val')
    [dMat, CritN] = nii_tab2mat(xlsname);
else
    [dMat, nuisance] = nii_xls2mat(xlsname , worksheetname,'', true);
    if isempty(dMat)
        [dMat, nuisance] = nii_xls2mat(xlsname , 'Data (2)','', true);
    end
end
if isempty(dMat)
   error('Unable to load file %s (maybe not Matlab compatible Excel format)\n', xlsname);
end
SNames = fieldnames(dMat);
numFields = length (SNames);
if numFields < 2
    error('File %s must have multiple columns (a column of file names plus a column for each behavior\n', xlsname);
end
% imgpath = nii_update_mat (fileparts(which(xlsname))); %use which for full path
imgpath = fileparts(xlsname);
numNII = 0; %number of NIfTI files
numMat = 0; %number of Mat files
numOK = 0;
if isempty (nuisance)
    nuisanceMat = [];
end
%designMat = [];
for i=1:size(dMat,2)
    matname = deblank( dMat(i).(SNames{1}));
    isValid = false;
    if numel(SNames) > 1
        for j = 2:numel(SNames)
            b = dMat(i).(SNames{j});
            if ~isempty(b) && isnumeric(b) && isfinite(b)
                isValid = true;
            end
        end
    end
    if ~isValid
        fprintf('Warning: no valid behavioral data for %s\n',matname);
        matname = '';
    end
    % added by GY: 
    % all nusiance variables must be present for all participants!
    % exclude participants with missing nuisance variables (i.e. NaNs)
    if ~isempty (nuisance) && ~isempty (matname)
        if sum (structfun (@isnan, nuisance(i))) 
            fprintf ('Warning: excluding %s because some nuisance regressors are missing\n', matname);
            matname = '';
        end
    end
    if ~isempty(matname)
        [matname] = findMatFileSub(matname,imgpath, xlsname);
        [~, ~, ext] = fileparts(matname);
        if strcmpi('.mat',ext) || strcmpi('.hdr',ext) || strcmpi('.nii',ext) || strcmpi('.voi',ext)
            if strcmpi('.mat',ext)
                numMat = numMat + 1;
            elseif strcmpi('.hdr',ext) || strcmpi('.nii',ext) || strcmpi('.voi',ext)
                numNII = numNII + 1;
            end
            dMat(i).(SNames{1}) = matname;
            numOK = numOK + 1;
            designMat(numOK) = dMat(i); %#ok<AGROW>
            if ~isempty (nuisance)
                nuisanceMat(numOK) = nuisance(i);%#ok<AGROW>
            end
        end
    end
end
if (numNII + numMat) == 0
    error('Unable to find any of the images listed in the file %s\n',imgpath);
end
if (numNII > 0) && (numMat >0) %mixed file
    error('Error: some images listed in %s are NIfTI format, others are Mat format. Use nii_nii2mat to convert NIfTI (.nii/.hdr) images.\n',imgpath);
end
if (numNII > 0)
    fprintf('Using NIfTI images. You will have more options if you use nii_nii2mat to convert NIfTI images to Mat format.\n');
    designUsesNiiImages = true;
end
%end nii_read_design()

function [fname] = findMatFileSub(fname, xpth, xlsname)
%looks for a .mat file that has the root 'fname', which might be in same
%folder as Excel file xlsname
fnameIn = fname;
[pth,nam,ext] = fileparts(fname);
if strcmpi('.nii',ext) || strcmpi('.hdr',ext) || strcmpi('.img',ext)%look for MAT file
    ext = '.mat';
    %fprintf('Excel file %s lists %s, but files should be in .mat format\n',xlsname,fnameIn);
else
    if exist(fname, 'file') == 2, return; end;
end
fname = fullfile(pth,[nam '.mat']);
if exist(fname, 'file'), return; end;
%next - check folder of Excel file
%[xpth,~,~] = fileparts(xlsname);
fname = fullfile(xpth,[nam ext]);
if exist(fname, 'file'), return; end;
fname = fullfile(xpth,[nam '.mat']);
if exist(fname, 'file'), return; end;
%next check for nii file:
fname = findNiiFileSub(fnameIn, xpth);
if exist(fname, 'file'), return; end;
fprintf('Unable to find image %s listed in %s: this should refer to a .mat (or .nii) file. (put images in same folder as design file)\n',fnameIn, xlsname);
fname = '';
%end findMatFileSub()

function [fname] = findNiiFileSub(fname, dir)
[pth,nam,~] = fileparts(fname);
fname = fullfile(pth,[nam '.nii']);
if exist(fname, 'file'), return; end;
fname = fullfile(pth,[nam '.hdr']);
if exist(fname, 'file'), return; end;
if exist(dir,'file') == 7
    pth = dir;
else
    [pth,~,~] = fileparts(dir);
end
fname = fullfile(pth,[nam '.nii']);
if exist(fname, 'file'), return; end;
fname = fullfile(pth,[nam '.hdr']);
if exist(fname, 'file'), return; end;
%findNiiFileSub