function  nii_val2mat(valname, modalityIndex, isNormalizedToStrokeControl)
%convert NPM/VLSM val file to TAB/MAT files for analysis with nii_stat
% valname    : name of text file that lists image names and respective
%Examples
% nii_stat_val2mat('behav.val');
if ~exist('valname','var')  
   [nam,pth] = uigetfile({'*.val;';'*.*'},'Select VAL file'); 
   if isequal(nam,0), return; end;
   valname =fullfile (pth, nam);
end
if exist(valname,'file') == 0
    fprintf('Unable to find a file named %s\n',valname);
    return;
end
[imgnames, beh_names, beh] = nii_val2matSub(valname);
%assume images are in same folder as val file
origpath = pwd;
[pth,nam,~] = spm_fileparts(deblank(valname));
tabname = fullfile(pth, [nam '.tab']);
if ~isempty(pth), cd(pth); end;
if ~exist('isNormalizedToStrokeControl','var')
    nii_nii2mat(imgnames); %convert image data
else
    nii_nii2mat(imgnames, modalityIndex, isNormalizedToStrokeControl); %convert image data
end
reportDesignSub (imgnames, beh_names, beh, tabname);
fprintf('Created design file named %s that can be analyzed with nii_stat\n',tabname);
cd(origpath); %back where we started
%end nii_stat_val2mat()

%--------SUBFUNCTIONS FOLLOW-------------

function reportDesignSub (imgnames, beh_names, beh, tabname)
%report design of data
if exist('tabname','var')
    fid = fopen(tabname, 'w'); %write to disk
else
    fid = 1; %write to screen
    
end
tab = sprintf('\t');
behLabel = 'ImageName';
for j = 1:size(beh_names,2)
        behLabel = [behLabel, tab, beh_names{j}];
end
fprintf(fid, '%s\n',behLabel);
for i =1:size(imgnames,1) %report design
    imgname = deblank(imgnames(i,:));
    [pth,nam,ext] = spm_fileparts(imgname);
    if strcmpi(ext,'.gz') %strip both extensions from .nii.gz
        imgname = fullfile(pth, nam );
        [pth,nam,ext] = spm_fileparts(imgname);
    end
    if strcmpi(ext,'.hdr') || strcmpi(ext,'.nii') || strcmpi(ext,'.img') || strcmpi(ext,'.voi') || strcmpi(ext,'.gz')
        imgname =fullfile(pth, nam); %strip extension: we will read .mat files, not NIfTI
    end;
    tabname = fullfile(pth, [nam '.tab']);
    behLabel = [];
    for j = 1:size(beh,2)
        	behLabel = [behLabel, tab, num2str(beh(i,j))];
    end
    fprintf(fid, '%s%s\n',imgname,behLabel);
end
if fid ~= 1
    fclose(fid);
end
%end reportDesignSub()

function [imgnames, behavLabels, behav] = nii_val2matSub(valname);
%Returns list of row headers (imgnames), row headers (behavLabels) and numeric matrix  
% valname: tab delimited text in NPM/VLSM format
behavLabels = [];
imgnames = [];
if (exist(valname) ~= 2)
    fprintf('Unable to find tab-delimited val file %s\n',valname);
    return;
end
fid = fopen(valname);
while(1)
    hdrline = fgetl(fid); % Get second row (first row of data)
    if(length(hdrline)==1) && (hdrline==-1) % Reached end of file, terminate
        fprintf('%s: file does not appear to be in val format %s\n', mfilename,valname);
        break
    end
    if hdrline(1)~='#'
        break
    end
end
tabLocs=findstr(char(9),hdrline); % find the tabs
tabN = length(tabLocs);
if tabN < 1
    fprintf('%s file %s is not tab-delimited table\n', mfilename,valname);
    fclose(fid);
    return;
end
header=textscan(hdrline,'%s',tabN+1,'delimiter','\t'); %+1 fencepost error, first item header{1}(1)
for i = 1:tabN
    behavLabels{i} =  CleanStrSubDot (char(header{1}(i+1)));
end
numSubj = 0;
while(1)
	datline = fgetl(fid); % Get second row (first row of data)
	if(length(datline)==1) && (datline==-1) % Reached end of file, terminate
        fprintf('%s found %d images each with %d behavioral variables in file %s\n', mfilename,numSubj,tabN, valname);
        fclose(fid);
        break
    end
    numSubj = numSubj + 1; 
	dat=textscan(datline,'%s',tabN+1,'delimiter','\t'); %+1 fencepost error
    imgnames = strvcat(imgnames, char(dat{1}(1)));
    for i = 1:tabN
        behav(numSubj,i) = nan; %for all empty cells, including those with no trailing tabs
    end
    nCol = min(tabN,size(dat{1},1)-1 ); %warning, Excel omits trailing tabs
    for i = 1:nCol
        v = CleanStrSub(dat{1}(i+1)); %+1 since first coumn is image name
        n = str2double(v);
        if isnan(n)
            fprintf('%s error, unable to convert "%s" to a number (row "%s", column "%s" in file %s)\n', mfilename,char(v),imgnames(numSubj,:), behavLabels{i}, valname);
        else
            behav(numSubj,i) = n; %numeric value, e.g. '0.82'
        end
    end	
end %for each line
%end nii_val2matSub()

function outStr = CleanStrSubDot(inStr)
%'Wab IQ: 2' -> 'WabIQ_2'
outStr =regexprep(inStr,' ','_');
outStr =regexprep(outStr,':','');
outStr =regexprep(outStr,'-','m');
outStr =regexprep(outStr,'\.','p'); %Illegal to have field names with a dot
%end CleanStrSubDot()

function outStr = CleanStrSub(inStr)
%'Wab IQ: 2' -> 'WabIQ_2'
outStr =regexprep(inStr,' ','_');
outStr =regexprep(outStr,':','');
%end CleanStrSub()