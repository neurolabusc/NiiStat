function nii_xl2roi (fnm, sheetname)
%
% fnm : filename of Excel spreadsheet
% sheetname : name of worksheet to convert
%Example
% nii_xl2roi %use GUI
% nii_xl2roi('components.xlsx');
% nii_xl2roi('components.xlsx', 'dti');
%Example worksheet (fill in a few JHU regions)
%   21|RG_L|gyrus rectus left|1 8   11
%   22|RG_R|gyrus rectus right|1    6   5
%   31|AG_L|angular gyrus left|1    4   2

if ~exist('fnm','var') || ~exist(fnm, 'file')  
    [A,Apth] = uigetfile({'*.xls;*.xlsx;';'*.*'},'Select Excel file');
    fnm = [Apth, A];
end;
[status, sheets] = xlsfinfo(fnm);
if isempty(sheets), error('Unable to read %s', fnm); return; end;
if exist('sheetname', 'var') && ~isempty(sheetname)
    idx = find(not(cellfun('isempty', strfind(sheets,sheetname))));
    if isempty(idx), error('Unable to find worksheet named %s', sheetname); end;
    xlSub(fnm, sheetname);
    return;
end;
for i = 1 : numel(sheets) %if sheetname not specified, generate once per worksheet
    xlSub(fnm, sheets{i});
end
%end nii_xl2roi()

function xlSub(xlsname, worksheetname)
%fprintf(' You will get an error if %s does not have a worksheet named "%s" (case sensitive)\n', xlsname, worksheetname);
try
    [~, ~, raw] = xlsread (xlsname, worksheetname,'','basic');
catch
    fprintf('Unable to read worksheet "%s" from "%s"\n',worksheetname, xlsname);
    return;
end;
nComponents = size(raw,2) - 1; %-1 as first column is region name
nRegions = size(raw,1);
if (nComponents < 1) or (nRegions < 2), error('Expected at least 2 regions (rows) and one component (columns)'); end;
%extract key to find correct region of interest atlas
str = raw{2,1}; %2nd row just in case 1st row is header name, e.g. 'i3mT1_AICHA_171|S_Sup_Temporal-2-L|S_Sup_Temporal-2-L'
idx = findstr(str,'|');
if numel(idx) < 2, error('Expected multiple "|" symbols ("171|S_Sup_Temporal-2-L|S_Sup_Temporal-2-L"), not "%s"',str); end;
str = str(idx(1)+1:idx(2)-1); % "171|S_Sup_Temporal-2-L|S_Sup_Temporal-2-L" -> "S_Sup_Temporal-2-L"
%determine correct ROI
[files, num] = nii_roi_list;
for i = 1 : size(files,1)
    roiname = [deblank(files(i,:)), '.txt'];
    if ~exist(roiname, 'file'), continue; end;
    roitxt = textread(roiname,'%s');
    pos = strfind(roitxt,str);
    row = find(~cellfun('isempty', pos));
    if ~isempty(row), break; end; %stop when roi contains string
end
if isempty(row), error('None of the ROIs contain the string "%s"', str); end;
fprintf('%s uses template %s and has %d components and %d regions\n', worksheetname, roiname, nComponents, nRegions);
niiname = [deblank(files(i,:)), '.nii'];
if ~exist(niiname, 'file'), error('Unable to find the file "%s"', niiname); end;
%load template
hdr = spm_vol(niiname);
img = spm_read_vols(hdr);
hdr.pinfo = [1;0;0]; %slope = 1, intercept = 0
hdr.dt    =[16,0]; %16=float32
for c = 1 : nComponents %make one image for each region
    imgPCA = zeros(size(img));
    roiRepeat = zeros(1,max(img(:))); %roiRepeat deals with the same roi being repeated on multiple rows
    for r = 1 : nRegions
       str = raw{r,1}; 
       idx = findstr(str,'|');
       if numel(idx) < 2, continue; end; %'Expected multiple "|" symbols ("171|S_Sup_Temporal-2-L|S_Sup_Temporal-2-L"), not "%s"',str); end;
       str = str(idx(1)+1:idx(2)-1); % "171|S_Sup_Temporal-2-L|S_Sup_Temporal-2-L" -> "S_Sup_Temporal-2-L"
       pos = strfind(roitxt,str);
       row = find(~cellfun('isempty', pos));
       if isempty(row), warning('Unable to find region %s', str); continue; end; %roi not found
       roiNum = row(1);
       %this next code is for sparse rois, where line number ~= index
       roiStr = roitxt{row(1)}; %e.g. '83|G_Parietal_Sup-4-L|G_Parietal_Sup-4-L'
       idx = strfind(roiStr,'|');
       if ~isempty(idx)
           roiNum2 = str2num(roiStr(1:idx(1)-1)); %e.g. '83|G...' -> '83'
           if roiNum ~= roiNum2
               warning('Assuming %s is region %d, not %d', roiStr, roiNum2, roiNum); 
               roiNum = roiNum2;
           end
       end
       %now fill in value
       imgPCA(img == roiNum) = imgPCA(img == roiNum) + raw{r, c+1};
       roiRepeat(roiNum) = roiRepeat(roiNum) + 1;
       %fprintf('%s component %d %s (%d) = %g\n', worksheetname, c, str, row(1), raw{r, c+1});
    end
    %average regions that are reported multiple times
    roiRepeat(roiRepeat == 0) = 1; %avoid divide by zero errors
    for r = 1 : numel(roiRepeat)
        imgPCA(img == r) = imgPCA(img == r)/roiRepeat(r);
    end
    %
    hdr.fname = [worksheetname, '_pca', num2str(c), '.nii'];
    spm_write_vol(hdr,imgPCA);
end
%end xlSub()

