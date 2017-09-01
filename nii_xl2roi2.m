function nii_xl2roi2 (fnm, sheetname)
%
% fnm : filename of Excel spreadsheet
% sheetname : name of worksheet to convert
%Example
% nii_xl2roi %use GUI
% nii_xl2roi('components.xlsx');
% nii_xl2roi('components.xlsx', 'dti');
%Example worksheet : Cell(1,1) is ROI name (e.g. 'JHU')
%Here the intensity "Addtl" will be stored for each region
% JHU	Magnitude	Addtl
% 25	0.176	0.176
% 43	0.215	0.039
% 31	0.242	0.027
% 186	0.267	0.025
%
%To create  a node file, first two columns provide node
%JHU	JHU2	Magnitude	Addtl
% 79	13	0.145	0.145
% 184	31	0.246	0.101
% 182	31	0.286	0.04
% 37	25	0.318	0.032
% 43	13	0.345	0.027
% 81	31	0.37	0.025
% 182	39	0.405	0.035
% 182	71	0.443	0.038
% 71	25	0.464	0.021
% 71	11	0.5	0.036
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
if size(raw,2) == 4
    xlSubNodz(xlsname, worksheetname)
    return;
end
nRegions = size(raw,1) - 1; %first row is header labels
if (nRegions < 1), error('Expected at least 1 region (2 rows: 1st is labels, 2nd is 1st region) in worksheet "%s" from "%s"\n',worksheetname, xlsname); end;
%extract key to find correct region of interest atlas
ROIshortName = raw{1,1}; %e.g. JHU
[ROIList, ROINumbers, ROIIndex] = nii_roi_list(ROIshortName, false) ;
if (ROIIndex < 1)
   error('Unable to find ROI named "%s"', raw{1,1});
end
roiname = deblank(ROIList(ROIIndex,:));
fprintf('%s uses template %s and with %d regions\n', worksheetname, roiname, nRegions);
niiname = [roiname, '.nii'];
if ~exist(niiname, 'file'), error('Unable to find the file "%s"', niiname); end;
%load template
hdr = spm_vol(niiname);
img = spm_read_vols(hdr);
hdr.pinfo = [1;0;0]; %slope = 1, intercept = 0
hdr.dt    =[16,0]; %16=float32
maxRoi = max(img(:));
imgPCA = zeros(size(img));
for r = 1 : nRegions
   roiNum = raw{r+1,1};
   if (isempty(roiNum) || (~isnumeric(roiNum)) || (roiNum < 1) || (roiNum > maxRoi))
      fprintf('Expected value 1..%d not "%s"', maxRoi, roiNum);
      continue;
   end
   val = raw{r+1,size(raw,2)};
   if (isempty(val) || (~isnumeric(val)) )
      fprintf('Expected numeric value not "%s"\n', raw{r+1,c});
      continue;
   end
   fprintf('region->val\t%d\t%g\n', roiNum, val);
   imgPCA(img == roiNum) = val;
   %fprintf('%s component %d %s (%d) = %g\n', worksheetname, c, str, row(1), raw{r, c+1});
end
hdr.fname = [worksheetname, '_', ROIshortName, '.nii'];
spm_write_vol(hdr,imgPCA);
%end xlSub()

function xlSubNodz(xlsname, worksheetname)
%fprintf(' You will get an error if %s does not have a worksheet named "%s" (case sensitive)\n', xlsname, worksheetname);
try
    [~, ~, raw] = xlsread (xlsname, worksheetname,'','basic');
catch
    fprintf('Unable to read worksheet "%s" from "%s"\n',worksheetname, xlsname);
    return;
end;
if size(raw,2) ~= 4
   error('xlSubNodz expects 4 columns [ROI1, ROI2, unused, link-weight]');
end
nRegions = size(raw,1) - 1; %first row is header labels
if (nRegions < 1), error('Expected at least 1 region (2 rows: 1st is labels, 2nd is 1st region) in worksheet "%s" from "%s"\n',worksheetname, xlsname); end;
%extract key to find correct region of interest atlas
ROIshortName = raw{1,1}; %e.g. JHU
[ROIList, ROINumbers, ROIIndex] = nii_roi_list(ROIshortName, false) ;
if (ROIIndex < 1)
   error('Unable to find ROI named "%s"', raw{1,1});
end
roiname = deblank(ROIList(ROIIndex,:));

fprintf('%s uses template %s with %d regions\n', worksheetname, roiname,  nRegions);
nodename = [roiname, '.node'];
if ~exist(nodename, 'file'), error('Unable to find the file "%s"', nodename); end;
n = load(nodename);
maxRoi = size(n,1);
edges = zeros(maxRoi,maxRoi);
for r = 1 : nRegions
   roiNum = raw{r+1,1};
   if (isempty(roiNum) || (~isnumeric(roiNum)) || (roiNum < 1) || (roiNum > maxRoi))
      fprintf('Expected value 1..%d not "%s"', maxRoi, roiNum);
      continue;
   end
   roiNum2 = raw{r+1,2};
   if (isempty(roiNum2) || (~isnumeric(roiNum2)) || (roiNum2 < 1) || (roiNum2 > maxRoi))
      fprintf('Expected value 1..%d not "%s"', maxRoi, roiNum2);
      continue;
   end
   val = raw{r+1,size(raw,2)};
   if (isempty(val) || (~isnumeric(val)) )
      fprintf('Expected numeric value not "%s"\n', raw{r+1,c});
      continue;
   end
   fprintf('regionxregion->val\t%d\t%d\t%g\n', roiNum,roiNum2, val);
   edges(roiNum,roiNum2) = val;
   edges(roiNum2,roiNum) = val;
end
nodzname = [worksheetname, '_', ROIshortName, '.nodz'];
nii_save_nodz('jhu',edges,nodzname);
%end xlSubNodz()


