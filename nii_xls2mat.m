function [stat] = nii_xls2mat(xlsname, worksheetname, tagname, firstColumnText)
%Parses behavioral data from a xls or xlsx format worksheet
%  xlsname : name of Excel file to open (.xls or .xlsx)
%  worksheetname : name of relevant sheet (page) in Excel file
%  tagname : (optional) this script searches for the row with this text in the first column
%            if tagname not specified all rows are read
%  firstColumnText : (optional) if true then 1st column always returned as
%              string. E.G. '748' will have class=char, not class=double
% Grigori Yourganov and Chris Rorden, 1/2014 distributed under GNU General Public Licence (version 2).
%Example
% s = nii_xls2mat('LIME_12_16_2013.xlsx','Data (2)','1001');
% s = nii_xls2mat('LIMEpf2.xlsx','Data (2)') %read all rows!

stat = [];
if exist(xlsname,'file') ~= 2
    fprintf('Unable to find file named "%s"\n',xlsname);
    return
end
if ~exist('firstColumnText', 'var')
    firstColumnText = false;
end
%if verLessThan('matlab', '7.14')
%	fprintf('You are running an old version of Matlab that has a fragile xlsread. If the script crashes, use Excel to save as an Excel 5.0/''95 format file\n');
%end
% read the Excel file into cell arrays
fprintf(' You will get an error if your Excel file does not have a worksheet named "%s" (case sensitive)\n', worksheetname);
try
    [~, txt, raw] = xlsread (xlsname, worksheetname);
catch
    fprintf('Unable to read worksheet "%s" from "%s"\n',worksheetname, xlsname);
    return;
end;
%[~, txt, raw] = xlsread (xlsname, 'NiiStat','','basic');

% prepare the list of subject ID's:
% they are presumed to be in the first column of the Excel file;
% if they are numeric, convert them to strings.
id_list_raw = raw (:, 1);
isnum = cellfun (@isnumeric, id_list_raw);
num_idx = find (isnum == 1); % indices of numeric IDs
str_idx = find (isnum == 0); % indices of string IDs
id_list_str (str_idx) = id_list_raw (str_idx);
id_list_str (num_idx) = cellfun (@num2str, id_list_raw(num_idx), 'UniformOutput',false);
% find the entry for the required ID
if exist('tagname','var') && ~isempty(tagname)
    subj_idx = find (strcmp (id_list_str, tagname));
    if isempty(subj_idx)
        fprintf('Unable to find 1st column named "%s" in worksheet "%s" of file "%s"\n',tagname, worksheetname, xlsname);
        return
    end
else
    %return all valid rows
    subj_idx = [];
    for i = 2: length(id_list_str)
        if ~strcmpi(id_list_str(i),'NaN')
            subj_idx = [subj_idx, i];
        end
    end
end
% get list of tag names from the first row of the Excel file
header_list = txt (1, :);
for i = 1:length (header_list)
    h = header_list{i}; %header
    if ~isempty(h)
        h = CleanStrSubDot(h);
        for j = 1: length(subj_idx)
            v = raw{subj_idx(j), i}; %value
            if i == 1 && firstColumnText %first first column to be char
                if isnumeric (v)
                    stat(j).(h) = num2str(v);
                else
                    stat(j).(h) = v; %not numeric, e.g. 'anomic'
                end
            else %attempt to convert string to number
                if isa(v,'double') %Windows computers return doubles, OSX returns strings?
                    n = v;
                else
                    n = str2double(v);
                end
                if isnan(n)
                    stat(j).(h) = v; %not numeric, e.g. 'anomic'
                else
                    stat(j).(h) = n; %numeric value, e.g. '0.82'
                end
            end
        end
    end
end
%end nii_xls2mat()

%  subfunctions follow
function outStr = CleanStrSubDot(inStr) %for field names
outStr = regexprep(inStr,' ','_');
outStr = regexprep(outStr,'\\','_');
outStr = regexprep(outStr,'/','_');
outStr = regexprep(outStr,'@','_');
outStr = regexprep(outStr,':','');
outStr = regexprep(outStr,'-','m');
outStr = regexprep(outStr,'\.','p'); %Illegal to have field names with a dot
%end CleanStrSubDot()
