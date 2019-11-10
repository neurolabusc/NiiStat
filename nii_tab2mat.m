function [stat, CritN] = nii_tab2mat(tabname, tagname)
%Given a tab-delimited file, returns all columns of row 'tagname', or all rows and columns if 'tagname' not provided 
% tabname: tab delimited text file with top row of labels and 1st column of tags
% tagname: desired row, e.g. 'A1', leave empty to return entire matrix
%	ID	WPM	DIAG
%	A1	1.2	Anomic			
%	B3	2.7	Brocas
%	C5	4.7	Brocas
% Note the software ignores initial lines that being with the # character:
%  #comments
%  PatientName	Behavior1
%  pat1.mat	32
%  pat2.mat	12
%  pat3.mat	10
%Example
% behavioralData2MatSub('~/lime/behav.tab','B3'); %returns B3|2.7|Brocas

CritN = [];
if ~exist('tabname','var')  %file not specified
   [fnm,pth] = uigetfile({'*.txt;*.tab;';'*.*'},'Select tab-delimited text file');
   tabname = [pth, fnm];
end;
if ~exist('tagname','var')  %file not specified
    tagname = [];
end
if isempty(tabname), return; end
if (exist(tabname,'file') ~= 2)
    fprintf('Unable to find tab-delimited text file %s\n',tabname);
    return;
end
stat = [];
fid = fopen(tabname);
while(1) %skip lines that begin with # (comments)
    hdrline = fgetl(fid); % Get second row (first row of data)
    key = '#CritPct';
    if strncmpi(hdrline,key,numel(key))
        parts = regexp(hdrline,'\t','split');
        if numel(parts) > 1
            CritPct = str2num(parts{2}); %#ok<ST2NM>
        end
    end

    if(length(hdrline)==1) && (hdrline==-1) % Reached end of file, terminate
        fprintf('%s: file does not appear to be in val format %s\n', mfilename,valname);
        break
    end
    if hdrline(1)~='#'
        break
    end
end
%hdrline = fgetl(fid); % Get second row (first row of data)
tabLocs=findstr(char(9),hdrline); % find the tabs
tabN = length(tabLocs);
if tabN < 1
    fprintf('%s file %s is not tab-delimited table\n', mfilename,tabname);
    fclose(fid);
    return;
end
header=textscan(hdrline,'%s',(tabN+1),'delimiter','\t');
outrow = 0;
while(1)
	datline = fgetl(fid); % Get second row (first row of data)
	if(length(datline)==1) && (datline==-1) % Reached end of file, terminate
        if ~isempty(tagname)
            fprintf('%s unable to find identifier %s in file %s\n', mfilename,tagname,tabname);
        elseif outrow == 0
            fprintf('%s unable to find rows in file %s\n', mfilename,tabname);
        end
        break
    end
	dat=textscan(datline,'%s',(tabN+1),'delimiter','\t');
    if isempty(tagname) || strcmpi(tagname, deblank(dat{1}{1})) %strcmp works regardless of string length, better than (tagname == deblank(dat{1}{1}))
        outrow = outrow + 1; %add another row
        for t = 1: size(dat{1},1) %excel does not put tabs for empty cells (tabN+1)
            h = CleanStrSubDot(header{1}{t}); %header
            v = CleanStrSub(dat{1}{t}); %value
            if ~isempty(v) 
                n = str2double(v);
                if (t == 1) || (isnan(n))
                    stat(outrow).(h) = v; %not numeric, e.g. 'anomic'
                else
                    stat(outrow).(h) = n; %numeric value, e.g. '0.82'
                end
            end
        end
        if strcmp(tagname, deblank(dat{1}{1})) %found our tag
            break  
        end
    end
end
fclose(fid);
if exist('CritPct','var')
    CritN = round(CritPct/100 * outrow);
end

%end nii_tab2mat()

%%%%  subfunctions follow

function outStr = CleanStrSubDot(inStr) %for field names
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