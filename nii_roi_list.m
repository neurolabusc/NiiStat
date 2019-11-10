function [file_list, number_list, idx] = nii_roi_list(name, isVerbose) 
%Reports names of text files in the template folder
% name : name (or number) of region of interest
% isVerbose : (optional) if TRUE (default) than location of ROIs repoted
%used by nii_nii2mat and nii_stat
%Examples
% [files, num] = nii_roi_list
% [file_list, number_list, idx] = nii_roi_list('jhu');
idx = 0;
pth = fileparts(which(mfilename));
sub = [pth  filesep 'roi' filesep];
if ~exist('isVerbose','var') || isVerbose
    fprintf('Using regions of interest from folder %s\n',sub);
end
s = dir( [sub '*.txt']); %look in 'templates' subdirectory
file_list = {s.name}';
%next, strip .txt
file_list = sort(file_list); %alphabetical order
% the block below added by GY, March 2018: discard filenames which start with underscore
underscored_idx = [];
for i = 1: numel(file_list)
    if file_list{i}(1) == '_'
        underscored_idx = [underscored_idx i];
    end
end 
file_list(underscored_idx) = [];
number_list='';
for i = 1: numel(file_list)
	[~,nam] = fileparts(char(file_list(i)));
	file_list(i)={[sub nam]};
	number_list=[number_list, sprintf('%d',i) '=', nam ' ']; %#ok<AGROW>
end
file_list = char(file_list); %convert to char array
if nargin < 1, return; end;
if isnumeric(name) %if "nii_roi_list(3)" set idx to 3 unless out of range
    if (name < 1) || (name > size(file_list,1))
        idx = 0;
    else
        idx = name;
    end
    return;
end
len = length(name);
for i = 1: size(file_list,1)
   nam = char(deblank(file_list(i,:)));
   if length(nam) >= len
    nam = nam((end-len+1):end);
    if strcmpi(nam, name), idx = i; end;
   end
end
%end nii_roi_list()
