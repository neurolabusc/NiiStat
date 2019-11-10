function [file_list, number_list, idx] = nii_modality_list(name)
%Reports names of modalities available for analysis
%
%Examples
% [file_list, number_list] = nii_modality_list;
% [file_list, number_list, idx] = nii_modality_list('fmri');

idx = 0;
file_list = strvcat('lesion','cbf','rest','i3mT1','i3mT2','fa','dti','md','ttp','mtt','cbv','dtifc','fmri','fmrib','alf','palf','fa_dki','md_dki','dax_dki','drad_dki','mk_dki','kax_dki','krad_dki','kfa_dki');
number_list='';
for i = 1: size(file_list, 1)
	nam = deblank(file_list(i,:));
	number_list=[number_list, sprintf('%d',i) '=', nam ' ']; %#ok<AGROW>
end
file_list = char(file_list); %convert to char array
if nargin < 1, return; end;
len = length(name);
for i = 1: size(file_list,1)
    nam = char(deblank(file_list(i,:)));
    if length(nam) >= len
        %nam = nam((end-len+1):end);
        if strcmpi(nam, name), idx = i; end;
    end
end
%end nii_listtemplates()