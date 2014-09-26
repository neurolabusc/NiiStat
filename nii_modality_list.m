function [file_list, number_list] = nii_modality_list() 
%Reports names of modalities available for analysis
%
file_list = strvcat('lesion','cbf','rest','i3mT1','i3mT2','fa','dti','md','ttp','cbv','dtifc');
number_list='';
for i = 1: size(file_list, 1)
	nam = deblank(file_list(i,:));
	number_list=[number_list, sprintf('%d',i) '=', nam ' ']; %#ok<AGROW>
end
file_list = char(file_list); %convert to char array
%end nii_listtemplates()