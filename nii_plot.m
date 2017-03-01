function nii_plot(fnm, roiIndex)
%Compute number of people required to replicate observed effects
%(optional) inputs:
% fnm : filename of power.mat file created by NiiStat
% roiIndex : index of region of interest
%Note
% uses 'power.mat' files created with 2017 and later versions of NiiStat
%

if (nargin < 2)  
    prompt = {'ROI number'};
    dlg_title = 'Plot analysis preferences';
    num_lines = 1;
    def = {'7'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer), return; end;
    roiIndex = str2double(answer{1});
end;
if ~exist('fnm','var'), fnm = 'power.mat'; end;
if ~exist(fnm,'file') 
   [A,Apth] = uigetfile({'*.mat;';'*.*'},'Select first image');
   fnm = [Apth, A];
end;
m = load(fnm);
nObs = size(m.les,1);
nVox = size(m.les,2);
nBeh = size(m.beh,2);
fprintf('%d participants, %d voxels/regions tested, %d behaviors tested\n', nObs, nVox, nBeh);
if size(m.beh,1) ~= nObs, error('behavioral data must have one entry per participant.\n'); end;
%if size(m.beh,2) ~= 1, error('%s only currently supports one behavior\n', mfilename); end;
roiGood = find(m.good_idx == roiIndex,1);
roiName = m.roi_names(roiIndex,:);
roiName = strrep(char(roiName), '_', ' ');
if isempty(roiGood) 
    error('"%s" (index %d) did not survive analysis (too few people with injury to this region?).', char(roiName), roiIndex);
end
les = m.les(:,roiGood);
%text report: you can copy into excel
str = sprintf('%s\t',roiName);
for i = 1 : nBeh
	str = [str  sprintf('%s\t',strrep(char(m.beh_names(i)), '_', ' '))];
end;
fprintf('%s\n', str);
for j = 1 : nObs
    str = sprintf('%g\t',les(j));
    for i = 1 : nBeh
        str = [str  sprintf('%g\t',m.beh(j,i) )];
    end;
    fprintf('%s\n', str);
end


for i = 1 : nBeh
    beh = m.beh(:,i);
    r = corrcoef(les, beh);
    behName = strrep(char(m.beh_names(i)), '_', ' ');
    fprintf('region "%s" (index %d) to behavior "%s" correlation r = %g\n', roiName, roiIndex, behName,  r(2));
    scatter( les, beh);
    title(sprintf('r = %g',r(2)));
    xlabel(roiName);
    ylabel(behName);
end
%end nii_plot()
