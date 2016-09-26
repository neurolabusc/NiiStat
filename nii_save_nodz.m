function TriBin = nii_save_nodz(roiname, matvals, nodzname, logicalMask) 
%logical mask is optional: if provided unused nodes appear smaller
% roiname : region used for connectome, e.g. 'jhu', 'AICHA'
% matvals : a NxN matrix with edge values, e.g. for JHU atlas 189x189 cells
% nodzname : name to save nodz file for viewing with surf ice
% logicalMask: optional: mask identifying tested nodes (untested nodes will appear small)
%Examples
% nii_save_nodz('jhu',rand(189,189),'test.nodz')

if min(matvals(:)) == max(matvals(:)), fprintf(' No variability, will not create %s\n', nodzname); return; end;
[~, ~, ROIIndex] = nii_roi_list(roiname, false);
if ROIIndex < 1, return; end; %unable to find ROI
if ~exist('logicalMask','var') %show all nodes with equal size
    str = nii_roi2mm (ROIIndex);
    fileID = fopen(nodzname,'w');
    fprintf(fileID, str);
    fprintf(fileID, '#ENDNODE\n');
    fclose(fileID);
    dlmwrite(nodzname,matvals,'delimiter','\t','-append');
    return;
end
%below: if mask provided, first expand mask from vector of upper triangle [size= (N*(N-1))/2] to full NxN matrix
[~, cstr] = nii_roi2mm (ROIIndex);
nLabel = numel(cstr);
TriBin = triu(ones(nLabel,nLabel), 1);
numROIorVox = sum(TriBin(:));
if (numel(logicalMask) ~= numROIorVox)
    error('%s size of logical (%d) mask is incompatible with the upper triangle of a %dx%d matrix', mfilename, numel(logicalMask), nLabel, nLabel);
end
TriBin(find(TriBin>0)) = logicalMask; %#ok<FNDSB>
%shrink diameter/size of any unused nodes
% labelsOK = (sum(TriBin') > 0); %label ok if either row OR column has an entry
labelsOK =  ((sum(TriBin) > 0) + (sum(TriBin') > 0));
for i = nLabel : -1 : 1
    if ~labelsOK(i)
        %cstr(i) = []; % <- to remove unused nodes
        parts = regexp(cstr{i},'\t','split');
        %parts{4} = '0.1';
        parts{5} = '0.1';
        cstr{i}=[sprintf('%s\t',parts{1:end-1}),parts{end}];
        %fprintf('UNTESTED REGION: %s\n',parts{end});
    elseif false %next lines only for testing!
        parts = regexp(cstr{i},'\t','split');
        fprintf('TESTED REGION: %s\n',parts{end});
    end
end
fileID = fopen(nodzname,'w');
for i = 1 : nLabel
    fprintf(fileID, sprintf('%s\n',cstr{i}));
end
fprintf(fileID, '#ENDNODE\n');
fclose(fileID);
dlmwrite(nodzname,matvals,'delimiter','\t','-append');
%end saveNodzSub()