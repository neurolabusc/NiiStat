function nii_nii2mat (niinames, modalityIndex, disknames, roi)
%Convert NIfTI images to mat files for analysis with nii_stat_xls
% niinames = (optional) filenames to convert
% modality : is this lesion, cbf, etc
% disknames : (optional) custom name for mat file
%Example
% nii_nii2mat(strvcat('wa.nii','wb.nii') );
% nii_nii2mat('vx.nii',1);
% nii_nii2mat('lesionP11.nii',1,'P11.mat');
% nii_nii2mat('lesionP11.nii','lesion','P11.mat');

if ~exist('niinames','var') %no files specified
 niinames = spm_select(inf,'^.*\.(gz|voi|img|nii)$','Select images to convert to mat');
end;
if length(niinames) < 1, return; end;
[kModalities, kModalityNumbers] = nii_modality_list();
if ~exist('modalityIndex','var') %get preferences
    prompt = {['Modality (' sprintf('%s',kModalityNumbers) ')']};
    dlg_title = 'Values for adjusting the image(s)';
    num_lines = 1;
    def = {'1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer), return; end;
    modalityIndex = str2double(answer{1});
end
if ischar(modalityIndex)
    [~, ~, idx] = nii_modality_list(modalityIndex);
    if idx < 1, error('Invalid modality name %s', modalityIndex); end;
    modalityIndex = idx;
end
%end modalityIndexSub()

%set modality
if (modalityIndex > size(kModalities,1)) || (modalityIndex < 1)
    fprintf('%s error: modalityIndex must be a value from 1..%d\n',mfilename,size(kModalities,1));
    return;
end
if modalityIndex == 1
    binarize = true;
    fprintf('Images will be binarized (0 or 1), if this is not desired please edit %s\n', mfilename);
else
    binarize = false;
end
%kROIs = strvcat('bro','jhu','aal','catani'); %#ok<*REMFF1>
% kROIs = strvcat('bro','jhu','aal','catani'); %#ok<*REMFF1>
% if isNormalizedToStrokeControl
%     kROIs = strvcat(kROIs, 'aalcat');
% end
kROIs = nii_roi_list();
if exist('roi', 'var')
    roiStart = roi;
    roiEnd = roi;   
    if roi > size(kROIs,1), error('roi must be a number 1..%d',size(kROIs,1)); end; 
else
    roiStart = 1;
    roiEnd = size(kROIs,1);
end;
Voxfield = deblank(kModalities(modalityIndex,:));
%convert images
for i=1:size(niinames,1)
    niiname = deblank(niinames(i,:));
    niiname = unGzSub (niiname);
    [p,n,~] =spm_fileparts(niiname);
    if exist('disknames','var') && size(disknames,1) >= i
    	diskName = deblank(disknames(i,:));
    else
		diskName = fullfile(p,[n '.mat']);
    end
    hdr = spm_vol(niiname);
    img = spm_read_vols(hdr);
    mn = min(img(:));
    mx = max(img(:));
    if mn == mx 
        error('No variability in image %s',niiname);
    end
    if binarize
        temp = img;
        mn = mn + ((mx-mn)/2); %threshold 50% between min and max
        img(temp < mn) = 0;
        img(temp >= mn) = 1;
    end
    stat = [];
    if i == 1, fprintf('Creating modality "%s"\n',Voxfield); end;
    if ndims(img) == 3
        stat.(Voxfield).hdr = hdr;
        stat.(Voxfield).dat = img;
    end
    for roiIndex = roiStart: roiEnd
        ROIname = deblank(kROIs(roiIndex,:));
        if i == 1
            [~, nam] = fileparts(ROIname);
            fprintf('Adding region of interest "%s"\n',[Voxfield '_' nam  ]); 
        end;
        %s = nii_roi2stats('jhu','LS_LM1001.nii','','','lesion_') %3d volume
        statfield = [Voxfield '_'];
        
        new = nii_roi2stats(ROIname,hdr,img,statfield);
        stat = nii_mergestruct(stat,new);
    end;
    if exist(diskName,'file')
        old = load(diskName);
        stat = nii_mergestruct(stat,old); %#ok<NASGU>
    end
    save(diskName,'-struct', 'stat');
end;

function fnm = unGzSub (fnm)
[pth,nam,ext] = spm_fileparts(fnm);
if strcmpi(ext,'.gz') %.nii.gz
    fnm = char(gunzip(fnm));    
elseif strcmpi(ext,'.voi') %.voi -> 
    onam = char(gunzip(fnm));
    fnm = fullfile(pth, [nam '.nii']);
    movefile(onam,fnm);
end;  
%end unGzSub()