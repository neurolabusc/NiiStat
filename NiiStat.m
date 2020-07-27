function NiiStat(xlsname, roiIndices, modalityIndices,numPermute, pThresh, minOverlap, regressBehav, maskName, GrayMatterConnectivityOnly, deSkew, doTFCE, doSVM)
%Analyze MAT files
% xlsname         : name of excel file where first column is mat name for
%                   participant and subsequent columns are behavioral values
% roiIndices      : [from nii_modality_list], will be numbers like 0=voxelwise,1=brodmann,2=jhu, 3=fox, 4=tpm, 5=aal
% modalityIndices : 1=lesion,2=cbf,3=rest,4=i3mT1,5=i3mT2,6=fa,7=dti
% numPermute    : -1=FDR,0=Bonferroni, else control for familywise error based on N iterations
%                  see nii_stat_core for more details
% pThresh       : 1-tailed statistical threshold
% minOverlap    : only examine voxels/regions non-zero in this many participants
% regressBehav  : (optional) if true use lesion volume to regress behavioral data
% maskName      : (optional) name of image to mask voxelwise data
% GrayMatterConnectivityOnly : (optional) if false, DTI and resting state
%                  will examine GM <-> WM, WM <-> CSF, GM <-> CSF connections
% deSkew        : Report and attempt to correct skewed behavioral data
% doTFCE        : Apply threshold-free cluster enhancement (voxelwise only)
%Examples
% NiiStat %use graphical interface
% NiiStat('LIME.xlsx',1,1,0,0.05,1)
%test

%Added by Roger to ensure right NiiStatGUI cfg file is opened and used

%[filepath, name,ext] = which('NiiStatGUI')

global usesGUI;

if which('NiiStatGUI')
    temp = which('niistatGUI_cfg.mat');
    [filepath, name, ext] = fileparts(temp);
    GUI = load(temp);
    if GUI.GUIdata.useClassicNiiStat == 1
        usesGUI = false;
        GUI = [];
    else
        usesGUI = true;
    end
else
    usesGUI = false;
    GUI = [];
end

%manually set usesGUI
%usesGUI = false;

% added by Roger and Grigori, Nov 2018:
% if Matlab version is older than 2016, do not use GUI
[~, temp] = version;
year = str2num (temp (length(temp)-4:length(temp)));
if year < 2016
    usesGUI = false;
    GUI = [];
end

fprintf('Version 3 March 2017 of %s %s %s\n', mfilename, computer, version);
ver; %report complete version information, e.g. "Operating System: Mac OS X  Version: 10.11.5 Build: 15F34"
if ~isempty(strfind(mexext, '32')), warning('Some features like SVM require a 64-bit computer'); end;
import java.lang.*;
hemiKey = 0;
interhemi = false; %% added by GY at RD's request
statname = '';
repopath=char(System.getProperty('user.home'));
checkForUpdate(fileparts(mfilename('fullpath')));
%checkForMostRecentMatFiles(repopath)

% updating SPM temportarily disabled by GY because SPM server is down
%if isempty(which('spm')) || ~strcmp(spm('Ver'),'SPM12'), error('SPM12 required'); end;
%if (spm_update ~= 0), warning('SPM is obsolete, run "spm_update(true)"'); end;
%%%

%%Added switch by Roger
if usesGUI
    xlsname = GUI.GUIdata.excelFile;
else
    if ~exist('xlsname','var')
        [file,pth] = uigetfile({'*.xls;*.xlsx;*.txt;*.tab','Excel/Text file';'*.txt;*.tab','Tab-delimited text (*.tab, *.txt)';'*.val','VLSM/NPM text (*.val)'},'Select the design file');
        if isequal(file,0), return; end;
        xlsname=[pth file];
    end
end


if (strcmpi('ver',xlsname)), return; end; %nii_stat('ver') cause software to report version and quit
if exist(xlsname,'file') ~= 2
    error('Unable to find Excel file named %s\n',xlsname);
end
[designMat, designUsesNiiImages, minOverlapValFile, nuisanceMat] = nii_read_design (xlsname);
if ~exist('minOverlap','var')
    if isempty(minOverlapValFile)
        minOverlap = 0;
    else
        minOverlap = minOverlapValFile;
    end
end
[~, xlsname, ~] = fileparts(xlsname);
if ~exist('regressBehav','var')
   regressBehav = false;
end
if ~exist('maskName','var')
   maskName = []; %no mask
end
if ~exist('GrayMatterConnectivityOnly','var')
    GrayMatterConnectivityOnly = true;
end
if ~exist('deSkew','var')
    deSkew = false;
end
if ~exist('customROI','var')
    customROI = false;
end
if ~exist('doTFCE','var')
    doTFCE = false;
end
if ~exist('reportROIvalues','var')
    reportROIvalues = false;
end
if ~exist('numPermute','var')
   numPermute = 0;
end
if ~exist('pThresh','var')
   pThresh = 0.05;
end
if ~exist('doSVM','var')
    doSVM = false;
end
doVoxReduce = false;
[kROIs, kROInumbers] = nii_roi_list();
[~, kModalityNumbers] = nii_modality_list();
if ~exist('modalityIndices','var') %have user manually specify settings
    cfg_filename = 'niistat_cfg.mat';
    prompt = {'Number of permutations (-1 for FDR, 0 for Bonferroni, large number for permute (3000), very small number for FreedmanLane(-3000):','Corrected P theshold:',...
        'Minimum overlap (1..numSubj):',...
        ['ROI (0=voxels ' sprintf('%s',kROInumbers) ' negative for correlations [multi OK]'],...
        ['Modality (' sprintf('%s',kModalityNumbers) ') [multiple OK]'],...
        'Special (1=explicit voxel mask, 2=control for lesion volume, 3=de-skew, 4=include WM/CSF connectivity, 5=customROI, 6=TFCE, 7=reportROImeans, 8=SVM, 9=LowRes, 10=LH only, 11=RH only; 12=interhemispheric) [multi OK]',...
        'Statistics name [optional]'
        };
    dlg_title = ['Options for analyzing ' xlsname];
    num_lines = 1;
    def = [];
    if exist(cfg_filename,'file')
        s = load(cfg_filename);
        if isfield(s,'answer')
            def = s.answer;
        end;
    end;
    if numel(def) ~= 7
      def = {'0','0.05','2','3','1','',''};
    end
    if minOverlap > 0
       def{3} = [num2str(minOverlap)];
    end
    if designUsesNiiImages
        def{4} = ['UNUSED (design file specifies ', num2str(numel(designMat)), ' voxelwise images)'];
        def{5} = def{4};

    end
    
    
   %Added switch by Roger
   
    if usesGUI
        numPermute = GUI.GUIdata.numPermutations;
        pThresh = GUI.GUIdata.correctedP;
        minOverlap = GUI.GUIdata.minOverlap;
        if ~designUsesNiiImages
            roiIndices = GUI.GUIdata.atlasChoice;
            modalityIndices = GUI.GUIdata.modalityChoice;
        end
        special = GUI.GUIdata.specialChoice;
    else 
        answer = inputdlg(prompt, dlg_title, num_lines, def);
        if isempty(answer), return; end;
        save(cfg_filename,'answer') % save user preferences
        numPermute = str2double(answer{1});
        pThresh = str2double(answer{2});minOverlap = str2double(answer{3});
        if ~designUsesNiiImages
             roiIndices = str2num(answer{4}); %#ok<ST2NM> - we need to read vectors
             modalityIndices = str2num(answer{5}); %#ok<ST2NM> - we need to read vectors
        end
        special = str2num(answer{6}); %#ok<ST2NM> - we need to read vectors
    end
   
    %reprocess ROIs
    if any(special == 13)
        
        
        
    end
        
    if any(special == 1) %select masking image
        
        %%Added switch by Roger
        if usesGUI
            maskName = GUI.GUIdata.chosenMask;
        else
            [mfile,mpth] = uigetfile('*.nii;*.hdr','Select the mask image');
            if isequal(mfile,0), return; end;
            maskName=[mpth mfile];
        end
        
    end
    if any(special == 2) %adjust behavior for lesion volume
        regressBehav = true;
    end
    if any(special == 3) %adjust behavior for skew
        deSkew = true;
    end
    if any(special == 4) %allow WM/CSF connections
        GrayMatterConnectivityOnly = false;
    end
    if any(special == 5) %allow user to specify custom ROIs
        customROI = true;
        %if (numel(roiIndices) ~= 1) || (roiIndices ~= 0)
        %    roiIndices = 0;
        %    fprintf('Custom ROIs require selecting the voxelwise modality\n');
        %end
    end
    if any(special == 6) %allow WM/CSF connections
        doTFCE = true;
    end
    if any(special == 7) %report values for each ROI
        reportROIvalues = true;
    end
    if any(special == 8) %report values for each ROI
        doSVM = true;
    end
    if any(special == 9)
        doVoxReduce = true;
    end
    if any(special == 10)
        hemiKey = 1;
    elseif any(special == 11)
        hemiKey = 2;
    end
    if any (special == 12)
        interhemi = true;
    end
    %Added switch by Roger
    if usesGUI
        statname = GUI.GUIdata.resultsFolderName;
    else
        statname = answer{7};
    end
   
end;
if designUsesNiiImages %voxelwise images do not have regions of interest, and are only a single modality
    roiIndices = 0;
    modalityIndices = 1;
end
for i = 1: length(modalityIndices) %for each modality
    modalityIndex = modalityIndices(i);
    for j = 1: length(roiIndices)
        roiIndex = roiIndices(j);
        specialStr = '';
        if exist('special','var') && ~isempty(special)
           specialStr = ['special=[', strtrim(sprintf('%d ',special)),'] '];
        end
        fprintf('Analyzing roi=%d, modality=%d, permute=%d, %sdesign=%s\n',roiIndex, modalityIndex,numPermute,specialStr, xlsname);
        %Roger added GUI as last argument
        processExcelSub(designMat, roiIndex, modalityIndex,numPermute, pThresh, minOverlap, regressBehav, maskName, GrayMatterConnectivityOnly, deSkew, customROI, doTFCE, reportROIvalues, xlsname, kROIs, doSVM, doVoxReduce, hemiKey, interhemi, statname,GUI, nuisanceMat); %%GY
    end
end
%end nii_stat_mat()

function nii = isNII (filename)
%returns true if filename is .nii or .hdr file
[~, ~, ext] = fileparts(filename);
nii = strcmpi('.voi',ext) || strcmpi('.hdr',ext) || strcmpi('.nii',ext);
%end isNII()

% function [designMat, designUsesNiiImages] = readDesign (xlsname)
% designUsesNiiImages = false;
% [~,~,x] = fileparts(xlsname);
% if strcmpi(x,'.tab') || strcmpi(x,'.txt')  || strcmpi(x,'.val')
%     dMat = nii_tab2mat(xlsname);
% else
%     dMat = nii_xls2mat(xlsname , 'Data (2)','', true);
% end
% SNames = fieldnames(dMat);
% numFields = length (SNames);
% if numFields < 2
%     error('File %s must have multiple columns (a column of file names plus a column for each behavior\n', xlsname);
% end
% numNII = 0; %number of NIfTI files
% numMat = 0; %number of Mat files
% numOK = 0;
% %designMat = [];
% for i=1:size(dMat,2)
%     matname = deblank( dMat(i).(SNames{1}));
%     isValid = false;
%     if numel(SNames) > 1
%         for j = 2:numel(SNames)
%             b = dMat(i).(SNames{j});
%             if ~isempty(b) && isnumeric(b) && isfinite(b)
%                 isValid = true;
%             end
%         end
%     end
%     if ~isValid
%         fprintf('Warning: no valid behavioral data for %s\n',matname);
%         matname = '';
%     end
%     if ~isempty(matname)
%
%         [matname] = findMatFileSub(matname,xlsname);
%
%         [~, ~, ext] = fileparts(matname);
%
%         if strcmpi('.mat',ext) || strcmpi('.hdr',ext) || strcmpi('.nii',ext)
%             if strcmpi('.mat',ext)
%                 numMat = numMat + 1;
%             elseif strcmpi('.hdr',ext) || strcmpi('.nii',ext)
%                 numNII = numNII + 1;
%             end
%             dMat(i).(SNames{1}) = matname;
%             numOK = numOK + 1;
%             designMat(numOK) = dMat(i); %#ok<AGROW>
%
%         end
%     end
% end
% if (numNII + numMat) == 0
%     error('Unable to find any of the images listed in the file %s\n',xlsname);
% end
% if (numNII > 0) && (numMat >0) %mixed file
%     error('Error: some images listed in %s are NIfTI format, others are Mat format. Use nii_nii2mat to convert NIfTI (.nii/.hdr) images.\n',xlsname);
% end
% if (numNII > 0)
%     fprintf('Using NIfTI images. You will have more options if you use nii_nii2mat to convert NIfTI images to Mat format.\n');
%     designUsesNiiImages = true;
% end
% %end readDesign()

%Roger added GUI as input argument
function processExcelSub(designMat, roiIndex, modalityIndex,numPermute, pThresh, minOverlap, regressBehav, mask_filename, GrayMatterConnectivityOnly, deSkew, customROI, doTFCE, reportROIvalues, xlsname, kROIs, doSVM, doVoxReduce, hemiKey, interhemi, statname,GUI, nuisanceMat) %%GY
%GrayMatterConnectivityOnly = true; %if true, dti only analyzes gray matter connections
%kROIs = strvcat('bro','jhu','fox','tpm','aal','catani'); %#ok<*REMFF1>
%kModalities = strvcat('lesion','cbf','rest','i3mT1','i3mT2','fa','dti','md'); %#ok<REMFF1> %lesion, 2=CBF, 3=rest
[kModalities, ~] = nii_modality_list();
if (modalityIndex > size(kModalities,1)) || (modalityIndex < 1)
    error('%s error: modalityIndex must be a value from 1..%d\n',mfilename,size(kModalities,1));
    return;
end
if roiIndex < 0
    kAnalyzeCorrelationNotMean = true;
    roiIndex = abs(roiIndex);
else
    kAnalyzeCorrelationNotMean = false;
end
if strcmpi('dtifc',deblank(kModalities(modalityIndex,:))) %read connectivity triangle
    kAnalyzeCorrelationNotMean = true;
end
if strcmpi('dti',deblank(kModalities(modalityIndex,:))) %read connectivity triangle
    kAnalyzeCorrelationNotMean = true;
end
if strcmpi('rest',deblank(kModalities(modalityIndex,:))) %read connectivity triangle
    kAnalyzeCorrelationNotMean = true;
end
if kAnalyzeCorrelationNotMean
   fprintf('analysis of connectivity between regions rather than mean intensity\n');
else
    interhemi = false; % interhemispheric analysis possible only for connectomes! --GY
end 
if roiIndex == 0 %voxelwise lesion analysis
   ROIfield = deblank(kModalities(modalityIndex,:));
else
    if doVoxReduce
        %fprintf('doVoxReduce disabled: only for voxelwise analyses\n');
        doVoxReduce = false;
    end
    if doTFCE
        fprintf('doTFCE disabled: threshold free cluster enhancement for voxelwise analyses only\n');
        doTFCE = false;
    end
    if (roiIndex > size(kROIs,1)) || (roiIndex < 0)
        fprintf('%s error: for modality %d the roiIndex must be a value from 1..%d\n',mfilename,modalityIndex,size(kROIs,1));
        return;
    end
    [~,nam] = fileparts(deblank(kROIs(roiIndex,:)));
    
    ROIfield = [deblank(kModalities(modalityIndex,:)) '_' nam];
    
end
if ~exist('statname','var') || isempty(statname)
    statname = [ROIfield '_' xlsname];%sprintf ('%s%s',deblank(kModalities(modalityIndex,:)),deblank(kROIs(roiIndex,:)));
end;
SNames = fieldnames(designMat);
matnames = [];
for i=1:size(designMat,2)
    matnames = strvcat(matnames, deblank( designMat(i).(SNames{1})) ); %#ok<REMFF1>
end
designMat = rmfield(designMat,SNames{1}); %remove first column - mat name
% read in the image data
if roiIndex == 0
    subfield = '.dat';
elseif kAnalyzeCorrelationNotMean
    subfield = '.r';
else
    subfield = '.mean';
end
subfield = [ROIfield subfield];
%explicit voxel mapping
requireVoxMask =  (roiIndex == 0) && exist('mask_filename','var') && ~isempty(mask_filename); %apply explicit masking image
if (requireVoxMask) && (doTFCE == 1)
   error('Explicit mask not (yet) compatible with TFCE');
end

%for large voxel datasets - first pass to find voxels that vary
voxMask = [];
%if false
matVer = inf;
if (requireVoxMask) || ((~customROI) && (roiIndex == 0) && (size(matnames,1) > 10) && (doTFCE ~= 1)) %voxelwise, large study
    fprintf('Generating voxel mask for large voxelwise statistics\n');
    idx = 0;
    for i = 1:size(matnames,1)
        [in_filename] = deblank(matnames(i,:));
        img = [];
        if isempty(in_filename)
            %warning already generated
        elseif isNII (in_filename)
            %error('Please use nii_nii2mat before conducting a large voxelwise statistics');
            %hdr = spm_vol (in_filename);
            %img = spm_read_vols (hdr);
            [hdr, img] = read_volsSub (in_filename);
        elseif (exist (in_filename, 'file'))
            dat = load (in_filename);
            matVer = matVerSub(dat, matVer);
            if  issubfieldSub(dat,subfield)
                hdr = dat.(ROIfield).hdr;
                img = dat.(ROIfield).dat;
            else
                fprintf('Warning: File %s does not have data for %s\n',in_filename,subfield);
            end
        end
        if ~isempty(img)
            if doVoxReduce
                [hdr, img] = resliceVolSub(hdr, img); %#ok<ASGLU>
            end
            %store behavioral and relevant imaging data for ALL relevant valid individuals
            idx = idx + 1;
            if idx == 1
                voxMask = zeros(size(img));
            end
            img(~isfinite(img)) = 0;
            img(img ~= 0) = 1;
            if numel(voxMask) ~= numel(img), s = dir(in_filename); error('Unexpected image dimensions vary %s %d', in_filename, s.bytes); end;
            voxMask  = voxMask + img;
        end;
    end %for each individual
    voxMask(voxMask < minOverlap) = 0;
    voxMask(voxMask > 0) = 1;
    if requireVoxMask
        %mask_hdr = spm_vol (mask_filename);
        %mask_img = spm_read_vols (mask_hdr);
        [mask_hdr, mask_img] = read_volsSub (mask_filename);
        mask_img(isnan(mask_img)) = 0; %exclude voxels that are not a number
        if ~isequal(mask_hdr.mat, hdr.mat) || ~isequal(mask_hdr.dim(1:3), hdr.dim(1:3))
            fprintf('WARNING: mask dimensions differ from data: attempting to reslice (blurring may occur)\n');
            inimg = mask_img; %reshape(mask_img,mask_hdr.dim(1:3)); %turn 1D vector into 3D
            imgdim = hdr.dim(1:3);
            mask_img = zeros(imgdim);
            for i = 1:imgdim(3)
                M = inv(spm_matrix([0 0 -i])*inv(hdr.mat)*mask_hdr.mat); %#ok<MINV>
                mask_img(:,:,i) = spm_slice_vol(inimg, M, imgdim(1:2), 1); % 1=linear interp; 0=nearest neighbor
            end %for each slice
        end %if dimensions differ
        fprintf('Explicit voxel mask includes %d of %d voxels\n',sum(mask_img(:) > 0), numel(mask_img));
        voxMask(mask_img == 0) = 0;
    end
    %voxMask = voxMask(:); %make 1d
    nOK = sum(voxMask(:) > 0);
    fprintf('%d of %d voxels (%g%%) show signal in at least %d participants\n',nOK, numel(voxMask),100*nOK/numel(voxMask), minOverlap );
    if nOK < 1
        error('No voxels survive in mask');
    end
end
idx = 0;
for i = 1:size(matnames,1)
    [in_filename] = deblank(matnames(i,:));
    if isempty(in_filename)
        %warning already generated
    elseif (exist (in_filename, 'file'))
        if isNII (in_filename)
            idx = idx + 1;
            data = [];
            %data.lesion.hdr = spm_vol (in_filename);
            %data.lesion.dat = spm_read_vols (data.lesion.hdr);
            [data.lesion.hdr, data.lesion.dat] = read_volsSub (in_filename);
            if doVoxReduce
                    [data.lesion.hdr, data.lesion.dat] = resliceVolSub(data.lesion.hdr, data.lesion.dat); %#ok<ASGLU>
            end
            if ~isempty(voxMask)
                data.lesion.dat = data.lesion.dat(voxMask == 1); %#ok<AGROW>
            end
            data.filename = in_filename;
            data.behav = designMat(i); % <- crucial: we inject behavioral data from Excel file!
            if ~isempty (nuisanceMat)
                data.nuisance = nuisanceMat(i);
            end
            subj_data{idx} = data; %#ok<AGROW>
        else
            %dat = load (in_filename, subfield);
            dat = load (in_filename);
            matVer = matVerSub(dat, matVer);
            [dat, cbfMean, cbfStd] = cbf_normalizeSub(dat, subfield);
            %if  issubfieldSub(dat,'lesion.dat')
            %	fprintf ('Volume\t%g\tfor\t%s\n',sum(dat.lesion.dat(:)), in_filename);
            %end
            %if  isfield(dat,subfield) % && ~isempty (data.behav)
            if (roiIndex > 0) && (~kAnalyzeCorrelationNotMean) && ~issubfieldSub(dat,subfield)
                voxField = [deblank(kModalities(modalityIndex,:)) '.dat'];
                if  issubfieldSub(dat,voxField) %we can generate ROI data from voxel data
                    fprintf('Creating %s for %s\n',subfield,in_filename);
                    %dat.(deblank(kModalities(modalityIndex,:))).hdr
                    roiName = deblank(kROIs(roiIndex,:)) ;
                    sn=[deblank(kModalities(modalityIndex,:)) '_'];
                    nii_roi2stats (roiName, dat.(deblank(kModalities(modalityIndex,:))).hdr, dat.(deblank(kModalities(modalityIndex,:))).dat, sn,in_filename);
                    dat = load (in_filename);
                end
            end
            if  issubfieldSub(dat,subfield)
                %store behavioral and relevant imaging data for ALL relevant valid individuals
                if doVoxReduce
                    [dat.(ROIfield).hdr, dat.(ROIfield).dat] = resliceVolSub(dat.(ROIfield).hdr, dat.(ROIfield).dat); %#ok<ASGLU>
                end

                idx = idx + 1;
                subj_data{idx}.filename = in_filename; %#ok<AGROW>
                subj_data{idx}.behav = designMat(i); %#ok<AGROW>
                if ~isempty (nuisanceMat)
                    subj_data{idx}.nuisance = nuisanceMat(i);%#ok<AGROW>
                end
                if isempty(voxMask)
                    %dat.(ROIfield).mean = normSub(dat.(ROIfield).mean, cbfMean, cbfStd);
                    subj_data{idx}.(ROIfield)  = dat.(ROIfield); %#ok<AGROW>
                else
                    dat.(ROIfield).dat = normSub(dat.(ROIfield).dat, cbfMean, cbfStd);
                    subj_data{idx}.(ROIfield).hdr = dat.(ROIfield).hdr; %#ok<AGROW>
                    subj_data{idx}.(ROIfield).dat = dat.(ROIfield).dat(voxMask == 1); %#ok<AGROW>
                end

                if regressBehav && isfield (dat.lesion, 'dat')
                    dat.lesion.dat(isnan(dat.lesion.dat(:)))=0; %zero NaNs: out of brain
                    subj_data{idx}.lesion.vol = sum(dat.lesion.dat(:)); %#ok<AGROW>
                end
                if (idx == 1) && (roiIndex < 1) %first image of voxelwise analyses
                    vox = numel(subj_data{idx}.(ROIfield).dat(:));
                    vox = vox * size(matnames,1); %worst case scenario: all individuals have image data
                    gb = (vox * 8)/ (1024^3); %doubles use 8-bytes
                    fprintf('The imaging data will require %.3f gb of memory\n',gb);
                end
            else
                fprintf('Warning: File %s does not have data for %s\n',in_filename,subfield);
            end
        end
    else
        fprintf('Unable to find file %s\n', in_filename);
    end
end
matVerCheckSub(matVer);
clear('dat'); %these files tend to be large, so lets explicitly free memory
n_subj = idx;
if n_subj < 3
    fprintf('Insufficient data for statistics: only found files %d with both "behav" and "%s" fields\n',n_subj,ROIfield);
    return;
end
% get the list of numeric fields of behavioural data
fields = fieldnames (subj_data{1}.behav); %fields = fieldnames (data.behav);
idx = 1;
for i = 1:length (fields)
    for s = 1:n_subj
        if isnumeric (subj_data{s}.behav.(fields{i}))
            beh_names{idx} = fields{i};  %#ok<AGROW>
            idx = idx + 1;
            break
        end
    end
end
n_beh = idx - 1;
if ~exist('beh_names','var')
    fprintf('No valid behavioral variables found\n');
    return
end
%beh_names = [];beh_names{1} = 'ASRS_total';n_beh = 1;fprintf('WARNING: Beta release (single behavior)#@\n');%#@
% make sure all the subjects have all numeric fields
beh = zeros(n_subj,n_beh);
beh(:) = nan;
nuisance = [];
for i = 1 : n_subj
    for j = 1:n_beh %length(beh_names)
        if isfield (subj_data{i}.behav, beh_names{j})
            if ~isnumeric(subj_data{i}.behav.(beh_names{j}) )
                fprintf ('Warning! Subject %s reports non-numeric data for field %s\n',subj_data{i}.filename, beh_names{j} );
            elseif strcmpi(subj_data{i}.behav.(beh_names{j}),'NaN') || (isnan(subj_data{i}.behav.(beh_names{j}) ))
                fprintf ('Warning! Subject %s reports NaN for field %s\n',subj_data{i}.filename, beh_names{j} );
            else
                beh(i, j) = subj_data{i}.behav.(beh_names{j});
                %fprintf('%d %d %f\n',i,j, beh(i, j));
                %class(beh(i, j))
            end
        else
            disp (['Warning! Subject ' subj_data{i}.filename ' does not have a field ' beh_names{j}]);
        end
    end
    
    if ~isempty (nuisanceMat)
        nuisance (i, :) = structfun (@(x) x, subj_data{i}.nuisance); %#ok<AGROW>
    end
    
end
   

if regressBehav
    vol = zeros(n_subj,1);
    vol(:) = nan;
    for i = 1:n_subj
        %subj_data{idx}.lesion.vol
        if isfield (subj_data{i}.lesion, 'vol')
            vol(i) = subj_data{i}.lesion.vol;
            fprintf ('Participant\t%s\tVolume\t%g\n',subj_data{i}.filename,vol(i));
        else
            fprintf ('Problem regressing for lesion volume! Subject %s does not have the field ".lesion.dat"\n', subj_data{i}.filename);
        end;
    end;
    if sum(~isnan(vol(:))) > 1
        nuisance = [nuisance vol];
         
%         for i = 1:n_beh
%             %beh_names1 = deblank(beh_names(i));
%             beh1 = beh(:,i);
%             good_idx = intersect (find(~isnan(beh1)), find(~isnan(vol)));
%             dat = beh1(good_idx)'; %behavior is the data
%             reg = vol(good_idx)'; %lesion volume is our regressor
%             preSD = std(dat);
%             if ~isnan(std(dat)) && (preSD ~= 0) && (std(reg) ~= 0) %both data and regressor have some variability
%                 G = ones (2, size(dat,2)); %constant
%                 G (2, :) = reg; % linear trend
%                 G_pseudoinv = G' / (G * G'); %aka: G_pseudoinv = G' * inv (G * G');
%                 Beta = dat * G_pseudoinv;
%                 dat = dat - Beta*G; %mean is zero
%                 fprintf('Regressing %s with lesion volume reduces standard deviation from %f to %f\n',char(deblank(beh_names(i))),preSD, std(dat) );
%                 beh(good_idx,i) = dat;
%             end
%         end
    end
end %if regressBehav - regress behavioral data using lesion volume

roiName = '';
if roiIndex == 0 %voxelwise lesion analysis
    les_names = [];
    hdr = subj_data{1}.(ROIfield).hdr;
    if ~isempty (voxMask) %10/16 added by CR
        logicalMask = logical (ones (numel(voxMask), 1));
    end
    for i = 1:n_subj
        if (i > 1) && (numel(subj_data{i}.(ROIfield).dat(:)) ~= numel(subj_data{1}.(ROIfield).dat(:)))
            %error('Number of voxels varies between images. Please reslice all images to the same dimensions');
            Interp = ~isBinSub(subj_data{i}.(ROIfield).dat); %interpolate continuous images, do not interpolate binary images
            fprintf('warning: reslicing %s to match dimensions of other images. Interpolation = %d\n',subj_data{i}.filename, Interp);
            [~, outimg] = nii_reslice_target(subj_data{i}.(ROIfield).hdr, subj_data{i}.(ROIfield).dat(:), subj_data{1}.(ROIfield).hdr, Interp) ;
            subj_data{i}.(ROIfield).dat = outimg; %#ok<AGROW>
            %fprintf('%d, %d\n',subj_data{i}.filename);

        end
        %fprintf('%d/%d= %d\n',i,n_subj, numel(subj_data{i}.(ROIfield).dat(:)));

        les(i, :) = subj_data{i}.(ROIfield).dat(:); %#ok<AGROW>
        if ~exist('logicalMask','var') || isempty (logicalMask)
            logicalMask = logical (ones (size (les, 2), 1)); %%% added by GY
        end %10/16 conditional by CR

    end
    nanIndex = isnan(les(:));
    if sum(nanIndex(:)) > 0
        les(nanIndex) = 0;
        fprintf('Warning: Not a number values in images replaced with zeros\n');
    end
else %if voxelwise else region of interest analysis
    if  exist('mask_filename','var') && ~isempty(mask_filename)
       error('Explicit mask only for voxelwise data');
    end
    %NEXT roi masking
    roiMaskI = []; %inclusive ROI mask - e.g. if [1 2 12] then only these 3 regions analyzed...
    les_names = cellstr(subj_data{1}.(ROIfield).label); %les_names = cellstr(data.(ROIfield).label);
   
    if customROI
        
        %Roger Added the switches below based on the GUI.GUIdata through
        
        if length(GUI) == 0
       
            answer = inputdlg(['Only include the following regions (1..', sprintf('%d',numel(les_names)),')'] ,'ROI inclusion criteria',1,{'1 2 6 8'});
            if isempty(answer), return; end;
            roiMaskI = str2num(answer{1});
       
        end
        
        if length(GUI) > 0
            
            if GUI.GUIdata.useClassicNiiStat == 0
              if strfind(ROIfield,char(GUI.GUIdata.atlas1name)) > 0
                     roiMaskI = GUI.GUIdata.atlas1picks; 
                  end
                   if strfind(ROIfield,char(GUI.GUIdata.atlas2name)) > 0
                     roiMaskI = GUI.GUIdata.atlas2picks; 
                  end
                   if strfind(ROIfield,char(GUI.GUIdata.atlas3name)) > 0
                     roiMaskI = GUI.GUIdata.atlas3picks; 
                  end
                   if strfind(ROIfield,char(GUI.GUIdata.atlas4name)) > 0
                     roiMaskI = GUI.GUIdata.atlas4picks; 
                   end
                   if strfind(ROIfield,char(GUI.GUIdata.atlas5name)) > 0
                     roiMaskI = GUI.GUIdata.atlas5picks; 
                  end
                   if strfind(ROIfield,char(GUI.GUIdata.atlas6name)) > 0
                     roiMaskI = GUI.GUIdata.atlas6picks; 
                   end
                   if strfind(ROIfield,char(GUI.GUIdata.atlas7name)) > 0
                     roiMaskI = GUI.GUIdata.atlas7picks; 
                  end
                   if strfind(ROIfield,char(GUI.GUIdata.atlas8name)) > 0
                     roiMaskI = GUI.GUIdata.atlas8picks; 
                  end    

%                 if strfind(ROIfield,'fox') > 0 %lesion_fox
%                     roiMaskI = GUI.GUIdata.fox_picks
%                 end
%                 if strfind(ROIfield,'aal') > 0 %lesion_aal
%                     roiMaskI = GUI.GUIdata.aal_picks
%                 end
%                 if strfind(ROIfield,'aalcat') > 0%lesion_aalcat
%                     roiMaskI = GUI.GUIdata.aalcat_picks
%                 end
%                 if strfind(ROIfield,'jhu') > 0 %lesion_jhu
%                     roiMaskI = GUI.GUIdata.jhu_picks
%                 end
%                 if strfind(ROIfield,'AICHA') > 0 
%                     roiMaskI = GUI.GUIdata.aicha_picks
%                 end
%                 if strfind(ROIfield,'bro') > 0 %lesion_bro
%                     roiMaskI = GUI.GUIdata.brodmann_picks
%                 end
%                 if strfind(ROIfield,'catani') > 0%lesion_catani
%                     roiMaskI = GUI.GUIdata.catani_picks
%                 end
%                 if strfind(ROIfield,'cat') > 0%lesion_cat
%                     roiMaskI = GUI.GUIdata.cat_picks
%                 end
                 %roiMaskI = GUI.GUIdata.atlas1picks;
            end
        end
                
    else
        global global_roiMask
        roiMaskI = global_roiMask; %inclusion mask
    end

    %find the appropriate ROI
    %[mpth,~,~] = fileparts( deblank (which(mfilename)));
    %roiName = fullfile(mpth,[deblank(kROIs(roiIndex,:)) '1mm.nii']);
    roiName = [deblank(kROIs(roiIndex,:)) '.nii'];
    if ~exist(roiName,'file')
        fprintf('No images created (unable to find image named %s\n',roiName);
        return;
    end
    hdr = spm_vol(roiName);
    %provide labels for each region
    %les_names = cellstr(subj_data{1}.(ROIfield).label); %les_names = cellstr(data.(ROIfield).label);
    %next: create labels for each region, add image values


    if hemiKey > 0 && ~interhemi
        if isempty (roiMaskI)
            roiMaskI = extract_hemi_idxSub (les_names, hemiKey);
            customROI = 1;
        else
            roiMaskI = intersect (roiMaskI, extract_hemi_idxSub (les_names, hemiKey));
        end
    end


    if kAnalyzeCorrelationNotMean %strcmpi('dti',deblank(kModalities(modalityIndex,:))) %read connectivity triangle

        if GrayMatterConnectivityOnly
            GM_mask = get_GM_Sub(les_names);
            GM_idx = find (GM_mask);
            if isempty (roiMaskI)
                roiMaskI = GM_idx;
            else
                roiMaskI = intersect (roiMaskI, GM_idx);
            end
        end


        labels = les_names;
        for i = 1:n_subj
            %http://stackoverflow.com/questions/13345280/changing-representations-upper-triangular-matrix-and-compact-vector-octave-m
            %extract upper triangle as vector
            A = subj_data{i}.(ROIfield).r;
            %%% commented out by GY
            %             A = shrink_matxCustomSub( A,  roiMaskI);
            %             if GrayMatterConnectivityOnly == true
            %                 [les_names,A] = shrink_matxSub(labels,A);
            %                 %fprintf('Only analyzing gray matter regions (%d of %d)\n',size(les_names,1),size(labels,1) );
            %             end
            B = triu(ones(size(A)),1);
            les(i, :) = A(B==1); %#ok<AGROW>
            if interhemi %% added by GY at RD's request
                C = zeros (size (A));
                L_idx = extract_hemi_idxSub (les_names, 1);
                R_idx = extract_hemi_idxSub (les_names, 2);
                roiMask_L = intersect (roiMaskI, L_idx);
                roiMask_R = intersect (roiMaskI, R_idx);
                C (roiMask_L, roiMask_R) = 1;
                C (roiMask_R, roiMask_L) = 1;
                if hemiKey == 1
                    C (roiMask_L, roiMask_L) = 1;
                elseif hemiKey == 2
                    C (roiMask_R, roiMask_R) = 1;
                end
                D = triu (ones (size (C)), 1);
                logicalMask = logical (C (D == 1));
            else
                if isempty (roiMaskI)
                    logicalMask = logical (ones (size (les, 2), 1));
                else
                    C = ones (size (A));
                    to_exclude = setdiff (1:length(labels), roiMaskI);
                    C (to_exclude, :) = 0; C (:, to_exclude) = 0;
                    D = triu (ones (size (C)), 1);
                    logicalMask = logical (C (D == 1));
                end
            end
            % A=[0 1 2 4; 0 0 3 5; 0 0 0 6; 0 0 0 0];  B = triu(ones(size(A)),1); v =A(B==1); v = 1,2,3,4,5,6
        end
        if GrayMatterConnectivityOnly
            fprintf('Connectivity only analyzing gray matter regions (%d of %d)\n',size(les_names,1),size(labels,1) );
        end
    else %not DTI n*n connectivity matrix
        for i = 1:n_subj
            les(i, :) = subj_data{i}.(ROIfield).mean;

            %%%% commented out by GY
            %             if ~isempty(roiMaskI)
            %                 A = subj_data{i}.(ROIfield).mean;
            %                 les(i, :) = A(roiMaskI); %#ok<AGROW> %%% NOTE TO GY
            %             else
            %                 les(i, :) = subj_data{i}.(ROIfield).mean;     %#ok<AGROW>
            %             end
        end
        %%% added by GY
        if isempty (roiMaskI)
            logicalMask = logical (ones (size (les, 2), 1));
        else
            logicalMask = logical (zeros (size (les, 2), 1));
            logicalMask (roiMaskI) = 1;
        end
    end

end %if voxelwise else roi

if customROI && roiIndex == 0
    %if roiIndex ~= 0, fprintf('Custom ROIs require selecting the voxelwise modality\n'); end;
    roiNames = spm_select(inf,'image','Select regions of interest');
    lesVox = les;
    les = zeros(n_subj, size(roiNames,1) );
    for i = 1:n_subj
        [les(i, :), les_names] = nii_nii2roi(roiNames,hdr,lesVox(i, :));
    end
    hdr = []; %no image for these regions of interest
    les_names = cellstr(les_names);
end %if custom ROI


if (numPermute < -2) && (numPermute >= -500)
    fprintf('Error: Current software can not understand %d permutations (reserved for future usage).\n', numPermute);
    return;
end
if ((size(beh,2) <= 1) || sum(isnan(beh(:)))) > 0 && (numPermute < -500)
    fprintf('Error: Freedman-Lane requires at least two columns of behavioral data and no empty cells.\n');
    return;
end
if deSkew
    for i =1:n_beh
        if isBinomialSub(beh(:,i))
            fprintf('Behavior %s is binomial\n',beh_names{i});
        else %if binomial else continuous
            sk = zskewSub(beh(:,i));
            %transform skewed data http://fmwww.bc.edu/repec/bocode/t/transint.html
            if abs(sk) < 1.96
                fprintf('Behavior %s has a Z-skew of %f\n',beh_names{i}, sk);
            else %if not skewed else transfrom
                mn = min(beh(:,i));
                beh(:,i) = beh(:,i) - mn; %24Sept2014 - previously would crash with negative values, e.g. sqrt(-3) fails isreal
                if sk > 1.96
                    beh(:,i) = sqrt(beh(:,i));
                else %negative skew
                    beh(:,i) = beh(:,i).^2;
                end
                skT = zskewSub(beh(:,i));
                fprintf('Behavior %s had a Z-skew of %f, after transform this became %f\n',beh_names{i}, sk, skT);
            end %if not significantly skewed else transform
        end %if binomial else continuous
    end %for each behavior
end %if de-Skew
if doTFCE
    hdrTFCE = hdr.dim;
else
    hdrTFCE = [];
end


if (reportROIvalues) && (numel(les_names) < 1)
    fprintf('Unable to create a ROI report [voxelwise analyses]\n');
elseif (reportROIvalues) && (kAnalyzeCorrelationNotMean)
    fprintf('Unable to create a ROI report [correlation matrix analyses]\n');
elseif reportROIvalues
    %note this next conditional removes regions with little variability.
    %  n.b. this same step is built into nii_stat_core, but we will do it here
    %  so reportROIvalues will match what will be computed
    if (minOverlap > 1) && (numel(les_names) > 0)
        nOK = 0;
        for j = 1:numel(les_names)
            if (sum ((les(:, j)) ~= 0) > minOverlap)
                nOK = nOK + 1;
                les(:, nOK) = les(:, j);
                les_names{nOK} = les_names{j};
            end
        end %for j
        if nOK < 1
           error('No regions non-zero in at least %d individuals', minOverlap);
        end
        if nOK < numel(les_names)
            fprintf('%d of %d regions non-zero in at least %d individuals\n', nOK, numel(les_names), minOverlap);
            les_names = les_names(1:nOK);
            les = les(:,1:nOK);
        end
    end %if minOverlap > 1 and not voxelwise
    %first row: column labels
    fprintf('filename\t');
    for j = 1:numel(les_names)
         fprintf('%s\t', les_names{j});
    end
    for j = 1:n_beh %length(beh_names)
        fprintf('%s\t', beh_names{j});
    end
    fprintf('\n');
    for i = 1:n_subj
        fprintf('%s\t',subj_data{i}.filename);
        for j = 1:numel(les_names)
             fprintf('%g\t',les(i, j));
        end
        for j = 1:n_beh %length(beh_names)
           if isnan(beh(i, j))
            fprintf('\t');
           else
            fprintf('%g\t',beh(i, j));
           end
        end
        fprintf('\n');
    end
    return; %no analysis - just report values
end


%%%% GY: moved min_overlap selection from nii_stat_core
%next: identify which voxels/regions should be analyzed
bad_idx = union (find (isnan (sum (abs(les), 1))), find(var(les,0,1)<eps)); %eliminate voxels/regions with no variability
if minOverlap > 0  %isBinomialLes
    bad_idx = union (bad_idx, find (sum ((les ~= 0), 1) < minOverlap)); %eliminate voxels/regions with no variability
end

% %%% PLEASE DELETE
% disp ('ZZZZZZHOPA');
% length (logical_mask)
% %%%



%%% the following line added by GY
logicalMask (bad_idx) = 0;
good_idx = setdiff (1:size(les, 2), bad_idx);
if length(good_idx) < 1 %no surviving regions/voxels
    if isBinomialLes
        error('%s error. no voxels damaged in at least %d participants.',mfilename,minOverlap);
    else
        error('%s error: no regions to analyze (voxels are either not-a-number or have no variability).',mfilename);
    end
end

% moved here from nii_stat_core by GY
chDirSub(statname);
diary ([deblank(statname) '.txt']);

if minOverlap > 0 %isBinomialLes
    fprintf('Only analyzing voxels non-zero in at least %d individuals.\n',minOverlap);
end

% the following fprintf's slightly modified by GY
if size(beh,2) == 1
    fprintf('**** Analyzing %s with %d participants for behavioral variable %s across %d (of %d) regions/voxels.\n',deblank(statname),size(les,1),deblank(beh_names{1}), sum(double(logicalMask)), size(les,2));
else
    fprintf('**** Analyzing %s with %d participants for %d behavioral variables across %d (of %d) regions/voxels.\n',deblank(statname),size(les,1),size(beh,2), sum(double(logicalMask)), size(les,2));
end

if ~isempty (voxMask)
    if numel (voxMask) ~= numel (logicalMask)
        error ('Something is very wrong: voxMask and logicalMask don''t match in size %d ~= %d', numel (voxMask),numel (logicalMask) );
    end
    zero_idx = find (voxMask == 0);
    logicalMask (zero_idx) = 0;
    voxMask = [];
end
if sum(isnan(beh(:))) > 0
    for i = 1 : n_beh
        fprintf('Behavior %d/%d: estimating behaviors one at a time (removing empty cells will lead to faster analyses)\n',i,n_beh);
        beh_names1 = deblank(beh_names(i));
        beh1 = beh(:,i);
        good_idx = find(~isnan(beh1));
        beh1 = beh1(good_idx);
        les1 = zeros(length(good_idx),size(les,2));
        for j = 1:length(good_idx)
            les1(j, :) = les(good_idx(j), :) ;
            %les1(j,1) = beh1(j); %to test analyses
        end
        %save(sprintf('nii_stat%d',i)); %<-troublshoot, e.g. load('nii_stat4');  nii_stat_core(les1, beh1, beh_names1,hdr, pThresh, numPermute, logicalMask,statname, les_names,hdrTFCE);
        %localMask = var(les1(:,logicalMask)) ~= 0;
        %localMask = var(les1(:,logicalMask)) ~= 0;
        localMask = var(les1) ~= 0;
        if minOverlap > 1
            %localMaskMinOverlap = sum(les1(:,logicalMask) ~= 0) > minOverlap;
            localMaskMinOverlap = sum(les1 ~= 0) > minOverlap;
            localMask(~localMaskMinOverlap) = false;
        end
        logicalMask1 = logicalMask;
        if numel(logicalMask1) == numel(localMask)
            logicalMask1(~localMask) = 0;
        elseif sum(logicalMask1) == numel(localMask)  %localMask generated on compressed dataset - expand it
            les1 = les1 (:, localMask); %squeeze data to only examine critical voxels
            localMaskExpanded = zeros(size(logicalMask1));
            localMaskExpanded(find(logicalMask1)) = localMask; %#ok<FNDSB>
            logicalMask1(~localMaskExpanded) = 0;
        else
            error ('Something is very wrong: logicalMask1 and localMask don''t match in size %d (or %d) ~= %d', numel(logicalMask1), sum(logicalMask1), numel(localMask) );
        end

        %if any(localMask == false) %regions that have variability overall do not have variability for this factor
        %    idx = find(logicalMask);
        %    logicalMask1(idx(find(localMask == false))) = false; %#ok<FNDSB>
        %end

        if doSVM
            nii_stat_svm(les1, beh1, beh_names1,statname, les_names, subj_data, roiName, logicalMask1, hdr, pThresh, numPermute, nuisance);
        else
            %nii_stat_core(les1, beh1, beh_names1,hdr, pThresh, numPermute, minOverlap,statname, les_names,hdrTFCE, voxMask);
            nii_stat_core(les1, beh1, beh_names1,hdr, pThresh, numPermute, logicalMask1,statname, les_names,hdrTFCE, nuisance);
        end

%         diary off
%         cd .. %leave the folder created by chDirSub

        %fprintf('WARNING: Beta release (quitting early, after first behavioral variable)#@\n');return;%#@
    end
else
    %for aalcat we may want to remove one hemisphere
    %les_names(:,1:2:end)=[]; % Remove odd COLUMNS: left in AALCAT: analyze right
    %les(1:2:end)=[]; % Remove odd COLUMNS: left in AALCAT: analyze right
    %les_names(2:2:end)=[]; % Remove even COLUMNS: right in AALCAT: analyze left
    %les(:,2:2:end)=[]; % Remove even COLUMNS: right in AALCAT: analyze left
    if doSVM
        nii_stat_svm(les, beh, beh_names, statname, les_names, subj_data, roiName, logicalMask, hdr, pThresh, numPermute, nuisance);
    else
        nii_stat_core(les, beh, beh_names,hdr, pThresh, numPermute, logicalMask,statname, les_names, hdrTFCE, nuisance);
    end
end

% moved here from nii_stat_core by GY
diary off
cd .. %leave the folder created by chDirSub

%end processMatSub()

% the following function is not used -- GY
function  mat = shrink_matxCustomSub( mat, roiMaskI)
%removes columns/rows as specified by the global "roiMask"
%this next bit allow us to remove ROIs
% global roiMask
% roiMask = [1 3 4]; %only include these regions of interest
if isempty(roiMaskI), return; end;
%indx = roiMask; %[1 2 5] preserves 1st, 2nd 5th
smallmat = mat(roiMaskI,:); %remove rows
mat = smallmat(:,roiMaskI); %remove columns
%if shrinkLabels == 1, labels = labels(roiMaskI,:); end;
%end shrink_matxCustomSub()


% the following function is not used -- GY
function [smalllabels, smallmat] = shrink_matxSub(labels, mat)
%removes columns/rows where label does not end with text '|1'
%  useful as the labels end with |1, |2, |3 for gray matter, white matter and CSF
% l = strvcat('ab|1', 'abbs|2','c|1')
% m = [1 2 3; 4 5 6; 7 8 9];
% [sl,sm] = shrink_matxSub(l,m);
% s now = [1 3; 7 9] - removes columsn and row without |1
index = strfind(cellstr(labels),'|1');
index = ~cellfun('isempty',index);
if (sum(index(:)) == 0)
    smallmat = mat;
    smalllabels = labels;
    fprintf(' Analysis will include all regions (this template does not specify white and gray regions)\n') %% GY
    return
end;
smallmat = mat(index,:);
smallmat = smallmat(:,index);
smalllabels = labels(index,:);
%end shrink_matxSub()


% modified version: return a binary mask of regions
% where 1 is a GM region, and 0 is WM or CSF -- GY
function GM_mask = get_GM_Sub(labels)
% returns
index = regexp (cellstr (labels), '\|1$'); % --GY
%index = strfind(cellstr(labels),'|1'); %% commented out by GY
index = ~cellfun('isempty',index);
if (sum(index(:)) == 0)
    smalllabels = labels;
    GM_mask = ones (length(labels), 1);
    fprintf(' Analysis will include all regions (this template does not specify white and gray regions)\n') %% GY
    return
end;
smalllabels = labels(index,:);
GM_mask = index;
%end get_GM_Sub()

% function added by GY: returns indices of left-hemisphere (or
% right-hemishere) regions
function hemi_idx = extract_hemi_idxSub (labels, hemiKey)
if hemiKey == 1
    hemi_regexp = {'_L$', '-L\|', '_L\|', '_Left\|', 'left\|'};
elseif hemiKey == 2
    hemi_regexp = {'_R$', '-R\|', '_R\|', '_Right\|', 'right\|'};
end
hemi_idx = [];
for i = 1:length (hemi_regexp)
    hemi_idx = [hemi_idx; find(cellfun(@length, regexp (labels, hemi_regexp{i})))];
end
hemi_idx = unique (hemi_idx);
%end extract_hemi_idxSub


% function [fname] = findMatFileSub(fname, xlsname)
% %looks for a .mat file that has the root 'fname', which might be in same
% %folder as Excel file xlsname
% fnameIn = fname;
% [pth,nam,ext] = fileparts(fname);
% if strcmpi('.nii',ext) || strcmpi('.hdr',ext) || strcmpi('.img',ext)%look for MAT file
%     ext = '.mat';
%     %fprintf('Excel file %s lists %s, but files should be in .mat format\n',xlsname,fnameIn);
% else
%     if exist(fname, 'file') == 2, return; end;
% end
% fname = fullfile(pth,[nam '.mat']);
% if exist(fname, 'file'), return; end;
% %next - check folder of Excel file
% [xpth,~,~] = fileparts(xlsname);
% fname = fullfile(xpth,[nam ext]);
% if exist(fname, 'file'), return; end;
% fname = fullfile(xpth,[nam '.mat']);
% if exist(fname, 'file'), return; end;
% %next check for nii file:
% fname = findNiiFileSub(fnameIn, xlsname);
% if exist(fname, 'file'), return; end;
% fprintf('Unable to find image %s listed in %s: this should refer to a .mat (or .nii) file. (put images in same folder as design file)\n',fnameIn, xlsname);
% fname = '';
% %end findMatFileSub()

% function [fname] = findNiiFileSub(fname, dir)
% [pth,nam,~] = fileparts(fname);
% fname = fullfile(pth,[nam '.nii']);
% if exist(fname, 'file'), return; end;
% fname = fullfile(pth,[nam '.hdr']);
% if exist(fname, 'file'), return; end;
% if exist(dir,'file') == 7
%     pth = dir;
% else
%     [pth,~,~] = fileparts(dir);
% end
% fname = fullfile(pth,[nam '.nii']);
% if exist(fname, 'file'), return; end;
% fname = fullfile(pth,[nam '.hdr']);
% if exist(fname, 'file'), return; end;
% %findNiiFileSub

function b = isBinomialSub(i)
%returns true if vector is binomial (less than three distinct values)
nMin = sum(i(:)==min(i(:)));
nMax = sum(i(:)==max(i(:)));
if (nMin+nMax) ~= length(i(:))
    b = false;
else
    b = true;
end
%end isBinomialSub()

function s = zskewSub(i)
%http://office.microsoft.com/en-us/excel-help/skew-HP005209261.aspx
%zSkew dividing the Skew by the Standard Error of the Skew
%standard error of Skewness http://en.wikipedia.org/wiki/Talk%3ASkewness,
%http://www.unesco.org/webworld/idams/advguide/Chapt3_1_3.htm
%http://jalt.org/test/PDF/Brown1.pdf -> Tabachnick and Fidell, 1996
n = numel(i);
mn = mean(i);
s = std(i);
if (n < 3) || (s == 0)
    s = 0;
    return
end
s = sum(((i-mn)/s).^3);
s = n/((n-1)*(n-2)) * s;
s = s/(sqrt(6/n)); %convert skew to Z-skew
%end zskewSub()

function img = normSub(img, mn, stdev)
if ~isfinite(mn) || ~isfinite(stdev), return; end;
%img = ((img-mn)/stdev)+1;
img = ((img-mn)/stdev)+100;
%normSub

function [dat, cbfMean, cbfStd] = cbf_normalizeSub(dat,subfield)
%normalize intensity values
% with SPM and ASLtbx voxels outside brain are zero or nan
cbfMean = nan;
cbfStd = nan;
if isempty(dat), return; end;
if isempty(strfind(subfield,'cbf')), return; end;
if ~issubfieldSub(dat,'cbf'), return; end;
if issubfieldSub(dat.cbf,'c1R') && issubfieldSub(dat.cbf,'c2R')
    cbfMean = (dat.cbf.c1R+dat.cbf.c2R)/2;
    cbfStd = 0.66;
    fprintf('Normalizing intensity using precomputed right hemisphere CBF estimates %g\n', cbfMean);
    return
end
fprintf('WARNING: DO NOT IGNORE! Please use modern CBF maps')
img = dat.cbf.dat(:);
img(~isfinite(img)) = 0; %remove nan, inf, -inf voxels
imgL=img;
imgL(1:floor(size(img,1)/2),:,:) = 0;
imgL = imgL(imgL ~= 0);
imgR=img;
imgR(ceil(size(img,1)/2):end,:,:) = 0;
imgR = imgR(imgR ~= 0);
if mean(imgR) > mean(imgL)
    cbfMean = mean(imgR);
    %cbfMean = median(imgR);
    %cbfStd = std(imgR);
    cbfStd = 30;
    str = 'LEFTdark';
else
    cbfMean = mean(imgL);
    %cbfMean = median(imgL);
    %cbfStd = std(imgL);
    cbfStd = 30;
    str = 'RIGHTdark';
end
fprintf('%s Lmean/Lstd/Lmdn/Rmean/Rstd/Rmdn:\t%g\t%g\t%g\t%g\t%g\t%g\n',str,mean(imgL), std(imgL), median(imgL), mean(imgR), std(imgR), median(imgR));
%end cbf_normalizeSub()

function [r] = issubfieldSub(s, f)
%https://fieldtrip.googlecode.com/svn/trunk/utilities/issubfield.m
try
  getsubfieldSub(s, f);    % if this works, then the subfield must be present
  r = true;
catch %#ok<CTCH>
  r = false;                % apparently the subfield is not present
end
%end issubfieldSub()

function [s] = getsubfieldSub(s, f)
% GETSUBFIELD returns a field from a structure just like the standard
% Matlab GETFIELD function, except that you can also specify nested fields
% using a '.' in the fieldname. The nesting can be arbitrary deep.
%
% Use as
%   f = getsubfield(s, 'fieldname')
% or as
%   f = getsubfield(s, 'fieldname.subfieldname')
%
% See also GETFIELD, ISSUBFIELD, SETSUBFIELD
% Copyright (C) 2005, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: getsubfield.m 7123 2012-12-06 21:21:38Z roboos $
if iscell(f)
  f = f{1};
end
if ~ischar(f)
  error('incorrect input argument for fieldname');
end
t = {};
while (1)
  [t{end+1}, f] = strtok(f, '.'); %#ok<AGROW,STTOK>
  if isempty(f)
    break
  end
end
s = getfield(s, t{:});
%end getsubfieldSub()

function [outhdr,outimg] = resliceVolSub(hdr, img)
%reslices image to new resolution and bounding box
%Example
% hdr = spm_vol('wT1_P098.nii')
% img = spm_read_vols(hdr);
% [hdr,img] = resliceVol(hdr, img)
if size(img,4) > 1, error('resliceVol fails with 4D data'); end;
outhdr = hdr;
outhdr.mat = [-2 0 0 80; 0 2 0 -114; 0 0 2 -52; 0 0 0 1]; %<- 2mm
outhdr.dim = [79 95 69]; %<-2mm,
outimg = zeros(outhdr.dim);
for i = 1: outhdr.dim(3)
    M = inv(spm_matrix([0 0 -i])*inv(outhdr.mat)*hdr.mat); %#ok<MINV>
    outimg(:,:,i) = spm_slice_vol(img, M, outhdr.dim(1:2), 1); % (0=nearest, 1=linear interp)
end
%end resliceVolSub()

function isBin = isBinSub(x)
nMin = sum(x(:)==min(x(:)));
nMax = sum(x(:)==max(x(:)));
%isBin = ((nMax + nMin) == numel(x));
isBin = ((nMax + nMin + sum(isnan(x(:)))) == numel(x));
%end isBinSub

function checkForUpdate(repoPath)
prevPath = pwd;
cd(repoPath);
if exist('.git','dir') %only check for updates if program was installed with "git clone"
    [s, r] = system('git fetch origin','-echo');
    if strfind(r,'fatal')
        warning('Unable to check for updates. Network issue?');
        cd(prevPath); %CR 8/2016
        return;
    end
    [~, r] = system('git status','-echo');
    if strfind(r,'behind')
        if askToUpdate
            system('git reset --hard HEAD');
            [~, r] = system('git pull','-echo');
            showRestartMsg
        end
    end
else %do nothing for now
    warning(sprintf('To enable updates run "!git clone git@github.com:neurolabusc/%s.git"\n or "!git clone git clone https://github.com/neurolabusc/%s.git"',mfilename,mfilename));
end
cd(prevPath);
%end checkForUpdate()

function showRestartMsg
uiwait(msgbox('The program must be restarted for changes to take effect. Click "OK" to quit the program. You will need to restart it just as you normally would','Restart Program'))
exit;
%end showRestartMsg()

function a = askToUpdate
% Construct a questdlg
choice = questdlg(sprintf('An update for %s is available. Would you like to update?',mfilename), ...
	'Auto update', ...
	'Yes','No','Yes');
% Handle response
switch choice
    case 'Yes'
        a = true;
    case 'No'
        a = false;
end
%end askToUpdate()

% function h = showDownloading
% h = msgbox('Downloading matfiles...','Downloading');
% %end showRestartMsg()
%
% function showUnzipMsg(pth)
% fprintf('Go to file: %s\n\nUnzip and enter password',pth);
% uiwait(msgbox(sprintf('Go to file: %s\n\nUnzip and enter password',pth),'Unzip files'));
% %end showRestartMsg()

% function a = askToUpdateMatFiles
% % Construct a questdlg
% choice = questdlg(sprintf('Would you like to download the most recent .mat files? You must also need a password to use the files'), ...
% 	'Auto update', ...
% 	'Yes','No','Yes');
% % Handle response
% switch choice
%     case 'Yes'
%         a = true;
%     case 'No'
%         a = false;
% end
% %end askToUpdateMatFiles()

% function checkForMostRecentMatFiles(repoPath)
% repo = 'NiiMatFiles';
% prevPath= pwd;
% if exist(repoPath,'dir')
%     cd(repoPath);
% end
% if askToUpdateMatFiles
%     dlh = showDownloading;
%     pthToFile = websave(repo,'http://people.cas.sc.edu/rorden/matfiles/current.zip');
%     delete(dlh);
%     showUnzipMsg(pthToFile);
% end
% cd(prevPath);
% %checkForMostRecentMatFiles()

function matVer = matVerSub(mat, matVer)
%return minimum mat.T1.lime across a series of matlab files
if ~isfield(mat,'T1') || ~isfield(mat.T1,'lime') || (matVer <= mat.T1.lime), return; end;
matVer = mat.T1.lime;
%matVerSub()

function matVerCheckSub(matVer)
if ~isfinite(matVer), return; end;
url = 'http://people.cas.sc.edu/rorden/matfiles/index.html';
try
	str = urlread(url);
catch
	warning('Unable to verify lime files: could not connect to connect to %s',url);
	return;
end
key = '<a href="M.';
pos = strfind(str,key);
if isempty(pos), return; end;
str = str((pos(1)+numel(key)):end);
key = '.dmg"';
pos = strfind(str,key);
if isempty(pos), return; end;
str = str(1:(pos(1)-1));
vers = str2num(str);
if isempty(vers), return; end;
if vers <= matVer, fprintf('Your LIME mat files are up to date\n'); return; end;
urlupdate = fullfile(fileparts(url), ['M.', str, '.dmg']);
warning('You have mat files %4.4f, the current version is %4.4f', matVer, vers);
fprintf('Go to <a href="%s">%s</a> for mat files %4.4f\n', urlupdate, urlupdate, vers);
%matVerCheckSub()

% moved from nii_stat_core by GY
function chDirSub(statname)
datetime=datestr(now);
datetime=strrep(datetime,':',''); %Replace colon with underscore
datetime=strrep(datetime,'-','');%Replace minus sign with underscore
datetime=strrep(datetime,' ','_');%Replace space with underscore
newdir = [statname '_' datetime ];
mkdir(newdir);
cd(newdir);
%chDirSub()


function [hdr, img] = read_volsSub (fnm)
[fnm, isGz] = unGzSub (fnm); %convert FSL .nii.gz to .nii
hdr = spm_vol(fnm); %load header data
img = spm_read_vols(hdr); %load image data
if (isGz), delete(fnm); end; %remove .nii if we have .nii.gz
%end read_volsSub()

function [fnm, isGz] = unGzSub (fnm)
[pth,nam,ext] = spm_fileparts(fnm);
isGz = false;
if strcmpi(ext,'.gz') %.nii.gz
    fnm = char(gunzip(fnm));
    isGz = true;
elseif strcmpi(ext,'.voi') %.voi ->
    onam = char(gunzip(fnm));
    fnm = fullfile(pth, [nam '.nii']);
    movefile(onam,fnm);
    isGz = true;
end;
%end unGzSub()
