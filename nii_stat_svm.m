function nii_stat_svm(les,beh, beh_names, statname, les_names, subj_data, roifname, logicalMask, hdr, pThresh, numPermute, nuisance)
analyzing_connectome = false;
voxelwise_analysis = isempty (les_names);

% multiple comparisons correction: default is FDR, Bonferroni is also an
% option, other ways of MC ocrrection are not implemented for SVR/SVM
if ~exist('numPermute','var')
    mcorr_method = -1; %FDR
else
    if numPermute == 0
        mcorr_method = 0; % Bonferroni
    else
        mcorr_method = -1; %FDR
    end
end
if ~exist('pThresh','var')
    pThresh = 0.05;
end
if ~exist('nuisance','var')
    nuisance = [];
end



if numel(les_names) ~= size(les,2) %for connectome analyses
    les_matrix = [];
    n = 0;
    % the for loops changed by GY
    for i = 2:length(les_names)
        for j = 1:(i-1) 
            area_name = [les_names{i} '*' les_names{j}];
            n = n+1;
            les_matrix{n} = area_name; %#ok<AGROW>
        end
    end
    analyzing_connectome = true;
    nLabel = numel(les_names);
    les_names = les_matrix;    
end

% added by GY
if ~voxelwise_analysis % if ROI analysis (as opposed to voxelwise)
    les = les (:, logicalMask);
    les_names = les_names (logicalMask);
    if numel(les_names) ~= size(les,2)
        fprintf('%s error: number of feature names does not match number of features',mfilename);
        return;
    end    
else
        % GY, Sept 2019: do nothing (the dTLVC section commented out)
%     bb = unique (les (:));
%     if length (bb) == 2 && bb(1) == 0 && bb(2) == 1 % voxelwise lesion maps
%         % Normalize the lesion maps: divide each (binary) voxel by sqrt(lesion size)
%         % So, each subject's lesion vector has unit norm;
%         % this could help reduce the effect of lesion size.
%         % See Zhang et al., "Multivariate Lesion-Symptom Mapping Using
%         % Support Vector Regression", HBM 2014, section 3.5
%         les = les ./ repmat (sqrt (sum (les, 2)), [1 size(les, 2)]);
%         les (find (isnan (les))) = 0;
%     end
end

orig_les = les;

for i = 1:length (les_names)
    fprintf ('%s\t', les_names{i});
end
fprintf ('\n');

    
if ~exist('statname','var')
    statname = 'anonymous';
end
if ~exist('subj_data','var')
    subj_data = [];
end;
chDirSub([statname '_svm']);
diary ([deblank(statname) 'svm.txt']);
for j = 1:size(beh_names,2) %for each behavioral variable
    beh_name1 = beh_names{j};
    beh1 = beh(:,j);
    
    if ~voxelwise_analysis
        [fnm, nOK] = tabFileSub(les,beh1, beh_name1,  les_names, subj_data);
    else
        nOK = 1;
    end
    
    % GY, Aug 2019: add nuisance regressors
    % NOT SURE HOW TO HANDLE THE DAMN THING.
    % For now: regress nuisances out of behavioural data for SVR,
    % and regress them out of the lesion data for SVM.
    % Not happy about this solution, will keep thinking.
    % see the article:
    % DeMarco, A. T., & Turkeltaub, P. E. (2018). Human brain mapping, 39(11), 4169-4182.
    if ~isempty (nuisance)
        if ~nii_isBinary(beh1)
            pred = [ones(size(beh, 1), 1) nuisance];
            beta = pred \ beh1;
            residuals = beh1 - pred*beta;
            beh1 = residuals;
            if ~voxelwise_analysis
                tabFileSub(les,beh1, ['residuals_' beh_name1], les_names, subj_data);
            end
        else
            les = orig_les;
            pred = [ones(size(beh, 1), 1) nuisance];
            beta = pred \ les;
            residuals = les - pred*beta;
            les = residuals;
        end
    end
    %%%%
    
    % another way: adding nuisance to the model
%     nuisance_idx = [];
%     if ~isempty (nuisance)
%         les = [les nuisance];
%         nuisance_idx = size(orig_les, 2)+1:size(les, 2);
%     end
%     
    
    
    if nOK < 1
        fprintf('Skipping SVM/SVR: no valid data\n');
    else
        % restructured by GY
        if nii_isBinary(beh1)
            % do 100 splits for voxelwise SVM, and 500 splits for ROI SVM
            if voxelwise_analysis
                [~, loadingMap{1}, ~, p] = nii_stat_svm_core(les, beh1, 100);
%                loadingMap{1}(nuisance_idx) = [];
            else
                [~, loadingMap{1}, ~, p] = nii_stat_svm_core(les, beh1, 500, 2); %do not specify thresholds: svm_core will select
%                loadingMap{1}(nuisance_idx) = [];
                
                out_name{1} = [statname '_' deblank(beh_name1) '_svm'];
                reportLoadingsSub (loadingMap{1}, les_names, deblank (beh_name1), p, 1);
            end
        else
            if voxelwise_analysis
                [~, loadingMap{1}, ~, p] = nii_stat_svr_core (les, beh1, deblank (beh_name1), 0);
%                loadingMap{1}(nuisance_idx) = [];
                
            else
                clipping_list = [0 1 -1];
                clipping_str = {'2tail' '1tailPOS' '1tailNEG'};
                for k = 1:length(clipping_list)
                    [~, loadingMap{k}, ~, p] = nii_stat_svr_core(les, beh1, deblank (beh_name1), clipping_list(k), 2); %compute regression
%                   loadingMap{k}(nuisance_idx) = [];
                    
                    out_name{k} = [statname '_' deblank(beh_name1) '_svr_' clipping_str{k}];
                    reportLoadingsSub (loadingMap{k}, les_names, deblank (beh_name1), p, 0);
                end
            end
        end
        if ~isempty (loadingMap{1}) % if analysis didn't work, loadingMap will be empty --GY
            if exist('roifname','var') && ~isempty(roifname)
                for k = 1:length(loadingMap) % length is either 1 for SVM and vox SVR, or 3 for ROI SVR
                    [threshMin, threshMax] = mc_corrected_threshold (loadingMap{k}, pThresh, mcorr_method);
                    unfolded_map = zeros (length (logicalMask), 1);
                    unfolded_map (logicalMask) = loadingMap{k};
                    thresh_map = threshold_map (unfolded_map, threshMin, threshMax);
                    suprathreshold = sum (thresh_map ~= 0 & ~isnan (thresh_map));
                    fprintf ('Thresholds: %.3g, %.3g. Range of feature weights: %.3g to %.3g (%d features pass the thresholod)\n', threshMin, threshMax, min(unfolded_map), max(unfolded_map), suprathreshold); 
                    if ~analyzing_connectome
                        nii_array2roi (unfolded_map, roifname, [out_name{k} '_unthreshZ.nii']);
                        if ~isempty (thresh_map)
                            nii_array2roi (thresh_map, roifname, [out_name{k} '_threshZ.nii']);
                        end
                    else
                        weight_matrix = zeros (nLabel, nLabel);
                        upper_triangle = logical (triu (ones (nLabel), 1));
                        weight_matrix (upper_triangle) = unfolded_map;
                        [~, atlas_name] = fileparts (roifname);
                        nii_save_nodz(atlas_name, weight_matrix, [out_name{k} '_unthreshZ.nodz'], logicalMask);
                        if ~isempty (thresh_map)
                            weight_matrix (upper_triangle) = thresh_map;
                            nii_save_nodz(atlas_name, weight_matrix, [out_name{k} '_threshZ.nodz'], logicalMask);
                        end
                    end
                end % for k = 1:length(loadingMap)
            end
            if voxelwise_analysis
                [threshMin, threshMax] = mc_corrected_threshold (loadingMap{1}, pThresh, mcorr_method);
                thresh_map = threshold_map (loadingMap{1}, threshMin, threshMax);
                out_name = [statname '_' deblank(beh_name1) '_svr'];
                save_voxelwise_loadings (loadingMap{1}, logicalMask, hdr, [out_name '_unthreshZ']);
                if ~isempty (thresh_map)
                    save_voxelwise_loadings (thresh_map, logicalMask, hdr, [out_name '_threshZ']);
                end 
            end
        end
        % /GY
    end
end
diary off %stop logging text
cd .. %leave the folder created by chDirSub
%end nii_stat_svm() LOCAL FUNCTIONS FOLLOW

function [fnm, nOK] = tabFileSub(les,beh1, beh_name1,  les_names, subj_data)  
if size(les,1) ~= size(beh1,1)
    error('nii_stat_svm confused');
end
fnm = [beh_name1  '.tab'];
fid = fopen(fnm, 'w');
n_subj = size(les,1);
fprintf(fid,'filename\t');
for j = 1:numel(les_names)
     fprintf(fid,'%s\t', les_names{j}); 
end
fprintf(fid,'%s\t', beh_name1);
fprintf(fid,'\n');
nOK = 0;
for i = 1:n_subj
    if  ~isfinite(std(les(i,:)))
        fprintf('%s WARNING: Skipping %s due to bogus data (NaN)\n', mfilename, subj_data{i}.filename);
    else
         if (std(les(i,:)) == 0) 
            fprintf('%s WARNING: No variability in imaging data for %s (all regions have an intensity of %g)\n', mfilename, subj_data{i}.filename, les(i,1));
        end
        if ~isempty(subj_data)
           fprintf(fid,'%s\t',subj_data{i}.filename); 
        else
            fprintf(fid,'%s\t',num2str(i));
        end
        for j = 1:numel(les_names)
             fprintf(fid,'%g\t',les(i, j));
        end
        fprintf(fid,'%g\t',beh1(i));
        fprintf(fid,'\n');
        nOK = nOK + 1;
    end
end
fclose(fid);
%end tabFileSub()

function chDirSub(statname)
datetime=datestr(now);
datetime=strrep(datetime,':',''); %Replace colon with underscore
datetime=strrep(datetime,'-','');%Replace minus sign with underscore
datetime=strrep(datetime,' ','_');%Replace space with underscore
newdir = [datetime statname];
mkdir(newdir);
cd(newdir);
%chDirSub()

% function saveNodzSub(roiname, matvals, nodzname) 
% if min(matvals(:)) == max(matvals(:)), fprintf(' No variability, will not create %s\n', nodzname); end;
% [kROI, kROINumbers, ROIIndex] = nii_roi_list(roiname, false);
% if ROIIndex < 1, return; end; %unable to find ROI
% str = nii_roi2mm (ROIIndex);
% fileID = fopen(nodzname,'w');
% fprintf(fileID, str);
% fprintf(fileID, '#ENDNODE\n');
% fclose(fileID);
% dlmwrite(nodzname,matvals,'delimiter','\t','-append')
% %saveNodzSub

function reportLoadingsSub (featureLoadings, les_names, beh_name, p, classification)
if p > 0.3
    fprintf ('Features are not reported because SVM/SVR accuracy is too poor\n');
    return;
end
if classification
    pos_str = 'class 0';
    neg_str = 'class 1';
else
    pos_str = 'higher score';
    neg_str = 'lower score';
end
threshZ = 1; % a rather arbitary threshold of "significance"
passed_thresh = abs(featureLoadings) > threshZ;
fprintf ('%s: %d features have weights greater than %d or less than -%d\n', beh_name, sum (passed_thresh), threshZ, threshZ);
featureLoadings = featureLoadings (passed_thresh);
les_names = les_names (passed_thresh);
[featureLoadings, sort_idx] = sort (featureLoadings, 'descend');
n_reported_pos = min (10, sum (featureLoadings > 0));
n_reported_neg = min (10, sum (featureLoadings < 0));
if n_reported_pos > 0
    if n_reported_pos == 10
        fprintf ('TOP 10 POSITIVE FEATURES (higher value -> %s):\n', pos_str);
    else
        fprintf ('POSITIVE FEATURES (higher value -> %s):\n', pos_str);
    end
    for i = 1:n_reported_pos
        fprintf ('%s (%g)\n', les_names{sort_idx(i)}, featureLoadings(i));
    end
end
if n_reported_neg > 0
    if n_reported_neg == 10
        fprintf ('TOP 10 NEGATIVE FEATURES (higher value -> %s):\n', neg_str);
    else
        fprintf ('NEGATIVE FEATURES (higher value -> %s):\n', neg_str);
    end
    for i = length(featureLoadings):-1:length(featureLoadings)-n_reported_neg+1        
        fprintf ('%s (%g)\n', les_names{sort_idx(i)}, featureLoadings(i));
    end
end
%reportLoadingsSub

function save_voxelwise_loadings (loadingMap, logicalMask, hdr, statName)
unfolded_map = zeros (hdr.dim);
unfolded_map (logicalMask) = loadingMap;
hdr.fname = [statName '.nii'];
hdr.pinfo = [1;0;0];
hdr.private.dat.scl_slope = 1;
hdr.private.dat.scl_inter = 0;
hdr.private.dat.dtype = 'FLOAT32-LE';%'INT16-LE', 'FLOAT32-LE';
hdr.dt    =[16,0]; %4= 16-bit integer; 16 =32-bit real datatype
spm_write_vol(hdr,unfolded_map);
% end save_voxelwise_loadings

function [threshMin, threshMax] = mc_corrected_threshold (zmap, pThresh, mcorr_method)
zmap = zmap (~isnan (zmap));
p2z = @(p) -sqrt(2) * erfcinv(p*2);
if mcorr_method == 0 % Bonferroni
    bonferroniP = pThresh / length (zmap);
    threshMin = -abs(p2z(bonferroniP));
    threshMax = abs(p2z(bonferroniP));
else % FDR
    p = spm_Ncdf(zmap);
    [~, crit_p, ~]=fdr_bh(p,pThresh,'pdep');
    threshMin = p2z(crit_p);
    p = spm_Ncdf(1-zmap); 
    [~, crit_p, ~]=fdr_bh(p,pThresh,'pdep');
    threshMax = -p2z(crit_p);
end

function thresh_map = threshold_map (unfolded_map, threshMin, threshMax)
thresh_map = unfolded_map;
thresh_map (thresh_map < 0 & thresh_map > threshMin) = 0;
thresh_map (thresh_map > 0 & thresh_map < threshMax) = 0;
if sum (thresh_map ~= 0 & ~isnan (thresh_map)) == 0 
    thresh_map = [];
end
