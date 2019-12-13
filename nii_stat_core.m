function nii_stat_core(les,beh, beh_names,hdr, kPcrit, kNumRandPerm, logicalMask,statname, roi_names, hdrTFCE, nuisance)
%Generates statistical tests. Called by the wrappers nii_stat_mat and nii_stat_val
% les    : map of lesions/regions one row per participant, one column per voxel
% beh    : matrix of behavior, one row per participant, one column per condition
% beh_names : text labels for behavioral data, one cell per condition
% hdr    : (optional) structure for header of NIfTI image data
% kPkrit : one-tailed statistical threshold
% kNumRandPerm : number of permutations to control for multiple comparisons
%                 set to a large value (2000-8000 to find peak derived threshold
%                 0 = Bonferroni Threshold (very fast, very conservative)
%                 -1 = False Discovery Rate (fast, different inference, depends on proportion of signal
%                 set to extreme negative value (e.g. -5000) for
%               Ged Ridgeway notes you can estimate confidence intervals for peak permutation, e.g. 8000 permutes with p < 0.05 [~, ci] = binofit(0.05*8000, 8000)
% logicalMask: ROIs/voxels are marked with 1 if they are to be analyzed, and with 0 if they are to be ignored (GY)
% kOnlyAnalyzeRegionsDamagedInAtleastNSubjects : (optional) only influences binomial imaging data
%             will only examine voxels where at least this many individuals
%             have a "1". This limits the number of tests to regions that
%             are frequently injured. Perhaps 20% of number of participants
%             is reasonable
% statname : (optional) results saved to disk will include this name. If not provided this becomes 'anonymous'
% roi_names : (optional) if one per column of les, then a labeled table is produced
% hdrTFCE   : (optional, voxelwise analyses only) - dimensions of volume
% nuisance: (optional) matrix of nuisance variables
%Examples:
%     five0five1 = [0 0 0 0 0 1 1 1 1 1]';
%     ascend1to10 = [1:10]';
%     mostly0mostly1 = [0 0 1 0 0 1 1 0 1 1]';
%     mostlyascend = [1 2 3 4 5 3 4 5 6 7]';
%     nii_stat_core(five0five1,mostly0mostly1,{'ZeroImagesTendToScoreZero_OneImagesTendToScoreOne'},[],0.05,0) %Liebermeister VLSM, z=1.750184
%     nii_stat_core(five0five1,mostlyascend,{'ZeroImagesScoreLowerThanOneImages'},[],0.05,0) %t-test VLSM, z=1.747699
%     nii_stat_core(ascend1to10,mostly0mostly1,{'darkerImagestendToBeInGroup0_brighterImagesTendToBeGroup1'},[],0.05,0); %t-test VBM,  reports z=1.548496
%     nii_stat_core(ascend1to10,mostlyascend,{'positiveCorrelationBetweenBrightnessAndScores'},[],0.05,0); %regression: reports z=3.595802
%      v = [6 5 4 3 2 1 2 3 4 5]';
%     factors = [[1:10]' v];
%     img = ([1:10]' + v);
%     noise = rand(10,1)*0.1;
%     imgnoise = img+noise;
%     nii_stat_core(imgnoise,factors,{'ascending'; 'v-shape'},[],0.05,0); %multiple independent regression
%     nii_stat_core(imgnoise,factors,{'ascending'; 'v-shape'},[],0.05,-1000); %freedman lane regression treating other variables as nuisance regressors

numFactors = size(beh,2); %how many behavioral variables are being tested?
numObs = size(beh,1); %number of observations
if ~exist('kOnlyAnalyzeRegionsDamagedInAtleastNSubjects','var')
    kOnlyAnalyzeRegionsDamagedInAtleastNSubjects = 1;
end
if ~exist('kPcrit','var')
    kPcrit = 0.05;
end
if ~exist('kNumRandPerm','var')
    kNumRandPerm = 0;  %Bonferroni correction
end
if ~exist('roi_names','var')
    roi_names = '';
end
if ~exist('statname','var')
    statname = 'anonymous';
end
if ~exist('hdrTFCE','var')
    hdrTFCE = [];
end
if ~exist('nuisance','var')
    nuisance = [];
end

if size(les,1) ~= size(beh,1)
    error('Error: les and beh must have the same number of rows (one image and one set of behavioral data per participant)');
end
if numFactors < 1
   error('There must be at least one behavioral factor (beh)');
end
if length(beh_names) ~= numFactors
   error('There must be one label (beh_names) for each column of behavior (beh)');
end
if size(les,1) < 3
    error('Error: data must have at least 3 rows (observations). Perhaps inputs need to be transposed?');
end
[isBinomialBehav, beh] = ifBinomialForce01Sub(beh,true); %is behavioral data binary?
[isBinomialLes, les] = ifBinomialForce01Sub(les);
if ((isBinomialLes ||isBinomialBehav) && (numel(hdrTFCE) == 3))
	fprintf('Both behavior and images must be continuous for TFCE.');
    hdrTFCE = [];
end

%%% GY: select good ROIs/voxels here
numROIorVox = size(les,2); %number of regions or voxels, regardless of being analyzed or ignored

good_idx = find (logicalMask);
%added by GY
if isempty (good_idx)
    error ('No ROIs or voxels pass selection criteria. Exiting...');
end
if size(les,2) == numel(logicalMask) %CR addded this conditional
    les = les (:, good_idx); %squeeze data to only examine critical voxels
end
    
%compute statistics
%The following code is no longer used: already encoded in logicalMask
%if exist('voxMask','var') && ~isempty(voxMask) %convert good_idx to unpacked voxels
%    %voxMask = [1 0 1; 0 1 0; 1 0 1]; good_idx = [1 3]; %<- illustrate logic, convert address from [1 3] to [1 5]
%    voxPos = find(voxMask ~= 0);
%    good_idx = voxPos(good_idx);
%end
if isempty(roi_names) %voxelwise analysis
    sumImg = zeros(hdr.dim(1), hdr.dim(2), hdr.dim(3));
    sumImg(good_idx(:)) = sum (les, 1);
    saveSumMapSub(hdr, sumImg,statname);%, voxMask); %create image showing sum of values
else
    sumImg = zeros(numROIorVox,1);
    sumImg(good_idx(:)) = sum (les, 1);
    saveSumMapROI(hdr,sumImg,statname); %create image showing sum of values
    %saveSumMapROI(hdr,sum (les, 1)/size(les,1),statname); %create image showing sum of values
end
numVox = size(les,2); %number of voxels/regions
z = zeros(numVox,numFactors); %pre-allocate observed z-scores
threshMax = zeros(numFactors,1); %pre-allocate statistical threhsold
threshMin = threshMax;  %pre-allocate statistical threhsold
startTime = tic;
if (kNumRandPerm < -500) ||  (kNumRandPerm > 0) %if will estimate thresholds with permutation
    rand ('seed', 6666);  %#ok<RAND> %for newer versions: rng(6666); set random number seed - call this to ensure precise permutations between runs
    %fprintf('NOTE: random number generator seed specified. Multiple runs will generate identical values.\n');
end
if (kNumRandPerm < -1) && (size(beh,2) <= 1) %special case: nuisance regressors which requires multiple behaviors
    error('Error: only one behavioral variable provided: unable to control for nuisance regressors (solution: either provide more behavioral data or specify positive number of permutations');
elseif (kNumRandPerm < -1) && (size(beh,2) > 1) %special case: nuisance regressors with permutations.
    X1 = [beh, ones(numObs,1)];%add constant column for intercept
    global global_flContrast
    if isempty(global_flContrast)
        fprintf('Computing regression with nuisance variables for %d regions/voxels, analyzing %d behavioral variables.\n',length(good_idx),size(beh,2));
        for i = 1: numFactors
            contrasts = zeros(numFactors+1,1); %+1 since last column is constant intercept
            contrasts(i) = 1;
            [z(:,i), threshMin(i), threshMax(i)] = glm_perm_flSub(les, X1, contrasts, abs(kNumRandPerm),kPcrit, good_idx, hdrTFCE);
        end;
    else
        fprintf('Computing freedman-lane regression with custom contrast for %d regions/voxels, analyzing %d behavioral variables.\n',length(good_idx),size(beh,2));
        contrasts = zeros(numFactors+1,1); %+1 since last column is constant intercept
        contrasts(1:length(global_flContrast)) = global_flContrast;
        if numel(global_flContrast) ~= numel(contrasts)
            fprintf(' Noted: contrast zero padded %s\n',mat2str(contrasts'));
        end
        %we will only analyse a single contrast
        numFactors = 1;
        z = zeros(numVox,numFactors); %pre-allocate observed z-scores
        threshMax = zeros(numFactors,1); %pre-allocate statistical threhsold
        threshMin = threshMax;  %pre-allocate statistical threhsold
        i = 1;
        [z(:,i), threshMin(i), threshMax(i)] = glm_perm_flSub(les, X1, contrasts, abs(kNumRandPerm),kPcrit, good_idx, hdrTFCE);
        %give the contrast a meaningful name
        s = mat2str(global_flContrast);
        s = strrep(s, ' ', '_');
        s = strrep(s, ']', '');
        s = strrep(s, '[', '');
        beh_names = [];
        beh_names{1} = ['CustomContrast_' s];
    end;

elseif isBinomialBehav && isBinomialLes %binomial data
    fprintf('Computing Liebermeister measures for %d regions/voxels with %d behavioral variables(positive Z when 0 voxels have behavior 0 and 1 voxels have behavior 1).\n',length(good_idx),size(beh,2));
    for i = 1: size(beh,2)
        [z(:,i), threshMin(i), threshMax(i)] = lieber_permSub(les,beh(:,i), kNumRandPerm, kPcrit, hdrTFCE);
    end;
else %behavior and/or lesions is continuous
	fprintf('Computing glm (pooled-variance t-test, linear regression) for %d regions/voxels with analyzing %d behavioral variables (positive Z when increased image brightness correlates with increased behavioral score).\n',length(good_idx),size(beh,2));
    
    % GY, Aug 2019: nuisance regressors
    if isempty (nuisance)
        % the good old way
        for i = 1: size(beh,2)
            [z(:,i), threshMin(i), threshMax(i)] = glm_permSub(les,beh(:,i), kNumRandPerm, kPcrit, good_idx, hdrTFCE);
        end;
    else
        % the new way, with nuisance regressors
        for i = 1: size(beh,2)
            [z(:,i), threshMin(i), threshMax(i)] = glm_permSub_nuisance(les,beh(:,i), nuisance, kNumRandPerm, kPcrit, good_idx, hdrTFCE);
        end;    
    end
end
%global global_powerMap %option to save parameters for power analysis
%if ~isempty(global_powerMap) && global_powerMap
    savePowerSub(les,beh(:,i),good_idx, hdr, roi_names, beh_names);
%end
%next: report thresholds
if (kNumRandPerm == -1) || (kNumRandPerm == -2) %report thresholds using FDR correction
    for i = 1:numFactors
        p2z = @(p) -sqrt(2) * erfcinv(p*2);
        p = spm_Ncdf(z(:,i));
        [~, crit_p, ~]=fdr_bh(p,kPcrit,'dep');
        threshMin(i) = p2z(crit_p);
        p = spm_Ncdf(1-z(:,i)); %convert z-scores to probabilities
        [~, crit_p, ~]=fdr_bh(p,kPcrit,'dep');
        threshMax(i) = -p2z(crit_p);
        fprintf('q=%.4f FDR correction for %s with %d tests is z<%.8f, z>%.8f\n',kPcrit,deblank(beh_names{i}),length(good_idx),threshMin(i),threshMax(i));
    end
elseif kNumRandPerm == 0 %report thresholds using bonferroni correction
    nTests = length(good_idx); %familywise: chance of error in single volume * size(beh,2);
    bonferroniP = kPcrit / nTests;
    p2z = @(p) -sqrt(2) * erfcinv(p*2);
    bonferroniZ = abs(p2z(bonferroniP));
    fprintf('p<%f Bonferroni correction for %d tests is z=%.8f, p<%.12f (one tailed, we assume lesions always impair performance)\n',kPcrit,nTests,bonferroniZ,bonferroniP);
    threshMin(:) = -bonferroniZ;
    threshMax(:) = bonferroniZ;
else %threshold using pre-computed permutations
    fprintf('%d permutations required %.4f seconds\n',abs(kNumRandPerm),toc(startTime)); %/3600 for hours
    fprintf('Thresholds are one tailed (we predict injured tissue will only cause poorer performance not better performance\n');
    for i = 1:numFactors
        fprintf('p<%.3f permutation correction for %s is z<%.5f z>%.5f\n',kPcrit,deblank(beh_names{i}),threshMin(i), threshMax(i));
    end;
end %else report permutation thresholds
%next: report results
if isempty(roi_names) %voxelwise
	reportResultsVoxel(z,[],threshMin,threshMax,good_idx,beh_names,hdr,statname);%, voxMask);
elseif numROIorVox ~= size(roi_names,1)
    reportResultsMatrix(z,[],threshMin,threshMax,good_idx,beh_names,roi_names,numROIorVox,hdr,statname, logicalMask);
else %if not voxelwise or matrix: must be region of interest analysis
	reportResultsROI(z,[],threshMin,threshMax,good_idx,beh_names,roi_names,hdr,statname);
end
%diary off
%cd .. %leave the folder created by chDirSub
%end nii_stat_core()

%%%%% SUBFUNCTIONS FOLLOW %%%%%%%


function reportResultsMatrix(z,c,threshMin,threshMax,good_idx,beh_names,roi_names, numROIorVox,hdr,statname, logicalMask)
%handshake problem, http://en.wikipedia.org/wiki/Triangular_number http://math.fau.edu/richman/mla/triangle.htm
nLabel = size(roi_names,1);
nTri = nLabel*(nLabel-1)*0.5;
if nTri ~= numROIorVox
    error('reportResultsMatrix error: expecting a matrix triangle');
end
%reflect from triangle matrix to full matrix
TriBin = triu(ones(nLabel,nLabel), 1);
Idx = 1:numROIorVox;
FullMat = zeros(nLabel,nLabel);
FullMat(TriBin==1) = Idx;
stat.label  = roi_names;
if isempty(hdr) || length(hdr.fname) < 1
    stat.roiname = '';
else
	[~,stat.roiname] = fileparts(hdr.fname);
end
for i = 1:length(beh_names)
    disp ('___');
    z_column = z(:, i);
    signif_idx = union(find(z_column(:) < threshMin(i)), find(z_column(:) > threshMax(i)) );
    fprintf('%s z=%f..%f, %d regions survive threshold\n', deblank(beh_names{i}), min(z_column(:)),max(z_column(:)),length (signif_idx) );
    stat.threshZ = zeros(nLabel,nLabel);
    stat.threshR = zeros(nLabel,nLabel);
    for j = 1:length (signif_idx)
        %area_name = roi_names{good_idx(signif_idx(j))};
        [xRow,yCol] = find(FullMat == good_idx(signif_idx(j)),1) ;
        area_name = [roi_names{xRow} '*' roi_names{yCol}];
        if ~isempty(c)
            C_column = c(:, i);
            stat.threshR(xRow,yCol) = C_column(signif_idx(j));
            fprintf ('%s (r = %f) z=%.7f\n', strtrim (area_name), C_column(signif_idx(j)),z_column(signif_idx(j)));
        else
            fprintf ('%s z=%.7f\n', strtrim (area_name), z_column(signif_idx(j)));
        end
        stat.threshZ(xRow,yCol) = z_column(signif_idx(j));
    end
    %unthresholded data
    stat.unthreshZ = zeros(nLabel,nLabel);
    for j = 1:length (z_column)
        [xRow,yCol] = find(FullMat == good_idx(j),1) ;
        stat.unthreshZ(xRow,yCol) = z_column(j);
    end
    %save data
    stat.logicalMask = logicalMask;
    matname = sprintf ('%s_%s.mat',statname, deblank(beh_names{i}));
    save(matname,'-struct', 'stat');
    %save node masks
    nodzname = sprintf ('%s_%sZ.nodz',statname, deblank(beh_names{i}));
    nii_save_nodz(stat.roiname, stat.threshZ, nodzname, logicalMask);
    nodzname = sprintf ('%s_%s_unthreshZ.nodz',statname, deblank(beh_names{i}));
    nii_save_nodz(stat.roiname, stat.unthreshZ, nodzname, logicalMask); 
    if ~isempty(c) %if we have correlation values
        nodzname = sprintf ('%s_%sR.nodz',statname, deblank(beh_names{i}));
        nii_save_nodz(stat.roiname, stat.threshR, nodzname, logicalMask);
    end;
end
%end reportResultsMatrix()

% function saveNodzSub(roiname, matvals, nodzname, good_idx) 
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

function reportResultsVoxel(z,c,threshMin,threshMax,good_idx,beh_names,hdr, statname) %,voxMask)
for i = 1:length(beh_names) %or size(beh,2)
    thresh_idx = union(find(z(:,i) < threshMin(i)), find(z(:,i) > threshMax(i)) );
    fprintf('%s z=%f..%f, %d voxels survive threshold\n', deblank(beh_names{i}),min(z(:,i)),max(z(:,i)),length(thresh_idx) );
    if ~isempty(hdr)
        zMask = z(:,i);
        %if isempty(voxMask)
            statImg = zeros(hdr.dim(1), hdr.dim(2), hdr.dim(3));
            statImg(good_idx(:)) = zMask;
        %else
        %    statImg = unmaskSub2(zMask, good_idx, voxMask);
        %end
        %statImg = reshape(statImg, hdr.dim(1), hdr.dim(2), hdr.dim(3));
        %statImg = reshape(statImg,hdr.dim(1),hdr.dim(2),hdr.dim(3));
        saveStatMapSub(hdr,statImg,[sprintf('Z%s', statname), deblank(beh_names{i})], threshMin(i),threshMax(i) );
        if ~isempty(thresh_idx)
            zMask = z(:,i);
            subThresh_idx = setdiff (1:size(zMask, 1), thresh_idx);
            zMask(subThresh_idx) = 0;
            %statImg = zeros(hdr.dim(1)*hdr.dim(2)*hdr.dim(3),1);
            %if isempty(voxMask)
                statImg = zeros(hdr.dim(1), hdr.dim(2), hdr.dim(3));
                statImg(good_idx(:)) = zMask;
            %else
            %    statImg = unmaskSub2(zMask, good_idx, voxMask);
            %end
            %statImg = reshape(statImg,hdr.dim(1),hdr.dim(2),hdr.dim(3));
            saveStatMapSub(hdr,statImg,[sprintf('threshZ%s', statname), deblank(beh_names{i})], threshMin(i),threshMax(i) );
            if ~isempty(c)
                    cMask = c(:,i);
                    cMask(subThresh_idx) = 0;
                    %statImg = zeros(hdr.dim(1)*hdr.dim(2)*hdr.dim(3),1);
                    %if isempty(voxMask)
                        statImg = zeros(hdr.dim(1), hdr.dim(2), hdr.dim(3));
                        statImg(good_idx(:)) = cMask;
                    %else
                    %    statImg = unmaskSub2(cMask, good_idx, voxMask);
                    %end
                    %statImg = reshape(statImg,hdr.dim(1),hdr.dim(2),hdr.dim(3));
                    saveStatMapSub(hdr,statImg,['threshC',statname, deblank(beh_names{i})], threshMin(i),threshMax(i) );
            end %if correlations
        end; %if some voxels survive threshold
    end; %if hdr exists
end %for each behavior
%end reportResultsVoxel()

function saveSumMapROI(hdr,sumImg,statname)
if isempty(hdr), return; end;
img = spm_read_vols(hdr);
vimg = zeros(size(img));
for i = 1: numel(sumImg)
    vimg(img == i) = sumImg(i);
end
if size(img(:)) ~= size(sumImg(:))
    fprintf('%s did not create a sum image, dimensions do not match template image %s\n',mfilename,fname);
    return
end
%sumImg = reshape(sumImg,size(img));
hdr.fname = [deblank(statname) 'sum.nii'];
hdr.pinfo = [1;0;0];
hdr.private.dat.scl_slope = 1;
hdr.private.dat.scl_inter = 0;
hdr.private.dat.dtype = 'FLOAT32-LE';%'INT16-LE', 'FLOAT32-LE';
hdr.dt    =[16,0]; %4= 16-bit integer; 16 =32-bit real datatype
spm_write_vol(hdr,vimg);
%end saveSumMapROI


function reportResultsROI(z,c,threshMin,threshMax,good_idx,beh_names,roi_names,hdr,statname)
for i = 1:length(beh_names)
    disp ('___');
    z_column = z(:, i);
    signif_idx = union(find(z_column(:) < threshMin(i)), find(z_column(:) > threshMax(i)) );
    fprintf('%s z=%f..%f, %d regions survive threshold\n', deblank(beh_names{i}), min(z_column(:)),max(z_column(:)),length (signif_idx) );
    for j = 1:length (signif_idx)
        area_name = roi_names{good_idx(signif_idx(j))};
        if ~isempty(c)
            C_column = c(:, i);
            fprintf ('%s (r = %f) z=%.7f\n', strtrim (area_name), C_column(signif_idx(j)),z_column(signif_idx(j)));
            %msg = sprintf ('%s (r = %f) z=%.7f', strtrim (area_name), C_column(signif_idx(j)),z_column(signif_idx(j)));
        else
            fprintf ('%s z=%.7f\n', strtrim (area_name), z_column(signif_idx(j)));
            %msg = sprintf ('%s z=%.7f', strtrim (area_name),z_column(signif_idx(j)));
        end
        %disp (msg);
    end
end
%next: (optional) write statistical maps as NIfTI format images
if isempty(hdr) || length(hdr.fname) < 1
    fprintf('No statistcal maps saved\n');
    return
elseif exist(hdr.fname,'file') ~= 2
    fprintf('No statistcal maps saved: unable to find image named %s\n',hdr.fname);
    return
end
img = spm_read_vols(hdr); %load ROI
%create new headers with same dimensions as ROI
hdrC = hdr;
hdrC.dt(1) = 16; %make 32 bit real
hdrC.private.dat.dtype = 'FLOAT32-LE';
hdrC.private.dat.scl_slope = 1;
hdrC.private.dat.scl_inter = 0;
hdrT = hdrC;
%Not sure how to set the NIfTI intent to Z-score with SPM....
%hdrT.private.intent.code = 'ZSCORE';
for i = 1:length(beh_names) %disp (beh_names{i});
    %P_column = t (:, i);
    %signif_idx = find (P_column < thresh(i)); % areas with significant correlation
    %signif_idx = union(find(t(i,:) < threshMin(i)), find(t(i,:) > threshMax(i)) );
    z_column = z(:, i);
    hdrT.fname = sprintf ('Z%s_%s.nii',statname, deblank(beh_names{i}));
    imgOutT = zeros(size(img));
    for j = 1:length (good_idx)
        imgOutT(img == good_idx(j)) = z_column(j);
    end %for all significant regions
    spm_write_vol(hdrT,imgOutT);
    %see if any voxels survive threshold
    signif_idx = union(find(z_column(:) < threshMin(i)), find(z_column(:) > threshMax(i)) );
    if ~isempty (signif_idx)
        hdrT.fname = sprintf ('threshZ%s_%s.nii',statname, deblank(beh_names{i}));
        imgOutT = zeros(size(img));
        for j = 1:length (signif_idx)
            zscore = z_column(signif_idx(j));
            imgOutT(img == good_idx(signif_idx(j))) = zscore;
        end %for all significant regions
        spm_write_vol(hdrT,imgOutT);
         if ~isempty(c)
            C_column = c(:, i);
            hdrC.fname = sprintf ('threshC_%s_%s.nii',statname, deblank(beh_names{i}));
            imgOutC = zeros(size(img));
            for j = 1:length (signif_idx)
                cV = C_column(signif_idx(j));
                imgOutC(img == good_idx(signif_idx(j))) = cV;
            end %for all significant regions
            spm_write_vol(hdrC,imgOutC);
        end
    end %if any significant regions
end %for all behavioral variables
%end reportResultsROI()

function [isBinomial, behav] = ifBinomialForce01Sub(behav,warnIfMixed)
%detect if data is binomial, if so force to be 0 and 1
mn = min(behav(:));
mx = max(behav(:));
if mn == mx
   fprintf('Warning: data has no variability. Not useful for analyses.\n');
end
nMin = sum(behav(:)==mn);
nMax = sum(behav(:)==mx);
if (nMin+nMax) ~= length(behav(:))
    isBinomial = false;
    if (sum(isnan(behav(:)))+nMin+nMax) == length(behav(:))
        fprintf('Warning: data would be binomial if not-a-numbers were reset to zero');
    end
    if exist('warnIfMixed','var') && (warnIfMixed == true) && (size(behav,2) > 1)
        for i = 1 : size(behav,2)
            mn = min(behav(:,i));
            mx = max(behav(:,i));
            nMin = sum(behav(:,i)==mn);
            nMax = sum(behav(:,i)==mx);
            if (nMin+nMax) == size(behav,1)
                fprintf('Warning: some columns are continuous and others are binomial (will analyze all variables as continuous, consider conducting separate analyses)\n');
                return;
            end %if column is binomial
        end %for each column
    end %if warning
    return;
end
isBinomial = true;
behav(behav(:)==mn) = 0;
behav(behav(:)==mx) = 1;
%end ifBinomialForce01Sub()

function saveSumMapSub(hdr,sumImg,statname) %, voxMask)
if isempty(hdr), return; end;
%sumImg = unmaskSub(sumImg, voxMask);
img = zeros(hdr.dim(1),hdr.dim(2),hdr.dim(3));
if size(img(:)) ~= size(sumImg(:))
    fprintf('%s did not create a sum image, dimensions do not match template image %s\n',mfilename,fname);
    return
end
sumImg = reshape(sumImg,size(img));
hdr.fname = [deblank(statname) 'sum.nii'];
hdr.pinfo = [1;0;0];
hdr.private.dat.scl_slope = 1;
hdr.private.dat.scl_inter = 0;
hdr.private.dat.dtype = 'FLOAT32-LE';%'INT16-LE', 'FLOAT32-LE';
hdr.dt    =[16,0]; %4= 16-bit integer; 16 =32-bit real datatype
spm_write_vol(hdr,sumImg);
%end saveSumMapSub()

% function img = unmaskSub2(img, good_idx, voxMask) %uncompress packed 1D vector to sparse 3D volume
% if isempty(voxMask)
%     return
% end
% old = img;
% img = zeros(sum(voxMask(:) ~= 0),1);
% img(:) = nan;
% img(good_idx(:)) = old;
% img = unmaskSub(img, voxMask);
% %end unmaskSub2()
%
% function img = unmaskSub(img, voxMask) %uncompress packed 1D vector to sparse 3D volume
% if isempty(voxMask)
%     return
% end
% old = img;
% img = zeros(size(voxMask));
% img(:) = nan;
% img(voxMask == 1) = old;
% %end unmaskSub()

function saveStatMapSub(hdr,statImg, statName,  minThreshold, maxThreshold)
%save map as NIfTI format image
img = zeros(hdr.dim(1),hdr.dim(2),hdr.dim(3));
if size(img(:)) ~= size(statImg(:))
    fprintf('%s did not create a sum image, dimensions do not match template image %s\n',mfilename,fname);
    return
end
%statImg = unmaskSub(statImg, voxMask);
statImg = reshape(statImg,size(img));
hdr.fname = [statName '.nii'];
hdr.descrip = sprintf('z U<%.3f U>%.3f', minThreshold, maxThreshold);
hdr.pinfo = [1;0;0];
hdr.private.dat.scl_slope = 1;
hdr.private.dat.scl_inter = 0;
hdr.private.dat.dtype = 'FLOAT32-LE';%'INT16-LE', 'FLOAT32-LE';
hdr.dt    =[16,0]; %4= 16-bit integer; 16 =32-bit real datatype
spm_write_vol(hdr,statImg);
%end saveStatMapSub()

function [uncZ, threshMin, threshMax] = lieber_permSub(Y, X, nPerms, kPcrit, hdrTFCE)
%returns uncorrected z-score for all voxels in Y given single column predictor X
% Y: volume data
% X: single column predictor
% nPerms: Number of permutations to compute
% kPcrit: Critical threshold
%Liebermeister measure: both Y and X are binomial
%Example
% [z] = lieber_perm([1 1 1 1 0 0 0 0; 0 1 0 0 1 1 1 1]',[1 1 1 1 0 0 1 0]')
% [z, mn, mx] = lieber_perm([1 1 1 1 0 0 0 0; 0 1 0 0 1 1 1 1]',[1 1 1 1 0 0 1 0]',2000,0.05)
%inspired by Ged Ridgway's glm_perm_flz
if ~exist('nPerms','var')
    nPerms = 0;
end
if ~exist('kPcrit','var')
    kPcrit = 0.05;
end
if ~exist('hdrTFCE','var')
    hdrTFCE = [];
end
if numel(hdrTFCE) == 3
	fprintf('WARNING: Liebermeister measure does not support TFCE (yet)\n');
end
[n p] = size(X); %rows is number of observations, columns should be 1
if (p ~= 1), error('glm_t is for one column of X at a time'); end;
tic
uncZ = lieberSub(Y,X);
if nPerms < 2
    threshMin = -Inf;
    threshMax = Inf;
    return
end
peak = nan(nPerms,1);
nadir = nan(nPerms,1);
peak(1) = max(uncZ);
nadir(1) = min(uncZ);
for p = 2:nPerms
    Xp  = X(randperm(n), :);
    Zp = lieberSub(Y,Xp);
    peak(p) = max(Zp);
    nadir(p) = min(Zp);
end
threshMin = permThreshLowSub (nadir, kPcrit);
threshMax = permThreshHighSub (peak, kPcrit);
%end lieber_permSub()

function z = lieberSub (groupVector, obsVector)
%Liebermeister measure
%[z] = lieber([1 1 1 1 0 0 0 0; 0 1 0 0 1 1 1 1]',[1 1 1 1 0 0 1 0]')
%[z] = lieber2([1 1 1 1 0 0 0 0]',[0 1 0 0 1 1 1 1]')
if size(obsVector,2) ~= 1
    error('Error: obsVector must be a single column\n');
end
M = size(obsVector,1); %number of observations
if (M ~= size(groupVector,1))
    error('%s error: group and observations matrices must have same number of rows\n',mfilename);
end
if M < 3
    error('Unabled to test voxels with less than 3 observations');
end
N = sum(obsVector(:)); %behavior = 1
a = sum(bsxfun(@min,groupVector,obsVector));
K = sum(groupVector,1);
pL = fexact(a,M,K,N,'test','l','tail','l');
pR = fexact(a,M,K,N,'test','l','tail','r');
indexL = find(pL(:) < pR(:));
p2z = @(p) -sqrt(2) * erfcinv(p*2);
z = pR; %vector Z now probability of right tail
z(indexL) = pL(indexL); %vector Z now probability of extreme tail
z = -p2z(z); %convert probabilities to Z-scores, 1/2014: minus so p0.05 is positive 1.64
z(indexL) = -z(indexL); %negative correlations have negative Z scores
%end lieberSub()

function [uncZ, threshMin, threshMax] = glm_permSub(Y, X, nPerms, kPcrit, good_idx, hdrTFCE)
%returns uncorrected z-score for all voxels in Y given single column predictor X
% Y: volume data
% X: single column predictor
%if either X or Y is binomial, results are pooled variance t-test, else correlation coefficient
% nPerms: Number of permutations to compute
% kPcrit: Critical threshold
%Example
% glm_t([1 1 0 0 0 1; 0 0 1 1 1 0]',[1 2 3 4 5 6]') %pooled variance t-test
%
%inspired by Ged Ridgway's glm_perm_flz
% save('glm_permSub'); %<-troublshoot, e.g. load('glm_permSub.mat'); glm_perm(Y, X, nPerms, kPcrit, good_idx, hdrTFCE)
if ~exist('nPerms','var')
    nPerms = 0;
end
if ~exist('kPcrit','var')
    kPcrit = 0.05;
end
% Basics and reusable components
[n f] = size(X); %rows=observations, columns=factors
[nY v] = size(Y); %#ok<NASGU> v is number of tests/voxels
if (f ~= 1), error('glm_permSub is for one column of X at a time (transpose?)'); end;
if nY ~= n, error('glm_permSub X and Y data sizes are inconsistent'); end
X = [X ones(size(X,1),1)];
c = [1 0]'; %contrast
df = n - rank(X);
pXX = pinv(X)*pinv(X)'; % = pinv(X'*X), which is reusable, because
pX  = pXX * X';         % pinv(P*X) = pinv(X'*P'*P*X)*X'*P' = pXX * (P*X)'
% Original design (identity permutation)
t = glm_quick_t(Y, X, pXX, pX, df, c);
if any(~isfinite(t(:)))
    warning('glm_permSub zeroed NaN t-scores'); %CR 3Oct2016
    t(~isfinite(t))= 0; %CR patch
end
uncZ = spm_t2z(t,df); %report Z scores so DF not relevant
if nPerms < 2
    threshMin = -Inf;
    threshMax = Inf;
    return
end
if numel(hdrTFCE) == 3
   img = zeros(hdrTFCE);
   img(good_idx) = t;
   img = tfceMex(img, 100);
   t = img(good_idx);
   uncZ = spm_t2z(t,df); %report Z scores so DF not relevant
end
% Things to track over permutations
peak = nan(nPerms,1);
nadir = nan(nPerms,1);
peak(1) = max(t);
nadir(1) = min(t);
perm5pct = round(nPerms * 0.05);
startTime = tic;
for p = 2:nPerms
    if p == perm5pct
        fprintf('Expected permutation time is %f seconds\n',toc(startTime)*20)
    end
    Xp  = X(randperm(n), :);
    pXX = pinv(Xp)*pinv(Xp)'; %??? supposedly not require- reusable?
    pXp = pXX * Xp'; % = pinv(Xp)
    tp  = glm_quick_t(Y, Xp, pXX, pXp, df, c);
    if numel(hdrTFCE) == 3
       img = zeros(hdrTFCE);
       img(good_idx) = tp;
       img = tfceMex(img, 100);
       tp = img(good_idx);
    end
    peak(p) = max(tp(:));
    nadir(p) = min(tp(:));
end
threshMin = permThreshLowSub (nadir, kPcrit);
threshMax = permThreshHighSub (peak, kPcrit);
threshMin = spm_t2z(threshMin,df); %report Z scores so DF not relevant
threshMax = spm_t2z(threshMax,df); %report Z scores so DF not relevant
%end glm_permSub()

function [uncZ, threshMin, threshMax] = glm_permSub_nuisance(Y, X, nuisance, nPerms, kPcrit, good_idx, hdrTFCE)
%returns uncorrected z-score for all voxels in Y given single column
%predictor X and a matrix of nuisance regressors in Z
% Y: volume data
% X: single column predictor
%if either X or Y is binomial, results are pooled variance t-test, else correlation coefficient
% nPerms: Number of permutations to compute
% kPcrit: Critical threshold
%Example
% glm_t([1 1 0 0 0 1; 0 0 1 1 1 0]',[1 2 3 4 5 6]') %pooled variance t-test
%
%inspired by Ged Ridgway's glm_perm_flz
if ~exist('nPerms','var')
    nPerms = 0;
end
if ~exist('kPcrit','var')
    kPcrit = 0.05;
end
% Basics and reusable components
[n f] = size(X); %rows=observations, columns=factors
[nY v] = size(Y); %#ok<NASGU> v is number of tests/voxels
if (f ~= 1), error('glm_permSub is for one column of X at a time (transpose?)'); end;
if nY ~= n, error('glm_permSub X and Y data sizes are inconsistent'); end
X = [X nuisance ones(size(X,1),1)];
c = [1 zeros(1, size(nuisance, 2)) 0]'; %contrast (incl. nuisance regressors)
df = n - rank(X);
pXX = pinv(X)*pinv(X)'; % = pinv(X'*X), which is reusable, because
pX  = pXX * X';         % pinv(P*X) = pinv(X'*P'*P*X)*X'*P' = pXX * (P*X)'
% Original design (identity permutation)
t = glm_quick_t(Y, X, pXX, pX, df, c);
if any(~isfinite(t(:)))
    warning('glm_permSub zeroed NaN t-scores'); %CR 3Oct2016
    t(~isfinite(t))= 0; %CR patch
end
uncZ = spm_t2z(t,df); %report Z scores so DF not relevant
if nPerms < 2
    threshMin = -Inf;
    threshMax = Inf;
    return
end
if numel(hdrTFCE) == 3
   img = zeros(hdrTFCE);
   img(good_idx) = t;
   img = tfceMex(img, 100);
   t = img(good_idx);
   uncZ = spm_t2z(t,df); %report Z scores so DF not relevant
end
% Things to track over permutations
peak = nan(nPerms,1);
nadir = nan(nPerms,1);
peak(1) = max(t);
nadir(1) = min(t);
perm5pct = round(nPerms * 0.05);
startTime = tic;
for p = 2:nPerms
    if p == perm5pct
        fprintf('Expected permutation time is %f seconds\n',toc(startTime)*20)
    end
    Xp  = X(randperm(n), :);
    pXX = pinv(Xp)*pinv(Xp)'; %??? supposedly not require- reusable?
    pXp = pXX * Xp'; % = pinv(Xp)
    tp  = glm_quick_t(Y, Xp, pXX, pXp, df, c);
    if numel(hdrTFCE) == 3
       img = zeros(hdrTFCE);
       img(good_idx) = tp;
       img = tfceMex(img, 100);
       tp = img(good_idx);
    end
    peak(p) = max(tp(:));
    nadir(p) = min(tp(:));
end
threshMin = permThreshLowSub (nadir, kPcrit);
threshMax = permThreshHighSub (peak, kPcrit);
threshMin = spm_t2z(threshMin,df); %report Z scores so DF not relevant
threshMax = spm_t2z(threshMax,df); %report Z scores so DF not relevant
%end glm_permSub()




function savePowerSub(les, beh ,good_idx, hdr, roi_names, beh_names)
%save components required for a power analysis
m.les = les;
m.beh = beh;
m.good_idx = good_idx;
m.hdr = hdr;
m.roi_names = roi_names;
m.beh_names = beh_names;
save('power.mat', '-struct', 'm');
%end savePowerSub()

function [uncZ, threshMin, threshMax] = glm_perm_flSub(Y, X, c, nPerms, kPcrit, good_idx, hdrTFCE)
% [UncZ threshLo threshHi] = glm_perm_fl(Y, X, c, nPerms, pClus)
%glm_perm_flz2: A simple Freedman-Lane permutation test for a t-contrast
% Usage:  [uncZ, threshLoZ, threshHiZ] = glm_perm_fl(Y, X, c, nPerms, pCrit)
%
%Example data from O'Gorman
%  O'Gorman TW (2005) The performance of randomization tests that use permutations of
%   independent variables. Communications in Statistics Simulation and Computation, 34(4): 895-908.
% Sales = [89.8 121.3 115.2 100.3 123 124.8 120 155 200.4 123.6 109.9 82.1 102.4 124.8 134.6 108.5 114 155.8 115.9 128.5 123.5 124.3 128.6 104.3 93.4 121.3 111.2 108.1 189.5 265.7 120.7 90 119 172.4 93.8 121.6 108.4 157 107.3 123.9 103.6 92.7 99.8 106.4 65.5 122.6 124.3 96.7 114.5 106.4 132.2];
% Price = [42.7 41.8 38.5 38.8 39.7 31.1 45.5 41.3 32.6 43.8 35.8 36.7 33.6 41.4 32.2 38.5 38.9 30.1 39.3 38.8 34.2 41 39.2 40.1 37.5 36.8 34.7 34.7 44 34.1 41.7 41.7 41.7 29.4 38.9 38.1 39.8 29 44.7 40.2 34.3 38.5 41.6 42 36.6 39.5 30.2 40.3 41.6 40.2 34.4];
% Income = [2948 4644 3665 2878 4493 3855 4917 4524 5079 3738 3354 4623 3290 4507 3772 3751 3853 3112 3090 3302 4309 4340 4180 3859 2626 3781 3500 3789 4563 3737 4701 3077 4712 3252 3086 4020 3387 3719 3971 3959 2990 3123 3119 3606 3227 3468 3712 4053 3061 3812 3815];
% Sex = [51.7 45.7 50.8 51.5 50.8 50.7 51.5 51.3 53.5 51.8 51.4 48 50.1 51.5 51.3 51.4 51 50.9 51.4 51.3 51.1 52.2 51 51 51.6 51.8 50 51.2 49.3 51.1 51.6 50.7 52.2 51 49.5 51.5 51.3 51 52 50.9 50.9 50.3 51.6 51 50.6 51.1 50.6 50.3 51.6 50.9 50];
% X = [Sex; Price; Income; ones(1,length(Income))]; %regressors including constant intercept
% c = [1 0 0 0]'; %model how first regressor (Sex) predicts Y, other columns are nuisance regressors
%[UncZ threshLo threshHi] = glm_perm_fl(Sales', X', c, 16000);
%
% Copyright 2010 Ged Ridgway
% http://www.mathworks.com/matlabcentral/fileexchange/authors/27434
if ~exist('nPerms','var')
    nPerms = 5000;
end
if ~exist('kPcrit','var')
    kPcrit = 0.05;
end
% Basics and reusable components
[n p] = size(X);
[nY ~] = size(Y);
if nY ~= n, error('Size of data and design are inconsistent'); end
c0 = eye(p) - c * pinv(c);
X0 = X * c0;
R0 = eye(n) - X0 * pinv(X0);
Y  = R0 * Y; % pre-regress, but keep whole design below, for Freedman-Lane
df = n - rank(X);
pXX = pinv(X)*pinv(X)'; % = pinv(X'*X), which is reusable, because
pX  = pXX * X';         % pinv(P*X) = pinv(X'*P'*P*X)*X'*P' = pXX * (P*X)'
% Original design (identity permutation)
t = glm_quick_t(Y, X, pXX, pX, df, c);
uncZ = spm_t2z(t,df); %report Z scores so DF not relevant
if nPerms < 10
    threshMin = -Inf;
    threshMax = Inf;
    return
end
if numel(hdrTFCE) == 3
   img = zeros(hdrTFCE);
   img(good_idx) = t;
   img = tfceMex(img, 100);
   t = img(good_idx);
   uncZ = spm_t2z(t,df); %report Z scores so DF not relevant
end
% Things to track over permutations
peak = nan(nPerms,1);
nadir = nan(nPerms,1);
peak(1) = max(t);
nadir(1) = min(t);
perm5pct = round(nPerms * 0.05);
startTime = tic;
for p = 2:nPerms
    if p == perm5pct
        fprintf('Expected permutation time is %f seconds\n',toc(startTime)*20)
    end
    Xp  = X(randperm(n), :);
    pXp = pXX * Xp'; % = pinv(Xp)
    tp  = glm_quick_t(Y, Xp, pXX, pXp, df, c);
    if numel(hdrTFCE) == 3
       img = zeros(hdrTFCE);
       img(good_idx) = tp;
       img = tfceMex(img, 100);
       tp = img(good_idx);
    end
    peak(p) = max(tp);
    nadir(p) = min(tp);
end
threshMin = permThreshLowSub (nadir, kPcrit);
threshMax = permThreshHighSub (peak, kPcrit);
threshMin = spm_t2z(threshMin,df); %report Z scores so DF not relevant
threshMax = spm_t2z(threshMax,df); %report Z scores so DF not relevant
%end glm_perm_flSub()

function s = glm_quick_t(y, X, pXX, pX, df, c)
b = pX * y;                     % parameters
y = y - X * b;                  % residuals
s = sum(y .* y);                % sum squared error
s = sqrt(s .* (c'*pXX*c) / df); % standard error
s = c'*b ./ s;                  % t-statistic
%end glm_quick_t()

function thresh = permThreshLowSub (permScores, kPcrit)
permScores = sort(permScores(:));
thresh =permScores(round(numel(permScores) * kPcrit));
%report next most significant score in case of ties
permScores = permScores(permScores < thresh);
if ~isempty(permScores)
    thresh = max(permScores(:));
end
%permThreshLowSub()

function thresh = permThreshHighSub (permScores, kPcrit)
permScores = sort(permScores(:),'descend');
thresh =permScores(round(numel(permScores) * kPcrit));
%report next most significant score in case of ties
permScores = permScores(permScores > thresh);
if ~isempty(permScores)
    thresh = min(permScores(:));
end
%permThreshHighSub()
