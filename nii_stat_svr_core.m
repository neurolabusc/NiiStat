function [r, z_map, predicted_labels, p] = nii_stat_svr_core (data, labels, beh_name, clipping, normRowCol, verbose, islinear, nuSVR)
% predict the scores using linear support vector regression
% INPUTS:
% data (#observations x #features) -- data matrix (ROI or voxelwise)
% labels -- scores to predict (i.e. values of dependent variable)
% beh_name -- string containing the name of the dependent variable
% clipping -- 1 if we only want to keep positive feature weights;
%            -1 if we only want to keep negative feature weights; 
%             0 if we want to keep positive and negative feature weights (default)
% normRowCol -- normalize none [0], rows [1, default], or columns [2]
% verbose -- if true, libsvm generates a lot of messages
% islinear -- should we use linear kernel for SVM (true by default, set to
%             false at your own risk)
% nuSVR -- set to 1 if, for some reason, you want to use nu-SVR formulation
% OUTPUTS:
% r -- prediction accuracy (correlation between actual and predicted scores)
% z_map -- map of Z-scored feature weights (features being ROIs or voxels)
% predicted_labels -- predicted score values per observation
% p -- p-value of classification accuracy, assuming binomial distribution

optimize_C = true; 
split_half = true; % use split-half resamplinbg to optimize C (recommended; otherwise, 8-fold CV will be used) 

if ~exist('normRowCol','var')
    normRowCol = 2; %[0, default], rows [1], or columns [2]
end
if ~exist('verbose','var')
    verbose = false;
end
if ~exist('minUnique','var')
   minUnique = 0; 
end
svmdir = [fileparts(which(mfilename))  filesep 'libsvm' filesep];
if ~exist(svmdir,'file') 
    error('Unable to find utility scripts folder %s',svmdir);
end
if ~exist('islinear','var')
    islinear = true;
end
if ~islinear
    fprintf('Warning: please do not use map when conducting non-linear svr!');
end
if ~exist('nuSVR','var') 
    nuSVR = false;
end
addpath(svmdir); %make sure we can find the utility scripts

if normRowCol ==  -2% rows [1], or columns [2], minus means exclude final column
    fprintf('Normalizing so each column has range 0..1 (last column excluded)\n');
    data(:,1:end-1) = normColSub(data(:,1:end-1));
elseif normRowCol ==  -1% rows [1], or columns [2], minus means exclude final column
    fprintf('Normalizing so each row has range 0..1 (last column excluded)\n');
    data(:,1:end-1) = normRowSub(data(:,1:end-1));
elseif normRowCol ==  1% rows [1], or columns [2]
    fprintf('Normalizing so each row has range 0..1\n');
    data = normRowSub(data);
elseif normRowCol ==  2% rows [1], or columns [2]
    %fprintf('Normalizing so each column has range 0..1\n');
    data = normColSub(data);
end
if clipping == 0
    fprintf ('Two-tailed analysis; keeping both positive and negative features\n');
elseif clipping == 1
    fprintf ('One-tailed analysis; keeping only positive features\n');
else
    fprintf ('One-tailed analysis; keeping only negative features\n');
end
labels = normColSub(labels); %CR 20April2015
n_subj = length (labels);
predicted_labels = zeros(n_subj,1); %pre-allocate memory
map = zeros(size(data)); %pre-allocate memory

if islinear 
    cmd = '-t 0';
else
    cmd = '-t 2';
end
if nuSVR
    cmd = [cmd ' -s 4'];
else
    cmd = [cmd ' -s 3'];
end

if optimize_C
    C_list = [0.001 0.0025 0.005 0.01 0.025 0.05 0.1 0.25 0.5 1 2.5 5];
    N = length (labels);
    if split_half
        N_splits = 20;
        for g = 1:N_splits
            shuffle = randperm (N);
            idx1 = shuffle (1:floor(N/2));
            idx2 = shuffle (floor(N/2)+1:N);
            data1 = data (idx1, :); data2 = data (idx2, :);
            scores1 = labels (idx1); scores2 = labels (idx2);
            for C_idx = 1:length(C_list)
                str = sprintf ('''%s -c %g''', cmd, C_list(C_idx));
                [out, subSVM1] = evalc (['svmtrain (scores1, data1, ' str ');']);
                [out, subSVM2] = evalc (['svmtrain (scores2, data2, ' str ');']); 
                [ww1, bb1] = subGetModelWeights (subSVM1, clipping);
                [ww2, bb2] = subGetModelWeights (subSVM2, clipping);
                sub_prediction{C_idx}(idx2) = ww1*data2' + bb1;
                sub_prediction{C_idx}(idx1) = ww2*data1' + bb2;
                sub_map{C_idx}(2*g-1, :) = ww1;
                sub_map{C_idx}(2*g, :) = ww2;                
            end
        end
    else
        % 8-fold cross validation
        s = floor (N/8);
        G = N - 8*s;
        for C_idx = 1:length(C_list)
            sub_prediction{C_idx} = zeros (size (labels));
        end        
        for g = 1:8
            if g <= G
                subtest_idx = g*s+g-s:g*s+g;
            else
                subtest_idx = g*s+G-s+1:g*s+G;
            end
            subtrain_idx = setdiff (1:N, subtest_idx);
            subtrain_data = data (subtrain_idx, :);
            subtrain_scores = labels (subtrain_idx);
            subtest_data = data (subtest_idx, :);
            subtest_scores = labels (subtest_idx);
            
            for C_idx = 1:length(C_list)
                str = sprintf ('''%s -c %g''', cmd, C_list(C_idx));
                [out, subSVM] = evalc (['svmtrain (subtrain_scores, subtrain_data, ' str ');']);
                [ww, bb] = subGetModelWeights (subSVM, clipping);
                sub_prediction{C_idx}(subtest_idx) = ww*subtest_data' + bb;
                sub_map{C_idx}(g, :) = ww;
            end
        end
    end
    
    for C_idx = 1:length(C_list)
        temp = corrcoef (sub_prediction{C_idx}, labels);
        sub_acc(C_idx) = temp (1, 2);
        temp = sub_map{C_idx};
        sub_repr(C_idx) = mean (mean (corrcoef (temp')));
    end
    cost = ((1+sub_acc)/2).^2 + sub_repr.^2;
    [~, opt_idx] = max (cost);
    C = C_list (opt_idx);
    fprintf ('Optimized value of C: %g\n', C);
else
    C = 0.01;
    fprintf ('No optimization of C; using default of %g\n', C);
end
cmd = sprintf ('%s -c %g', cmd, C); 

% run SVR on full data set
[~, SVM] = evalc(sprintf('svmtrain (labels, data, ''%s'')',cmd)');
[fulldata_map, ~] = subGetModelWeights (SVM, clipping);

for subj = 1:n_subj
    train_idx = setdiff (1:n_subj, subj);
    train_data = data (train_idx, :);
    train_labels = labels (train_idx);
    if verbose
        SVM = svmtrain (train_labels, train_data, cmd);       
        % added by Grigori Yourganov: scale correction
        % step 1: predict training labels
        pred_train_labels = svmpredict (train_labels, train_data, SVM);                  
    else %if verbose else silent
        [~, SVM] = evalc(sprintf('svmtrain (train_labels, train_data, ''%s'')',cmd)'); 
        [out, pred_train_labels] = evalc ('svmpredict (train_labels, train_data, SVM);');        
    end %if verbose else silent
    % step 2 of scale correction: estimate scale&offset
    y = train_labels; %regression line: y = a*x + b
    x = pred_train_labels;
    m = length (train_labels);
    c = (m+1)*sum(x.^2) - sum(x)*sum(x);
    a = ((m+1)*sum(x.*y) - sum(x)*sum(y)) / c;
    b = (sum(x.^2)*sum(y) - sum(x)*sum(x.*y)) / c;
    % predict the test labels
    [ww, bb] = subGetModelWeights (SVM, clipping);
    predicted_labels(subj) = ww*data(subj, :)' + bb;
    % step 3 of scale correction: rescale using estimated scale&offset
    predicted_labels(subj) = a*predicted_labels(subj) + b;
    map (subj, :) = ww; % used to be SVM.sv_coef' * SVM.SVs;
end
predicted_labels = normColSub(predicted_labels); % GY

[r, p] = corrcoef (predicted_labels', labels);
r = r(1,2); %r = correlation coefficient
p = p(1,2) /2; %p = probability, divide by two to make one-tailed 
if (r < 0.0) %one tailed: we predict predicted scores should positively correlate with observed scores
    p = 1.0-p; %http://www.mathworks.com/matlabcentral/newsreader/view_thread/142056
end
fprintf('r=%g r^2=%g, p=%g numObservations=%d numPredictors=%d DV=%s\n',r,r^2, p,size (data, 1),size (data, 2), beh_name);
% standardize the maps using jaccknife estimates; GY, Jan 2018; updated May 2018
t_map = mean (map, 1) ./ (std (map, 1, 1) * sqrt(n_subj-1));
z_map = zeros (size (t_map));
z_map (:) = nan;
z_map (~isnan (t_map)) = spm_t2z (t_map(~isnan (t_map)), length(labels) - 1);
% compute reproducibility of maps; unused for now, might be useful in future 
% "pseudomap" is analogous to jackknife "pseudovalue", see e.g. Efron & Tibshirani
pseudomap = n_subj*repmat (fulldata_map, [n_subj 1]) - (n_subj-1)*map;
cc = corrcoef (pseudomap');
reproducibility = mean (cc (find (triu (ones (size (cc)), 1))));
%mean_map = mean (map, 1);
%z_map = mean_map ./ std (mean_map);
if exist('good_idx','var') %insert NaN for unused features
    z_mapOK = zeros(n_dim,1);
    z_mapOK(:) = nan;
    z_mapOK(good_idx) = z_map;
    z_map = z_mapOK;
end
%plot results
figure;
plot (labels, predicted_labels, 'o');
axis ([min(labels(:)) max(labels(:)) min(labels(:)) max(labels(:))]);
%set (gca, 'XTick', [0 1 2 3 4]);
xlabel ('Actual score');
ylabel ('Predicted score');
plot_title = beh_name;
plot_title (plot_title == '_') = '-';
title (sprintf ('%s', plot_title));
%end nii_stat_svr_core()


function x = normColSub(x)
%normalize each column for range 0..1
% x = [1 4 3 0; 2 6 2 5; 3 10 2 2] -> x = [0 0 1 0; 0.5 0.333 0 1; 1 1 0 0.4]
if size(x,1) < 2 %must have at least 2 rows
    fprintf('Error: normalizing columns requires multiple rows\n');
    return
end
x = bsxfun(@minus,x,min(x,[],1)); %translate so minimum = 0
x = bsxfun(@rdivide,x,max(x,[],1)); %scale so range is 1
%end normColSub()

function x = normRowSub(x)
%normalize each column for range 0..1
% x = [1 4 3 0; 2 6 2 5; 3 10 2 2] -> x = [0.25 1 0.75 0; 0 1 0 0.75; 0.125 1 0 0]
if size(x,2) < 2 %must have at least 2 rows
    fprintf('Error: normalizing rows requires multiple columns');
    return
end
if min(max(x,[],2)-min(x,[],2)) == 0
    fprintf('Error: unable to normalize rows: some have no variability');
    return;
end
x = bsxfun(@minus,x,min(x,[],2)); %translate so minimum = 0
x = bsxfun(@rdivide,x,max(x,[],2)); %scale so range is 1
%end normRowSub()

% function x = columnUniqueThreshSub(x, minUnique)
% %removes columns that have fewer than minUnique different values
% n_col = size(x,2);
% for c = n_col:-1:1
%     if numel(unique(x(:,c))) < minUnique
%         x(:,c) = [];
%     end
% end
% if n_col ~= size(x,2)
%     fprintf('%d columns of %d had at least %d unique values\n',size(x,2), n_col, minUnique);
% end
%end columnUniqueThreshSub()

function [good_dat, good_idx] = columnUniqueThreshSub(dat, minUnique)
n_col = size(dat,2);
good_idx=[];
for col = 1:n_col    
    %how many items have the most frequent value
    s = sum(dat(:,col) == mode(dat(:,col)));
    if (size(dat,1) - s) >= minUnique
        good_idx = [good_idx, col];  %#ok<AGROW>
    end
end
if numel(good_idx) ~= size(dat,2)
     fprintf('%d columns of %d had at least %d items in the smaller class\n',numel(good_idx), n_col, minUnique);
end
good_dat = dat(:,good_idx);
%end columnUniqueThreshSub()

function [ww, bb] = subGetModelWeights (model, clipping)
ww = model.sv_coef' * model.SVs; % model weights
bb = -model.rho; % model offset
if (clipping == 1)
    ww (find (ww < 0)) = 0;
elseif (clipping == -1)
    ww (find (ww > 0)) = 0;
end
% end subGetModelWeights

% function x = columnUniqueThreshSubOld(x, minUnique)
% %removes columns that have fewer than minUnique different values
% % x = [ 1 2 3 3 5 6; 2 2 2 2 3 3]'; ->minUnique(3) -> [1 2 3 3 5 6]'
% n_col = size(x,2);
% for c = size(x,2):-1:1    
%     %how many items have the most frequent value
%     s = sum(x(:,c) == mode(x(:,c)));
%     if (size(x,1) - s) < minUnique
%         x(:,c) = [];
%     end
% end
% if n_col ~= size(x,2)
%      fprintf('%d columns of %d had at least %d items in the smaller class\n',size(x,2), n_col, minUnique);
% end
% %end columnUniqueThreshSub()