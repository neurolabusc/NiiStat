function [acc, z_map, acc_per_class, p] = nii_stat_svm_core (data, class_labels, maxNsplits, normRowCol, verbose, islinear)
% classify the data using linear support vector machine
% INPUTS:
% data (#observations x #features) -- data matrix (ROI or voxelwise)
% class_labels -- zeros and ones for observations in class 1 and 2
% maxNsplits -- number of splits in the cross-validation procedure
%               (recommended: 100 for voxelwise analysis, 500 for ROI analysis)
% normRowCol -- normalize none [0], rows [1, default], or columns [2]
% verbose -- if true, libsvm generates a lot of messages
% islinear -- should we use linear kernel for SVM (true by default, set to
%             false at your own risk)
% OUTPUTS:
% acc -- classification accuracy (proportion of correct out-of-sample assignments)
% z_map -- map of Z-scored feature weights (features being ROIs or voxels)
% acc_per_class -- accuracy for, separately, class 1 and class 2
% p -- p-value of classification accuracy, assuming binomial distribution 

% should we optimize C parameter for SVM? it's worth it, but takes time
optimize_C = true;

if ~exist('maxNsplits','var')
    maxNsplits = 100; 
end
if ~exist('normRowCol','var')
    normRowCol = 1; %none [0], rows [1], or columns [2]
end
if ~exist('verbose','var') %vebosity not specified
    verbose=false;
end
svmdir = [fileparts(which(mfilename))  filesep 'libsvm' filesep];
if ~exist(svmdir,'file') 
    error('Unable to find utility scripts folder %s',svmdir);
end
if ~exist('islinear','var')
    islinear = true;
end
if ~islinear
    fprintf('Warning: please do not use loading map when conducting non-linear svm!');
end
addpath(svmdir); %make sure we can find the utility scripts


%for aalcat we may want to remove one hemisphere
%data(:,1:2:end)=[]; % Remove odd COLUMNS: left in AALCAT: analyze right
%data(:,2:2:end)=[]; % Remove even COLUMNS: right in AALCAT: analyze left
[data, good_idx]  = requireVarSub(data);
%data = requireNonZeroSub(data); %CBF only

%normalize each predictor for range 0..1 to make magnitudes similar
%
%data  = normColSub(data);
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
    fprintf('Normalizing so each column has range 0..1\n');
    data = normColSub(data);
end
%data  = normRowSub(data);
%binarize class_label: either 0 or 1
if (min(class_labels) == max(class_labels))
    error('No variability in class labels (final column of file)');
end
if ~exist('thresholdHi','var') %if Excel file not specified, have user select one
    fprintf('Class values from %g..%g (median %g)\n',min(class_labels(:)), max(class_labels(:)), median(class_labels(:)));
    mdn = median(class_labels(:));
    if mdn > min(class_labels(:))
        thresholdHi = mdn;
        %thresholdLo = thresholdHi-eps;
        thresholdLo = max(class_labels(class_labels < mdn));
    else
        thresholdLo = mdn;
        thresholdHi = min(class_labels(class_labels > mdn));
    end
    %answer = inputdlg({'Class1 no more than','Class2 no less than'}, sprintf('Threshold (input range %g..%g)',min(class_labels),max(class_labels)), 1,{num2str(min(class_labels)),num2str(max(class_labels))});
    %thresholdLo = str2double (cell2mat(answer(1)));
    %thresholdHi = str2double (cell2mat(answer(2)));
end
%fprintf('Class labels range from %g to %g, values of %g or less will be group0, values of %g or more will be group1\n', min(class_labels), max(class_labels), thresholdLo, thresholdHi);
%fprintf('Processing the command line: \n');
%fprintf(' %s (''%s'', %d, %g, %g, %g, %g);\n', mfilename, xlsname, normRowCol, thresholdLo, thresholdHi, verbose, islinear);

[class_labels , data] = binarySub(class_labels, data, thresholdLo, thresholdHi);

if islinear 
    cmd = '-t 0';
else
    cmd = '-t 2';
end
    
class1_idx = find (class_labels == 1)';
class0_idx = find (class_labels == 0)';
n0 = length (class0_idx); 
n1 = length (class1_idx); % # examples per class
N = n0 + n1;
min_n = min (n0, n1);
if (n1 < 2) || (n0 < 2)
    fprintf('Each group must have at least 2 observations: please use a different class threshold\n');
    z_map = []; acc = []; acc_per_class = []; p = [];
    return;
end
data = data'; 

if optimize_C
    N_splits = 20; % quick and dirty split-half cross-validation to select C
    C_list = [0.001 0.0025 0.005 0.01 0.025 0.05 0.1 0.25 0.5 1 2.5 5];
    for g = 1:N_splits
        shuffle = randperm (n0);
        idx01 = class0_idx (sort(shuffle (1:floor(min_n/2))));
        idx02 = class0_idx (sort(shuffle (floor(min_n/2)+1:min_n)));
        shuffle = randperm (n1);
        idx11 = class1_idx (sort(shuffle (1:floor(min_n/2))));
        idx12 = class1_idx (sort(shuffle (floor(min_n/2)+1:min_n)));
        
        data1 = data (:, [idx01 idx11]); data2 = data (:, [idx02 idx12]);
        class_labels1 = class_labels ([idx01 idx11]); class_labels2 = class_labels ([idx02 idx12]);
        for C_idx = 1:length(C_list)
            str = sprintf ('''%s -c %g''', cmd, C_list(C_idx));
            [out, subSVM1] = evalc (['svmtrain (class_labels1, data1'', ' str ');']);
            [out, subSVM2] = evalc (['svmtrain (class_labels2, data2'', ' str ');']);
            ww1 = subSVM1.sv_coef' * subSVM1.SVs;
            ww2 = subSVM2.sv_coef' * subSVM2.SVs;
            temp = corrcoef (ww1, ww2);
            repr(C_idx, g) = temp (1, 2);
            [~, ~, temp, ~] = evalc ('svmpredict (class_labels2, data2'', subSVM1)');
            acc1 = temp(1)/100;
            [~, ~, temp, ~] = evalc ('svmpredict (class_labels1, data1'', subSVM2)');
            acc2 = temp(1)/100;
            acc(C_idx, g) = mean ([acc1 acc2]);
        end  
    end
    cost = ((1+acc)/2).^2 + repr.^2;    
    [~, optC_idx] = max (mean (cost, 2));
    C = C_list (optC_idx);
    fprintf ('Optimized value of C: %g\n', C);
else
    C = 0.01;
    fprintf ('No optimization of C; using default of %g\n', C);
end
cmd = sprintf ('%s -c %g', cmd, C);
clear acc;

n_train = min_n - 1; % number of training examples per class
n_test0 = n0 - n_train; % number of test examples in class 1
n_test1 = n1 - n_train; % ... and in class 2
% maximum number of training-test splits
warning ('OFF', 'MATLAB:nchoosek:LargeCoefficient');
theoretical_max_n_splits = nchoosek (n1, n_test1) * nchoosek (n0, n_test0);
n_splits = min (theoretical_max_n_splits, maxNsplits); 
correct = zeros (size (class_labels));
ntrials = zeros (size (class_labels)); % number of times each example was tested
prev_splits = zeros (1, 2*n_train);
curr_split = zeros (1, 2*n_train); % reset at first iteration
for split = 1:n_splits
    % make sure a newly created split is unique
    while ~isempty (find (ismember (prev_splits, curr_split, 'rows'))) %#ok<EFIND>
        shuffle = randperm (n1);
        train_idx1 = class1_idx(sort (shuffle (1:n_train)));
        test_idx1 = class1_idx(sort (shuffle (n_train+1:n1)));
        shuffle = randperm (n0);
        train_idx0 = class0_idx(sort (shuffle (1:n_train)));
        test_idx0 = class0_idx(sort (shuffle (n_train+1:n0)));
        curr_split = [train_idx1 train_idx0];
    end
    prev_splits (split, :) = curr_split;
    test_idx = [test_idx1 test_idx0];
    train_idx = [train_idx1 train_idx0];
    test_data = data (:, test_idx);
    test_labels = class_labels (test_idx)';
    train_data = data (:, train_idx);
    train_labels = class_labels (train_idx)';  
   	if verbose
        model = svmtrain (train_labels', train_data', cmd);
        [assignments, ~, ~] = svmpredict (test_labels', test_data', model);
    else %if verbose else silent
        [~, model] = evalc ('svmtrain (train_labels'', train_data'', cmd)'); %-t 0 = linear
        [~, assignments, ~, ~] = evalc ('svmpredict (test_labels'', test_data'', model)');
    end %if verbose...
    
    map (split, :) = model.sv_coef' * model.SVs; %#ok<AGROW>
    ntrials (test_idx) = ntrials (test_idx) + 1;
    correct_idx = test_idx (find (test_labels == assignments'));
    correct (correct_idx) = correct (correct_idx) + 1;    
end %for each split  
% mean_map = mean (map, 1);
% z_map = mean_map ./ std (mean_map);
% Z-scoring of maps: a version of delete-d jackknife
d = abs (n0 - n1) + 2; 
t_map = mean (map, 1) ./ (std (map, 1, 1) * sqrt(N/d-1));
z_map = zeros (size (t_map));
z_map (:) = nan;
z_map (~isnan (t_map)) = spm_t2z (t_map(~isnan (t_map)), length(class_labels) - 1);
if exist('good_idx','var')  %insert NaN for unused features
    z_mapOK = zeros(size(data, 1), 1);
    z_mapOK(:) = nan;
    z_mapOK(good_idx) = z_map;
    z_map = z_mapOK;
end
temp = correct ./ ntrials;
acc = mean (temp (find (~isnan (temp)))); %#ok<*FNDSB>
acc_per_class(1) = mean (temp (intersect (class1_idx, find (~isnan (temp)))));
acc_per_class(2) = mean (temp (intersect (class0_idx, find (~isnan (temp)))));

%report results
prob = max(n0,n1)/(n0+n1);
p = bipSub(acc*n_splits, n_splits, prob);
fprintf('Observed %d in group0 and %d in group1 (prob = %g%%) with %d predictors\n', n0, n1,max(n0,n1)/(n0+n1), size(data,1));
fprintf('BinomialProbality nHits= %g, nTotal= %g, IncidenceOfCommonClass= %g, p< %g\n', acc*n_splits, n_splits,prob, p);
fprintf('Overall Accuracy %g (%g for group0, %g for group1)\n', acc, acc_per_class(2),acc_per_class(1));
fprintf('Thresh0\t%g\tThresh1\t%g\tn0\t%d\tn1\t%d\tProb\t%g\tAcc\t%g\tAcc0\t%g\tAcc1\t%g\tp<\t%g\n',...
    thresholdLo,thresholdHi,n0,n1,prob,acc,acc_per_class(2),acc_per_class(1),p);
% ----- MAIN FUNCTION ENDS

function prob = bipSub(x,n,p)
%report extreme tail: chance for equal or more extreme values
if x > (p*n)
    h = n - x;
    px = 1-p;
else
   h = x; 
   px = p;
end
prob = binocdfSub(h, n, px);
%end

function cdf = binocdfSub (x, n, p)
% probability of at least x or more correct in n attempts, where probability of correct is p
%  x : number of successes, in range 0..n (for accuracy of 0..100%)
%  n : number of attempts (integer)
%  p : chance probability of success
% http://en.wikipedia.org/wiki/Binomial_distribution
%Normal Approximation when values will overflow factorial version
% http://en.wikipedia.org/wiki/Binomial_distribution#Normal_approximation
% http://code.metager.de/source/xref/gnu/octave/scripts/statistics/distributions/binocdf.m
if  ((n*p) > 9) && ((n * (1-p)) > 9) %use approximation
    k = (x >= 0) & (x < n) & (n == fix (n)) & (p >= 0) & (p <= 1);
    tmp = floor (x(k));
    cdf = betainc (1 - p, n - tmp, tmp + 1);
    return
end
cdf = 0.0;
for i = 0:x
        cdf = cdf + nchoosek(n, i)*(power(p,i) )*(power(1.0-p,(n-i) ) );
end

% function p = bipSub(nHits, nTotal, prob)
% % probability of at least nHits correct classifications in nTotal attempts
% % nHits is in range 0..n (for accuracy of 0..100%)
% % http://en.wikipedia.org/wiki/Binomial_distribution
% %Approximation suitable for large values that would overflow factorial version
% % http://www.dummies.com/how-to/content/how-to-find-the-normal-approximation-to-the-binomi.html
% np = nTotal * prob;
% n_not_p = nTotal * (1-prob);
% if (np > 10) && (n_not_p > 10)
%     var = sqrt(np * (1.0 - prob));
%     z = (nHits-np)/var;
%     p = 1 - spm_Ncdf(z);
%     return
% end
% if nHits > (prob*nTotal)
%     nHitsSmallTail = nTotal - nHits;
% else
%    nHitsSmallTail = nHits; 
% end
% p = 0.0;
% for lHitCnt = 0:nHitsSmallTail
%     p = p + nchoosek(nTotal, lHitCnt)*(power(prob,lHitCnt) )*(power(1.0-prob,(nTotal-lHitCnt) ) );
% end
% %end bipSub()

% function x = normColSub(y)
% %normalize each column for range 0..1
% x = y';
% mn = min(x); %minimum
% rng = max(x) - mn; %range
% rng(rng == 0) = 1; %avoid divide by zero
% rng = 1 ./ rng; %reciprocal: muls faster than divs
% for i = 1 : numel(mn)
%     x(:,i)=(x(:,i)-mn(i)) * rng(i);
% end
% x = x';
% %normColSub

function x = normColSub(x)
%normalize each column for range 0..1
% x = [1 4 3 0; 2 6 2 5; 3 10 2 2] -> x = [0 0 1 0; 0.5 0.333 0 1; 1 1 0 0.4]
if size(x,1) < 2 %must have at least 2 rows
    fprintf('Error: normalizing columns requires multiple rows\n');
    return
end
if min(max(x,[],1)-min(x,[],1)) == 0
    fprintf('Error: unable to normalize columns: some have no variability\n');
    return;
end
x = bsxfun(@minus,x,min(x,[],1)); %translate so minimum = 0
x = bsxfun(@rdivide,x,max(x,[],1)); %scale so range is 1
%end normColSub()

function x = normRowSub(x)
%normalize each column for range 0..1
% x = [1 4 3 0; 2 6 2 5; 3 10 2 2] -> x = [0.25 1 0.75 0; 0 1 0 0.75; 0.125 1 0 0]
if size(x,2) < 2 %must have at least 2 rows
    fprintf('Error: normalizing rows requires multiple columns\n');
    return
end
if min(max(x,[],2)-min(x,[],2)) == 0
    fprintf('Error: unable to normalize rows: some have no variability\n');
    return;
end
x = bsxfun(@minus,x,min(x,[],2)); %translate so minimum = 0
x = bsxfun(@rdivide,x,max(x,[],2)); %scale so range is 1
%end normRowSub()

function num = tabreadSub(tabname)
%read cells from tab based array. 
fid = fopen(tabname);
num = [];
row = 0;
while(1) 
	datline = fgetl(fid); % Get second row (first row of data)
	%if (length(datline)==1), break; end
    if(datline==-1), break; end %end of file
    if datline(1)=='#', continue; end; %skip lines that begin with # (comments)
    tabLocs= strfind(datline,char(9)); %findstr(char(9),datline); % find the tabs
    row = row + 1;
    if (row < 2) , continue; end; %skip first row 
    if (tabLocs < 1), continue; end; %skip first column
    dat=textscan(datline,'%s',(length(tabLocs)+1),'delimiter','\t');
    for col = 2: size(dat{1},1) %excel does not put tabs for empty cells (tabN+1)
    	num(row-1, col-1) = str2double(dat{1}{col}); %#ok<AGROW>
    end
end %while: for whole file
fclose(fid);
%end tabreadSub()

% function [good_dat, good_idx] = requireNonZeroSub (dat)
% %remove columns with zeros
% good_idx=[];
% for col = 1:size(dat,2)
%     if sum(dat(:,col) == 0) > 0
%        %fprintf('rejecting column %d (non-numeric data')\n',col) %
%     elseif min(dat(:,col)) ~= max(dat(:,col))
%         good_idx = [good_idx, col];  %#ok<AGROW>
%     end
% end %for col: each column
% if numel(good_idx) ~= size(dat,2)
%     fprintf('Some predictors have zeros (analyzing %d of %d predictors)\n',numel(good_idx), size(dat,2));
% end
% good_dat = dat(:,good_idx);
% %end requireNonZeroSub()

function [good_dat, good_idx] = requireVarSub (dat)
good_idx=[];
for col = 1:size(dat,2)
    if sum(isnan(dat(:,col))) > 0
       %fprintf('rejecting column %d (non-numeric data')\n',col) %
    elseif min(dat(:,col)) ~= max(dat(:,col))
        good_idx = [good_idx, col];  %#ok<AGROW>
    end
end %for col: each column
if sum(isnan(dat(:))) > 0
    fprintf('Some predictors have non-numeric values (e.g. not-a-number)\n');
end
if numel(good_idx) ~= size(dat,2)
    fprintf('Some predictors have no variability (analyzing %d of %d predictors)\n',numel(good_idx), size(dat,2));
end
good_dat = dat(:,good_idx);
%end requireVarSub()

%function dataBin = binarySub(data, threshold)
%dataBin = zeros(size(data));
%dataBin(data >= threshold) = 1;
%end binarySub

function [classBin,data] = binarySub(classIn, data, thresholdLo, thresholdHi)
%rows where classIn is <= thresholdLo are assigned group 0 in classBin
%rows where classIn is >= thresholdHi are assigned group 1 in classBin
%rows where classIn is between thresholdLo and ThresholdHi are deleted
classBin = zeros(numel(classIn),1);
classBin(classIn > thresholdLo) = nan;
classBin(classIn >= thresholdHi) = 1;
bad = isnan(classBin);
data(bad,:) = [];
classBin(bad)=[];
%end binarySub

