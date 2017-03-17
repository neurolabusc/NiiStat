function nii_stat_svm(les,beh, beh_names, statname, les_names, subj_data, roifname)

    
if numel(les_names) ~= size(les,2) %for correlation analyses
 les_matrix = [];
 n = 0;
 for i = 1:length(les_names)
    for j = (i+1):length (les_names)
       area_name = [les_names{i} '*' les_names{j}];
       %fprintf('%s\n',area_name);
       n = n+1;
       les_matrix{n} = area_name; %#ok<AGROW>
    end
 end
 les_names = les_matrix;
end
if numel(les_names) ~= size(les,2)
    fprintf('%s error: number of feature names does not match number of features',mfilename);
    return;
end

if ~exist('statname','var')
    statname = 'anonymous';
end
if ~exist('subj_data','var')
    subj_data = [];
end;
chDirSub([statname '_svm']);
diary ([deblank(statname) 'svm.txt']);
for j = 1:size(beh_names,2) %for each beahvioral variable
    beh_name1 = beh_names{j};
    beh1 = beh(:,j);
    [fnm, nOK] = tabFileSub(les,beh1, beh_name1,  les_names, subj_data);
    if nOK < 1
        fprintf('Skipping SVM/SVR: no valid data\n');
    else
        [~, loadingMap] = nii_stat_svm_core(fnm); %do not specify thresholds: svm_core will select
        if exist('roifname','var') && ~isempty(roifname)
           nii_array2roi (loadingMap, roifname, [statname '_svm.nii']) 
        end
        if ~nii_isBinary(beh1)
            [~, loadingMap] = nii_stat_svr_core(fnm); %compute regression
            if exist('roifname','var') && ~isempty(roifname)
                nii_array2roi (loadingMap, roifname, [statname '_svr.nii']) 
            end
        end
    end
    %nii_stat_svm_core(fnm, min(beh1(:)), 0.5+min(beh1(:)) );
    %nii_stat_svm_core(fnm, min(beh1(:)), max(beh1(:)) );
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
        if ~isempty('subj_data')
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