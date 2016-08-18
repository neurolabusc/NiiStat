function pth = nii_update_mat (pth)
%returns path to images. 
% pth : folder with Excel spreadsheet
%If nii or mat files exit in pth, then assume standard analysis
% otherwise, assume images come from a Github repository.
% the software will look for a repository, and if one does not exist user
% will be prompted to link to one.
%Examples
%  nii_update_mat; %Excel file in current directory
%  nii_update_mat('~\dir'); %path to excel file
%  nii_update_mat(1); %create new NiiMatM repo

repos = { 'M' 'LIME' 'DEMO'};
defaultRepo = 3; %DEMO is 3rd item of repos

if ~exist('pth','var'), pth = pwd; end;
if isnumeric(pth) %e.g. nii_update_mat(1);
    defaultRepo = pth;
    if (pth < 1) || (pth > numel(repos))
        error('Repository number must be between 1..%d', numel(repos));
    end
    pth = uigetdir(pwd,'Select location for repository');
    if isnumeric(pth), return; end; %user pressed cancel
else
    checkForUpdate(pth);
    %if images exist in the path assume we are not using a repository
    if ~isempty(dir(fullfile(pth,'*.nii')))
        fprintf('Reading images from %s\n', pth);
        return; 
    end
    if ~isempty(dir(fullfile(pth,'*.mat')))
        if hasNiiStatMatFilesSub(pth)
            checkForUpdate(pth);
            fprintf('Reading images from %s\n', pth);
            return; 
        end
    end
    %otherwise, see if we can find a github repository
    nameFolds=subFolderSub(pth);
    if ~isempty(nameFolds)
        for i = 1:numel(nameFolds)
            pthx = fullfile(pth, nameFolds{i});
            if ~isempty(dir(fullfile(pthx,'*.mat')))
                if hasNiiStatMatFilesSub(pthx)
                    pth = pthx;
                    fprintf('Reading images from %s\n', pth);
                    checkForUpdate(pth)
                    return; 
                end
            end
        end
    end
end; %if pth is numeric else check folder
repoStr = sprintf('%s;', repos{:});
prompt = {sprintf('Repo (%s):',repoStr),...
        'Username:',...
        'Password:'};
dlg_title = ['Download images from GitHub?'];
num_lines = [1 50];
def = {repos{defaultRepo},'',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answer), return; end;
pth = fullfile(pth, ['NiiMat', answer{1}]);
if isdir(pth)
    rmdir(pth,'s');
end
usr_passwd = '';
if ~isempty(answer{2})
    usr_passwd = answer{2};
    if ~isempty(answer{3})
        usr_passwd = [usr_passwd, ':', answer{3}];
    end
    usr_passwd = [usr_passwd, '@'];
end
gitConfigCmd = ['git config --global http.sslVerify false'];
gitConfigCmdUndo = ['git config --global http.sslVerify true'];
fprintf('\nChanging your git configuration to allow bypassing ssl certificate verification...\n%s\n\n',gitConfigCmd);
fprintf('\nTo undo this in the future run the following command in terminal:\n%s\n\n',gitConfigCmdUndo);
system(gitConfigCmd);
repourl = ['https://',usr_passwd,'gitlab.com/neurolabusc/NiiMat', answer{1}, '.git'];
cmd = sprintf('%s clone %s %s %s', findGitSub(), repourl, pth);
fprintf('%s\n', cmd);
[s, r] = system(cmd,'-echo');
if strfind(r,'fatal')
	warning('Unable to check for updates. Network issue? Matlab blocking commands?');
    clipboard('copy', cmd);
    uiwait(msgbox('A command has been copied. Use CTRL + V to paste it in Terminal, then press enter. Click OK to open Terminal'));
    system('open -a Terminal');
    %error('Paste ');

    return;
end
if ismac 
    pthlfs = fullfile(fileparts(pth), 'lfs');
    if isdir(pthlfs)
        rmdir(pthlfs,'s')
    end
end
if ~exist(pth,'file'), error('Unable to find %s\n',pth); pth = []; end;
%end nii_update_mat()

function checkForUpdate(repoPath)
prevPath = pwd;
cd(repoPath);
if exist('.git','dir') %only check for updates if program was installed with "git clone"
    [s, r] = system([findGitSub(), ' fetch origin'],'-echo');
    if strfind(r,'fatal')
        warning('Unable to check for updates. Network issue?');
        cd(prevPath); %CR 8/2016
        return;
    end
    [~, r] = system([findGitSub(), ' status'],'-echo');
    if strfind(r,'behind')
        if askToUpdateSub
            [~, r] = system([findGitSub(), ' pull'],'-echo');
            %showRestartMsg %no need: we reload images
        end
    end
else %do nothing for now
    %warning(sprintf('To enable updates run "!git clone git@github.com:neurolabusc/%s.git"',mfilename));
    warning('Reading images from a folder that is not set up for updating (perhaps delete folder and restart):',repoPath);
end
cd(prevPath);
%end checkForUpdate()

function nameFolds=subFolderSub(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
%end subFolderSub()

function isNiiStat = hasNiiStatMatFilesSub(pth)
isNiiStat = false;
if ~exist('pth','var') || isempty(pth), pth = pwd; end;
%find all mat files
f = dir(fullfile(pth,'*.mat'));
%if any file has valid NiiStat modality than exit
mod = nii_modality_list;
for i = 1: numel(f)
    fnm = fullfile(pth, f(i).name);
    mat = load(fnm);
    for m = 1 : size(mod,1)
        if isfield(mat, mod(m,:))
            isNiiStat = true;
            return;
        end
    end
end
%end hasNiiStatMatFilesSub()

function a = askToUpdateSub
% Construct a questdlg
choice = questdlg(sprintf('An update for %s is available. Would you like to update?',fileparts(pwd)), ...
	'Auto update', ...
	'Yes','No','Yes');
% Handle response
switch choice
    case 'Yes'
        a = true;
    case 'No'
        a = false;
end
%end askToUpdateSub()

function fnm = findGitSub()
fnm = 'git';
if ismac 
    pth = fileparts(mfilename('fullpath'));
    fnmlfs = fullfile(pth, 'gitlfs', 'git-lfs');
    if exist(fnmlfs, 'file')
        fnm = fnmlfs ;
    end
end
%end findGitSub()