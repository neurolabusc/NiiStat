function pth = nii_update_mat (pth)
%returns path to images. 
% pth : folder with Excel spreadsheet
%If nii or mat files exit in pth, then assume standard analysis
% otherwise, assume images come from a Github repository.
% the software will look for a repository, and if one does not exist user
% will be prompted to link to one.

if ~exist('pth','var'), pth = pwd; end;
%if images exist in the path assume we are not using a repository
if ~isempty(dir(fullfile(pth,'*.mat'))) || ~isempty(dir(fullfile(pth,'*.mat')))
    fprintf('Reading images from %s\n', pth);
    return; 
end
%otherwise, see if we can find a github repository
nameFolds=subFolderSub(pth);
if ~isempty(nameFolds)
    for i = 1:numel(nameFolds)
        pthx = fullfile(pth, nameFolds{i});
        if ~isempty(dir(fullfile(pthx,'*.mat'))) || ~isempty(dir(fullfile(pthx,'*.mat')))
            pth = pthx;
            
            fprintf('Reading images from %s\n', pth);
            checkForUpdate(pth)
            return; 
        end
    end
end
prompt = {'Repo:',...
        'Username:',...
        'Password:'};
dlg_title = ['Unable to find *>mat/*.nii files: get from GitHub?'];
num_lines = 1;
def = {'LIME','crorden6',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
pth = fullfile(pth, ['NiiMat', answer{1}]);
if isdir(pth)
    rmdir(pth,'s');
end
repourl = ['https://',answer{2},':',answer{3},'@gitlab.com/neurolabusc/NiiMat', answer{1}, '.git'];
cmd = sprintf('git clone %s %s', repourl);
[s, r] = system(cmd,'-echo');
if strfind(r,'fatal')
	warning('Unable to check for updates. Network issue?');
    return;
end
if ~exist(pth,'file'), error('Unable to find %s\n',pth); pth = []; end;
%end nii_update_mat()

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
            [~, r] = system('git pull','-echo');
            showRestartMsg
        end
    end
else %do nothing for now
    warning(sprintf('To enable updates run "!git clone git@github.com:neurolabusc/%s.git"',mfilename));
end
cd(prevPath);
%end checkForUpdate()

function nameFolds=subFolderSub(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
%end subFolderSub()
