function roifileresults = GetROIFilesGivenDirectory(roidirname)

            onefile = which('jhu.txt');
            onefilelength = length(onefile);
            roifolder = fullfile(onefile(1:onefilelength-11),roidirname);
            
            %find roi files in that folder
            files = dir(roifolder);
            L = length(files);
            index = false(1, L);
            for k = 1:L
                M = length(files(k).name);
                if M > 4 && strcmp(files(k).name(M-2:M), 'txt')
                    index(k) = true;
                end
            end
            roifileresults = files(index);

end
