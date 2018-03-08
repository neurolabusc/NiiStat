function roifiledir = GetROIPathGivenDirectoryName(roidirname)

            onefile = which('jhu.txt');
            onefilelength = length(onefile);
            roifiledir = fullfile(onefile(1:onefilelength-11),roidirname);

end