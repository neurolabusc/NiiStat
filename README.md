##### Introduction

NiiStat is a set of Matlab scripts for analyzing neuroimaging data from clinical populations

This GitHub archive contains the latest release. For usage see the [NITRC website](http://www.nitrc.org/projects/niistat/)

##### Versions

9-October-2016 :
 - Saves node/edge map for connectome analyses (DTI/Rest, GLM/SVM)
 - Supports most NPM val files directly, unifies logical mask code

30-June-2016 :
 - Main program name changed from "nii_stat" to "NiiStat" to match NITRC
 - Behavioral data should be in excel spreadsheet tab named "NiiStat" (old format expected "Data (2)", which still works).
 - Fix bug specific to 6June2016 version where region of interest maps were drawn scrambled.
 - Script automatically checks GitHub for updates.

3-Mar-2015 :
 - Added modalities ALF, fMRI

29-Aug-2014 :
 - Fixed bug where nii_stat would crash when 3 conditions were met: (1) conducting correlation analysis while (2) asked to remove white matter regions from (3) ROI that does not include any white matter regions

3-Sept-2014:
 - New modality: dtifc is the fiber count for diffusion tractography (while dti is the density)

8-Sept-2014:
 - Previous versions did not allow .XLS files where Excel interpreted the filename as a number (e.g. "748" failed, "LM748" was fine)
 - New feature: creates mean intensity map for ROI analyses (similar to sum maps generated for voxelwise analysis).
 - nii_nii2mat will binarize images to 0 or 1 when the modality "lesion" is selected (helps when lesion criteria vary, e.g. for some files 0/1 and others 0/255)
 - nii_nii2mat now adds support for .nii.gz and .voi images (in addition to .hdr and .nii)
 - Excel worksheets no longer crash nii_stat if they have behavioral variables with a "/" or "\" in their names - these are now detected and converted to "_"

19-Sept-2014:
 - Adds SVR, improves SVM support

24-Sept-2014:
 - For voxelwise analyses, no longer crashes if voxel dimensions vary between subjects (images resliced to match first using nii_reslice_target)
 - Added nii_mat2nii that extracts voxelwise image(s) from mat files (reversing nii_nii2mat)
 - Previous versions would crash if asked to asked to transform data with z-skew > 1.96 and negative predictors (data translated to be positive prior to sqrt transform)

17-Oct-2014
 - Uploaded to github

##### Installation

There are two methods to install this software:

 - Press the green "Clone or Download" button on the [GitHub page](https://github.com/neurolabusc/NiiStat). Extract the downloaded zip file. This is the easiest method to install the software
 - From the command line, type "git clone https://github.com/neurolabusc/NiiStat.git" to create a repository. This requires you to have the git software installed (it comes with OSX and Linux, but will have to be manually installed fore Windows-based PCs). The advantage of this method is that the software will automatically keep up to date.

In either case, you will want to make sure that the "NiiStat" folder that you installed is in your Matlab path.

##### Notes

 - [NiiStat Wiki](https://www.nitrc.org/plugins/mwiki/index.php/niistat:MainPage)
 - [NiiStat Tutorial](https://www.nitrc.org/plugins/mwiki/index.php/niistat:TutorialPage)
