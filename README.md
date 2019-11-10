# NiiStat     
A version of NiiStat which runs entirely in Octave, a free version of MATLAB.  

## Introduction
NiiStat is a set of MATLAB/Octave scripts for analyzing neuroimaging data from clinical populations.       
Note, the code in this repo will only work in Octave.

This GitHub archive contains the latest release for Octave. For usage see the [NITRC website](http://www.nitrc.org/projects/NiiStat/).
     
## Installation     
1) Download Octave 4.4.1 and install it     
MacOS: https://octave-app.org/Download.html     
Linux: https://www.gnu.org/software/octave/download.html     
Windows: https://ftpmirror.gnu.org/octave/windows/     
Source: https://ftpmirror.gnu.org/octave        
          
2) Run Octave and install the io package       
In the Octave Command Window, run the following command:       
+ `pkg install -forge io`       
          
3) Download NiiStat     
Pick __one__ of the following methods:       
a) Press the green "Clone or Download" button on the [GitHub](https://github.com/AnthonyAndroulakis/NiiStat) page. Extract the downloaded zip file. Renamed the extracted folder (NiiStat-master) to NiiStat. This is the easiest method to install the software.     
__OR__       
b) From the command line, type `git clone https://github.com/AnthonyAndroulakis/NiiStat`. This requires you to have the git software installed. The advantage of this method is that the software will be kept up-to-date automatically.      
         
4) Download SPM and Compile it     
In the Octave Command Line,        
__cd into the NiiStat folder (that you just downloaded)__      
__Download SPM12 r4787__       
+ `unzip('https://github.com/spm/spm12/archive/r7487.zip',pwd);`         
__Patch SPM12__      
+ `urlwrite('https://raw.githubusercontent.com/spm/spm-docker/master/octave/spm12_r7487.patch','spm12_r7487.patch');`     
+ `system('patch -p3 -d spm12-r7487 < spm12_r7487.patch');`      
__Compile MEX files__      
+ `cd spm12-r7487/src`      
+ `system('make PLATFORM=octave');`      
+ `system('make PLATFORM=octave install');`      
__Go back 2 levels to access NiiStat__         
+ `cd ../../`     
           
5) Run NiiStat
+ `NiiStat`     

## Versions

9-October-2016 (MATLAB):
 - Saves node/edge map for connectome analyses (DTI/Rest, GLM/SVM)
 - Supports most NPM val files directly, unifies logical mask code

30-June-2016 (MATLAB):
 - Main program name changed from "nii_stat" to "NiiStat" to match NITRC
 - Behavioral data should be in excel spreadsheet tab named "NiiStat" (old format expected "Data (2)", which still works).
 - Fix bug specific to 6June2016 version where region of interest maps were drawn scrambled.
 - Script automatically checks GitHub for updates.

3-Mar-2015 (MATLAB):
 - Added modalities ALF, fMRI

29-Aug-2014 (MATLAB):
 - Fixed bug where nii_stat would crash when 3 conditions were met: (1) conducting correlation analysis while (2) asked to remove white matter regions from (3) ROI that does not include any white matter regions

3-Sept-2014 (MATLAB)
 - New modality: dtifc is the fiber count for diffusion tractography (while dti is the density)

8-Sept-2014 (MATLAB):
 - Previous versions did not allow .XLS files where Excel interpreted the filename as a number (e.g. "748" failed, "LM748" was fine)
 - New feature: creates mean intensity map for ROI analyses (similar to sum maps generated for voxelwise analysis).
 - nii_nii2mat will binarize images to 0 or 1 when the modality "lesion" is selected (helps when lesion criteria vary, e.g. for some files 0/1 and others 0/255)
 - nii_nii2mat now adds support for .nii.gz and .voi images (in addition to .hdr and .nii)
 - Excel worksheets no longer crash nii_stat if they have behavioral variables with a "/" or "\" in their names - these are now detected and converted to "_"

19-Sept-2014 (MATLAB):
 - Adds SVR, improves SVM support

24-Sept-2014 (MATLAB):
 - For voxelwise analyses, no longer crashes if voxel dimensions vary between subjects (images resliced to match first using nii_reslice_target)
 - Added nii_mat2nii that extracts voxelwise image(s) from mat files (reversing nii_nii2mat)
 - Previous versions would crash if asked to asked to transform data with z-skew > 1.96 and negative predictors (data translated to be positive prior to sqrt transform)

10-Nov-2019 (Octave):
 - edited NiiStat code to work in Octave, a free version of MATLAB

## Notes

 - [NiiStat Wiki](https://www.nitrc.org/plugins/mwiki/index.php/niistat:MainPage)
 - [NiiStat Tutorial](https://www.nitrc.org/plugins/mwiki/index.php/niistat:TutorialPage)
 
 ###### This NiiStat version was last edited by Anthony Androulakis on November 10, 2019.
