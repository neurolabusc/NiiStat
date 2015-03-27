##### Introduction

NiiStat is a set of Matlab scripts for analyzing neuroimaging data from clinical populations

This GitHub archive contains the latest release. For usage see the [NITRC website](http://www.nitrc.org/projects/niistat/)

##### Versions

3Mar2015 : 
 - Added modalities ALF, fMRI
 
29Aug2014 : 
 - Fixed bug where nii_stat would crash when 3 conditions were met: (1) conducting correlation analysis while (2) asked to remove white matter regions from (3) ROI that does not include any white matter regions

3Sept2014: 
 - New modality: dtifc is the fiber count for diffusion tractography (while dti is the density)

8Sept2014:  
 - Previous versions did not allow .XLS files where Excel interpreted the filename as a number (e.g. "748" failed, "LM748" was fine)
 - New feature: creates mean intensity map for ROI analyses (similar to sum maps generated for voxelwise analysis).
 - nii_nii2mat will binarize images to 0 or 1 when the modality "lesion" is selected (helps when lesion criteria vary, e.g. for some files 0/1 and others 0/255)
 - nii_nii2mat now adds support for .nii.gz and .voi images (in addition to .hdr and .nii)
 - Excel worksheets no longer crash nii_stat if they have behavioral variables with a "/" or "\" in their names - these are now detected and converted to "_"

19Sept2014:  
 - Adds SVR, improves SVM support

24Sept2014:  
 - For voxelwise analyses, no longer crashes if voxel dimensions vary between subjects (images resliced to match first using nii_reslice_target)
 - Added nii_mat2nii that extracts voxelwise image(s) from mat files (reversing nii_nii2mat)
 - Previous versions would crash if asked to asked to transform data with z-skew > 1.96 and negative predictors (data translated to be positive prior to sqrt transform)

17Oct2014
 - Uploaded to github