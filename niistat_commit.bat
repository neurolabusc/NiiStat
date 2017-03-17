#!/bin/bash
#go to nii_preprocess folder
if [[ ! -e /home/crlab/NiiStat ]]; then
            mkdir /home/crlab/NiiStat
fi
cd /home/crlab/NiiStat
#commit with message
git add .
d=`date`
h=`hostname`
git commit -a -m "$h-$d"
#git commit -a
#add files in NiiStat folder
#git add .
git push origin master
