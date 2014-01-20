#!/bin/bash

set -e #stop if any error occurs
set -x #print out progress on the script

# directory where data is, organized .../luna_id/scan_date/dti, with dicom (raw) files in dti directory. within parent, should only be subject folders
cd /Volumes/Governator/DTI_STUDY/subjects
#mkdir /Volumes/Governator/DTI_STUDY/tbss # will put FA files in this folder for later TBSS analysis

##################
# SINGLE SUBJECT #
##################

# pulls all subject folders, then loops analysis for each subject
subjects=$( find $PWD -mindepth 1 -maxdepth 1 -type d )
for subject in $subjects
do
    cd $subject
    
    # extracts luna id from data path
    lunaID=$( echo $subject | cut -d "/" -f 6 )
	echo $lunaID
    # finds all scans for subject, then loops analysis for each scan
    scans=$( find $PWD -mindepth 1 -maxdepth 1 -type d )
    for scan in $scans
    do
        cd $scan
        
        # extracts scan date for each scan
        scanDate=$( echo $scan | cut -d "/" -f 7 )
        
        # create preprocessing and analysis directories
        if [ ! -d $scan/preprocessing ]; then
        	mkdir $scan/preprocessing
        fi
        
        if [ ! -d $scan/analysis ]; then
        	mkdir $scan/analysis
        fi
        
        # converts dicoms to fsl zipped nifti format and moves image and bval/bvec files to preprocessing folder
        cd "dti"
        dcm2nii *.dcm
        echo ${lunaID}_${scanDate}
        mv *.nii.gz ../preprocessing/${lunaID}_${scanDate}.nii.gz
        mv *.bvec ../preprocessing/${lunaID}_${scanDate}.bvec
        mv *.bval ../preprocessing/${lunaID}_${scanDate}.bval
        cd ../preprocessing

        #################
        # PREPROCESSING #
        #################

        # skull stripping, also creates a brain mask file
        bet ${lunaID}_${scanDate}.nii.gz ${lunaID}_${scanDate}_brain.nii.gz -f 0.3 -m

        # correct for eddy current distortion and motion
        eddy_correct ${lunaID}_${scanDate}.nii.gz ${lunaID}_${scanDate}_eddy.nii.gz 0

        # should rotate vectors to account for motion correction; however, this is commented out as most people (including susumu mori) say it doesn't matter as long as the rotations aren't huge (and maybe you shouldn't be using your data if they are)
        # rotate_bvecs ${lunaID}_${scanDate}_eddy.ecclog ${lunaID}_${scanDate}.bvec
        # mkdir rotate_bvecs
        # mv *.mat ${lunaID}_${scanDate}.bvec_old rotate_bvecs

        ###################
        # TENSOR ANALYSIS #
        ###################

        # tensor calculations
        dtifit --data=${lunaID}_${scanDate}_eddy.nii.gz --out=../analysis/${lunaID}_${scanDate} --mask=${lunaID}_${scanDate}_brain_mask.nii.gz --bvecs=${lunaID}_${scanDate}.bvec --bvals=${lunaID}_${scanDate}.bval

        # calculates radial diffusivity, since this is not done by dtifit
        cd ../analysis
        fslmaths ${lunaID}_${scanDate}_L2.nii.gz -add ${lunaID}_${scanDate}_L3.nii.gz -div 2 ${lunaID}_${scanDate}_RD.nii.gz

        cp *FA.nii.gz /Volumes/Governator/DTI_STUDY/tbss

    done

done

#########################
# GROUP ANALYSIS - TBSS #
#########################
