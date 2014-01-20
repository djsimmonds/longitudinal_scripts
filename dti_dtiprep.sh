#!/bin/bash

set -e #stop if any error occurs
set -x #print out progress on the script

# directory where data is, organized .../luna_id/scan_date/dti, with dicom (raw) files in dti directory. within parent, should only be subject folders
cd /Volumes/Governator/DTI_STUDY/subjects
#mkdir /Volumes/Governator/DTI_STUDY/tbss_dtiprep # will put FA files in this folder for later TBSS analysis

##################
# SINGLE SUBJECT #
##################

# pulls all subject folders, then loops analysis for each subject
subjects=$( find $PWD -mindepth 1 -maxdepth 1 -type d )
for subject in $subjects
do
    cd $subject
    
    # extracts luna id from data path
    lunaID=$( echo $subject | cut -d "/" -f 7 )
	echo $lunaID
	
    # finds all scans for subject, then loops analysis for each scan
    scans=$( find $PWD -mindepth 1 -maxdepth 1 -type d )
    for scan in $scans
    do
        cd $scan
        
        # extracts scan date for each scan
        scanDate=$( echo $scan | cut -d "/" -f 8 )
        
        # create directory for dtiprep
        # only runs if directory doesn't exist, so if you want to rerun an analysis, delete this folder
		if [ ! -d $scan/dtiprep ]; then
			mkdir $scan/dtiprep
			
			# converts dicoms to nrrd format and runs DTIPrep pipeline
			~/Dani/DTIPrep/DTIPrep-build/bin/DicomToNrrdConverter --inputDicomDirectory $scan/dti --outputDirectory $scan/dtiprep --outputVolume ${lunaID}_${scanDate}.nrrd
			cd $scan/dtiprep
			~/Dani/DTIPrep/DTIPrep-build/bin/DTIPrep --DWINrrdFile *.nrrd --xmlProtocol ${lunaID}_${scanDate}.xml --default --check --outputFolder .
		fi
		
    done
    
done
