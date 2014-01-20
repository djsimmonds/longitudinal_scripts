## extract mean FA values for each tract across subjects

## data path
path=/Volumes/Governator/DTI_STUDY/tbss/stats/
atlaspath=/usr/local/fsl/data/atlases/JHU/

cd $path

## take atlas labels and mask with skeleton
fslmaths $atlaspath/JHU-WhiteMatter-labels-1mm.nii.gz -mas mean_FA_skeleton_mask.nii.gz JHU-WhiteMatter-labels-1mm_skeleton.nii.gz

## make individual tract masks
mkdir skel_tracts_1mm
mv JHU-WhiteMatter-labels-1mm_skeleton.nii.gz skel_tracts_1mm/
cd skel_tracts_1mm

tracts=$(echo 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48)
for t in $tracts; do fslmaths JHU-WhiteMatter-labels-1mm_skeleton.nii.gz -thr $t -uthr $t ${t}.nii.gz; done

## extract mean fa values for each tract for all subjects
cd ..
mkdir tracts_FA

for t in $tracts; do fslstats -t all_FA_skeletonised.nii.gz -k skel_tracts_1mm/${t}.nii.gz -M >> tracts_FA/${t}.txt; done

## code in R to combine text files into one matrix file and save
#path.data<-"/mnt/Schwarzenagger/Governator/DTI_STUDY/dtitk/icbm_dti81/tbss/stats/tracts_FA"
#path.study<-"/home/danisimmonds/Dani/dti_0511"
#files<-dir(path.data)
#temp<-scan(paste(path.data,files[1],sep="/"))
#data.tract<-matrix(NA,length(temp),length(files))
#for(i in 1:length(files)) data.tract[,i]<-scan(paste(path.data,files[i],sep="/"))
#save(data.tract,file=paste(path.study,"data.tract.Rframe",sep="/"))

