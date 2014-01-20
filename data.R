## required libraries
library(Rniftilib) ## for reading and writing nifti files

## paths
paths<-list()
paths$study<-"/home/danisimmonds/Dani/dti_0511"
paths$data<-"/mnt/Schwarzenagger/Governator/DTI_STUDY/dtitk/icbm_dti81/tbss/stats"

## skeleton mask subsampled to 2mm, thresholded at 0.2 and binarized (winds up ~double volume of 1mm mask)
mask<-nifti.image.read(paste(paths$data,"mean_FA_skeleton_mask_2mm_thr02_bin.nii.gz",sep="/"))

## all_FA_skeletonised from tbss subsampled to 2mm, divided by resampled mask to adjust FA values from resampling, masked, smoothed with sigma=1.7mm (4mm FWHM; ~2x acquisition resolution), divided by smoothed resampled mask to again adjust FA values, remasked
subjs<-nifti.image.read(paste(paths$data,"all_FA_skeletonised_2mm_div_mas_s_fwhm4mm_div_mas.nii.gz",sep="/"))

## load vox.R and tracts.R if script has already been run (run script if not)
#load(paste(paths$study,"vox.Rframe",sep="/"))
#load(paste(paths$study,"tracts.Rframe",sep="/"))

## tract means matrix
data.tract<-matrix(NA,dim(subjs)[4],length(tracts))
for(i in 1:length(tracts)){
	temp<-matrix(NA,dim(subjs)[4],length(tracts[[i]]))
	for(j in 1:length(tracts[[i]])) temp[,j]<-subjs[vox$i[tracts[[i]][j]],vox$j[tracts[[i]][j]],vox$k[tracts[[i]][j]],]
	data.tract[,i]<-rowMeans(temp)
}
save(data.tract,file=paste(paths$study,"data.vox.tract.Rframe",sep="/"))

## voxel matrix
data.vox<-matrix(NA,dim(subjs)[4],dim(vox)[1])
for(i in 1:dim(vox)[1]) data.vox[,i]<-subjs[vox$i[i],vox$j[i],vox$k[i],]
save(data.vox,file=paste(paths$study,"data.vox.Rframe",sep="/"))

