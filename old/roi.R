## rois
ynames <- c("all","proj","assoc","assoc.limb","cal","cer.c","cer.p","par","occ","sm","fron","temp","mt","bg","thal","ML","PCT","CST","CP","IC.A","IC.P","IC.R","CR.A","CR.S","CR.P","PTR","SFOF","SS","EC","UF","SLF","FOR.CB","FOR.C","CIN.CG","CIN.H","CAL.G","CAL.B","CAL.S","CAL.T","CER.I","CER.M","CER.S")
Y <- sapply(1:length(ynames), function(y) scan(paste("/mnt/connor/DTI_STUDY/tbss0711/roi_data_FA",ynames[y],sep="/")))
save(Y,file="/home/danisimmonds/Dani/dti_0811/analysis/roi/Y")

## L/R
Y <- sapply(1:length(ynames), function(y) c(scan(paste("/mnt/connor/DTI_STUDY/tbss0711/roi_data_FA/",ynames[y],".R",sep="")),scan(paste("/mnt/connor/DTI_STUDY/tbss0711/roi_data_FA/",ynames[y],".L",sep=""))))
save(Y,file="/home/danisimmonds/Dani/dti_0811/analysis/roiLR/Y")

## vox
## !#
fslmaths mean_FA_skeleton_mask -subsamp2 mean_FA_skeleton_mask_2mm
fslmaths mean_FA_skeleton_mask_2mm -thr 0.2 -bin mean_FA_skeleton_mask_2mm_bin
fslmaths mean_FA_skeleton_mask_2mm_bin -s 1.7 mean_FA_skeleton_mask_2mm_bin_s
fslmaths all_FA_skeletonised -subsamp2 -div mean_FA_skeleton_mask_2mm -mas mean_FA_skeleton_mask_2mm_bin all_FA_skeletonised_2mm
fslmaths all_FA_skeletonised_2mm -s 1.7 -div mean_FA_skeleton_mask_2mm_bin_s -mas mean_FA_skeleton_mask_2mm_bin all_FA_skeletonised_2mm_s

