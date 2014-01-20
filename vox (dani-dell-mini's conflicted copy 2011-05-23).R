## required libraries
library(Rniftilib) ## for reading and writing nifti files

## paths
paths<-list()
paths$study<-"/home/danisimmonds/Dani/dti_0511"
paths$data<-"/mnt/Schwarzenagger/Governator/DTI_STUDY/dtitk/icbm_dti81/tbss/stats"
paths$labels<-"/usr/share/fsl/data/atlases"

## skeleton mask subsampled to 2mm, thresholded at 0.2 and binarized (winds up ~double volume of 1mm mask)
mask<-nifti.image.read(paste(paths$data,"mean_FA_skeleton_mask_2mm_thr02_bin.nii.gz",sep="/"))

## labels
## labels from template
labels.txt<-list()
labels1<-nifti.image.read(paste(paths$labels,"JHU/JHU-WhiteMatter-labels-2mm.nii.gz",sep="/"))
labels.txt[[1]]<-c(
	"Middle cerebellar peduncle",
	"Pontine crossing tract (a part of MCP)",
	"Genu of corpus callosum",
	"Body of corpus callosum",
	"Splenium of corpus callosum",
	"Fornix (column and body of fornix)",
	"Corticospinal tract R",
	"Corticospinal tract L",
	"Medial lemniscus R",
	"Medial lemniscus L",
	"Inferior cerebellar peduncle R",
	"Inferior cerebellar peduncle L",
	"Superior cerebellar peduncle R",
	"Superior cerebellar peduncle L",
	"Cerebral peduncle R",
	"Cerebral peduncle L",
	"Anterior limb of internal capsule R",
	"Anterior limb of internal capsule L",
	"Posterior limb of internal capsule R",
	"Posterior limb of internal capsule L",
	"Retrolenticular part of internal capsule R",
	"Retrolenticular part of internal capsule L",
	"Anterior corona radiata R",
	"Anterior corona radiata L",
	"Superior corona radiata R",
	"Superior corona radiata L",
	"Posterior corona radiata R",
	"Posterior corona radiata L",
	"Posterior thalamic radiation (include optic radiation) R",
	"Posterior thalamic radiation (include optic radiation) L",
	"Sagittal stratum (include inferior longitidinal fasciculus and inferior fronto-occipital fasciculus) R",
	"Sagittal stratum (include inferior longitidinal fasciculus and inferior fronto-occipital fasciculus) L",
	"External capsule R",
	"External capsule L",
	"Cingulum (cingulate gyrus) R",
	"Cingulum (cingulate gyrus) L",
	"Cingulum (hippocampus) R",
	"Cingulum (hippocampus) L",
	"Fornix (cres) / Stria terminalis (can not be resolved with current resolution) R",
	"Fornix (cres) / Stria terminalis (can not be resolved with current resolution) L",
	"Superior longitudinal fasciculus R",
	"Superior longitudinal fasciculus L",
	"Superior fronto-occipital fasciculus (could be a part of anterior internal capsule) R",
	"Superior fronto-occipital fasciculus (could be a part of anterior internal capsule) L",
	"Uncinate fasciculus R",
	"Uncinate fasciculus L",
	"Tapetum R",
	"Tapetum L"
)
## sub-labels
tract.type<-c("cer","cer","cal","cal","cal","assL","pro","pro","pro","pro","cer","cer","cer","cer","pro","pro","pro","pro","pro","pro","pro","pro","pro","pro","pro","pro","pro","pro","pro","pro","ass","ass","ass","ass","assL","assL","assL","assL","assL","assL","ass","ass","ass","ass","ass","ass","cal","cal")

## additional labels from harvard/MGH unthresholded probabilistic gray matter template to identify voxels which are not in the JHU template via their proximity to gray matter (first will try cortical, then subcortical) --> this is really a hack, a better way would probably be tractography to identify, but this would take a long time
labels2<-nifti.image.read(paste(paths$labels,"HarvardOxford/HarvardOxford-cort-maxprob-thr0-2mm.nii.gz",sep="/"))
labels.txt[[2]]<-c(
	"Frontal Pole",
	"Insular Cortex",
	"Superior Frontal Gyrus",
	"Middle Frontal Gyrus",
	"Inferior Frontal Gyrus, pars triangularis",
	"Inferior Frontal Gyrus, pars opercularis",
	"Precentral Gyrus",
	"Temporal Pole",
	"Superior Temporal Gyrus, anterior division",
	"Superior Temporal Gyrus, posterior division",
	"Middle Temporal Gyrus, anterior division",
	"Middle Temporal Gyrus, posterior division",
	"Middle Temporal Gyrus, temporo-occipital part",
	"Inferior Temporal Gyrus, anterior division",
	"Inferior Temporal Gyrus, posterior division",
	"Inferior Temporal Gyrus, temporooccipital part",
	"Postcentral Gyrus",
	"Superior Parietal Lobule",
	"Supramarginal Gyrus, anterior division",
	"Supramarginal Gyrus, posterior division",
	"Angular Gyrus",
	"Lateral Occipital Cortex, superior division",
	"Lateral Occipital Cortex, inferior division",
	"Intracalcarine Cortex",
	"Frontal Medial Cortex",
	"Juxtapositional Lobule Cortex (formerly Supplementary Motor Cortex)",
	"Subcallosal Cortex",
	"Paracingulate Gyrus",
	"Cingulate Gyrus, anterior division",
	"Cingulate Gyrus, posterior division",
	"Precuneous Cortex",
	"Cuneal Cortex",
	"Frontal Orbital Cortex",
	"Parahippocampal Gyrus, anterior division",
	"Parahippocampal Gyrus, posterior division",
	"Lingual Gyrus",
	"Temporal Fusiform Cortex, anterior division",
	"Temporal Fusiform Cortex, posterior division",
	"Temporal Occipital Fusiform Cortex",
	"Occipital Fusiform Gyrus",
	"Frontal Operculum Cortex",
	"Central Opercular Cortex",
	"Parietal Operculum Cortex",
	"Planum Polare",
	"Heschl's Gyrus (includes H1 and H2)",
	"Planum Temporale",
	"Supracalcarine Cortex",
	"Occipital Pole"
)
labels3<-nifti.image.read(paste(paths$labels,"HarvardOxford/HarvardOxford-sub-maxprob-thr0-2mm.nii.gz",sep="/"))
## for some reason, some numbers are skipped for the labels, so need to adjust
labels3.ind=c(2:4,10:13,16:18,26,41:43,49:54,58)
labels.txt[[3]]=c(
	"Left Cerebral White Matter",
	"Left Cerebral Cortex",
	"Left Lateral Ventrical",
	"Left Thalamus",
	"Left Caudate",
	"Left Putamen",
	"Left Pallidum",
	"Brain-Stem",
	"Left Hippocampus",
	"Left Amygdala",
	"Left Accumbens",
	"Right Cerebral White Matter",
	"Right Cerebral Cortex",
	"Right Lateral Ventricle",
	"Right Thalamus",
	"Right Caudate",
	"Right Putamen",
	"Right Pallidum",
	"Right Hippocampus",
	"Right Amygdala",
	"Right Accumbens"
)
## getting voxels for analysis
temp<-which(mask[,,]>0)
vox<-data.frame(
	ind=temp,
	i=numeric(length(temp)),
	j=numeric(length(temp)),
	k=numeric(length(temp)),
	x=numeric(length(temp)),
	y=numeric(length(temp)),
	z=numeric(length(temp))
)
atlas=character(length(temp))
label=character(length(temp))

## log
sink("log.vox.R")
for(i in 1:length(vox$ind)){
	cat(i,"")
	vox$i[i]<-vox$ind[i]%%dim(mask)[1]
	vox$j[i]<-(vox$ind[i]%%(dim(mask)[1]*dim(mask)[2]))%/%dim(mask)[1]+1
	vox$k[i]<-vox$ind[i]%/%(dim(mask)[1]*dim(mask)[2])+1
	vox$x[i]<-mask$qto.xyz[1,4]+vox$i[i]*mask$qto.xyz[1,1]
	vox$y[i]<-mask$qto.xyz[2,4]+vox$j[i]*mask$qto.xyz[2,2]
	vox$z[i]<-mask$qto.xyz[3,4]+vox$k[i]*mask$qto.xyz[3,3]
	if(labels1[vox$i[i],vox$j[i],vox$k[i]]>0){
		atlas[i]<-"JHU"
		label[i]<-labels.txt[[1]][labels1[vox$i[i],vox$j[i],vox$k[i]]]
	}else if(labels2[vox$i[i],vox$j[i],vox$k[i]]>0){
		atlas[i]<-"MGH_cort"
		label[i]<-labels.txt[[2]][labels2[vox$i[i],vox$j[i],vox$k[i]]]
	}else if(labels3[vox$i[i],vox$j[i],vox$k[i]]>0){
		atlas[i]<-"MGH_sub"
		label[i]<-labels.txt[[3]][which(labels3.ind==labels3[vox$i[i],vox$j[i],vox$k[i]])]
	}else{
		atlas[i]<-"NA"
		label[i]<-"NA"
	}
}
sink()
vox$ind<-as.factor(vox$ind)
vox$atlas<-as.factor(atlas)
vox$label<-as.factor(label)

## saves file with voxel indices
save(vox,file=paste(paths$study,"vox.Rframe",sep="/"))

## identify tracts for analysis
tracts<-list()
for(i in 1:length(labels.txt[[1]])) tracts[[i]]<-which(vox$label==labels.txt[[1]][i])
save(tracts,file=paste(paths$study,"tracts.Rframe",sep="/"))
