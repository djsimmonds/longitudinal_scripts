## required libraries
library(Rniftilib) ## for reading and writing nifti files

## paths
paths<-list()
paths$study<-"/home/danisimmonds/Dani/dti_0511/tbss"
paths$data<-"/mnt/Schwarzenagger/Governator/DTI_STUDY/tbss/stats"
paths$labels<-"/usr/share/fsl/data/atlases"

## skeleton mask subsampled to 2mm, thresholded at 0.2 and binarized (winds up ~double volume of 1mm mask)
mask<-nifti.image.read(paste(paths$data,"mean_FA_skeleton_mask_2mm_thr02_bin.nii.gz",sep="/"))
d<-mask$dim
c<-mask$qto.xyz

## labels
## labels from template
l1.img<-nifti.image.read(paste(paths$labels,"JHU/JHU-WhiteMatter-labels-2mm.nii.gz",sep="/"))
l1.atlas<-c(
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
l1.type<-c("Ce","Ce","Ca","Ca","Ca","AL","P","P","P","P","Ce","Ce","Ce","Ce","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","A","A","A","A","AL","AL","AL","AL","AL","AL","A","A","A","A","A","A","Ca","Ca")
l1.tract<-c("MCP","PCT","CCG","CCB","CCS","FOR","CST","CST","ML","ML","ICP","ICP","SCP","SCP","CP","CP","IC","IC","IC","IC","IC","IC","CR","CR","CR","CR","CR","CR","PTR","PTR","SS","SS","EC","EC","CG","CG","CG","CG","FOR","FOR","SLF","SLF","SFOF","SFOF","UF","UF","CCT","CCT")
l1.sub<-c("<NA>","<NA>","<NA>","<NA>","<NA>","CB","<NA>","<NA>","<NA>","<NA>","<NA>","<NA>","<NA>","<NA>","<NA>","<NA>","A","A","P","P","R","R","A","A","S","S","P","P","<NA>","<NA>","<NA>","<NA>","<NA>","<NA>","CG","CG","H","H","C","C","<NA>","<NA>","<NA>","<NA>","<NA>","<NA>","<NA>","<NA>")

## additional labels from harvard/MGH unthresholded probabilistic gray matter template to identify voxels which are not in the JHU template via their proximity to gray matter (first will try cortical, then subcortical) --> this is really a hack, a better way would probably be tractography to identify, but this would take a long time
l2.img<-nifti.image.read(paste(paths$labels,"HarvardOxford/HarvardOxford-cort-maxprob-thr0-2mm.nii.gz",sep="/"))
l2.atlas<-c(
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
l2.type<-c("F","F","F","F","F","F","SM","T","T","T","T","T","T","T","T","T","SM","P","P","P","P","O","O","O","F","F","F","F","F","P","P","O","F","T","T","O","T","T","T","O","F","SM","P","T","T","T","O","O")


l3.img<-nifti.image.read(paste(paths$labels,"HarvardOxford/HarvardOxford-sub-maxprob-thr0-2mm.nii.gz",sep="/"))
l3.atlas<-c(
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
## for some reason, some numbers are skipped for the labels, so need to adjust
l3.ind<-c(2:4,10:13,16:18,26,41:43,49:54,58)
l3.type<-c("C","C","V","T","BG","BG","BG","B","MT","MT","BG","C","C","V","T","BG","BG","BG","MT","MT","BG")

## additional labels from harvard/MGH unthresholded probabilistic gray matter template to identify voxels which are not in the JHU template via their proximity to gray matter (first will try cortical, then subcortical) --> this is really a hack, a better way would probably be tractography to identify, but this would take a long time
l4.img<-nifti.image.read(paste(paths$labels,"Cerebellum/Cerebellum-MNIfnirt-maxprob-thr0-2mm.nii.gz",sep="/"))
l4.atlas<-c(
	"Left I-IV",
	"Right I-IV",
	"Left V",
	"Right V",
	"Left VI",
	"Vermis VI",
	"Right VI",
	"Left Crus I",
	"Vermis Crus I",
	"Right Crus I",
	"Left Crus II",
	"Vermis Crus II",
	"Right Crus II",
	"Left VIIb",
	"Vermis VIIb",
	"Right VIIb",
	"Left VIIIa",
	"Vermis VIIIa",
	"Right VIIIa",
	"Left VIIIb",
	"Vermis VIIIb",
	"Right VIIIb",
	"Left IX",
	"Vermis IX",
	"Right IX",
	"Left X",
	"Vermis X",
	"Right X"
)
l4.type<-c("C","C","C","C","C","V","C","C","V","C","C","V","C","C","V","C","C","V","C","C","V","C","C","V","C","C","V","C")
l4.sub<-c("I-IV","I_IV","V","V","VI","VI","VI","Crus I","Crus I","Crus I","Crus II","Crus II","Crus II","VIIb","VIIb","VIIb","VIIIa","VIIIa","VIIIa","VIIIb","VIIIb","VIIIb","IX","IX","IX","X","X","X")

## getting voxels for analysis
ind<-which(mask[,,]>0)

## voxel info
i<-numeric(length(ind))
j<-numeric(length(ind))
k<-numeric(length(ind))
x<-numeric(length(ind))
y<-numeric(length(ind))
z<-numeric(length(ind))
atlas<-character(length(ind))
L1.ind<-numeric(length(ind))
L1.atlas<-character(length(ind))
L1.type<-character(length(ind))
L1.tract<-character(length(ind))
L1.sub<-character(length(ind))
L2.ind<-numeric(length(ind))
L2.atlas<-character(length(ind))
L2.type<-character(length(ind))
L3.ind<-numeric(length(ind))
L3.ind.<-numeric(length(ind))
L3.atlas<-character(length(ind))
L3.type<-character(length(ind))
L4.ind<-numeric(length(ind))
L4.atlas<-character(length(ind))
L4.type<-character(length(ind))
L4.sub<-character(length(ind))

## log
sink("log.vox.R")

for(a in 1:length(ind)){
	cat(a,"")

	## voxel indices
	i[a]<-ind[a]%%d[1]
	j[a]<-(ind[a]%%(d[1]*d[2]))%/%d[1]+1
	k[a]<-ind[a]%/%(d[1]*d[2])+1
	x[a]<-c[1,4]+i[a]*c[1,1]
	y[a]<-c[2,4]+j[a]*c[2,2]
	z[a]<-c[3,4]+k[a]*c[3,3]

	## labels
	if(l1.img[i[a],j[a],k[a]]>0){
		atlas[a]<-"JHU"
		L1.ind[a]<-l1.img[i[a],j[a],k[a]]
		L1.atlas[a]<-l1.atlas[L1.ind[a]]
		L1.type[a]<-l1.type[L1.ind[a]]
		L1.tract[a]<-l1.tract[L1.ind[a]]
		L1.sub[a]<-l1.sub[L1.ind[a]]

	}else if(l2.img[i[a],j[a],k[a]]>0){
		atlas[a]<-"MGH_cort"
		L2.ind[a]<-l2.img[i[a],j[a],k[a]]
		L2.atlas[a]<-l2.atlas[L2.ind[a]]
		L2.type[a]<-l2.type[L2.ind[a]]

	}else if(l3.img[i[a],j[a],k[a]]>0){
		atlas[a]<-"MGH_sub"
		L3.ind[a]<-l3.img[i[a],j[a],k[a]]
		L3.ind.[a]<-which(l3.ind==L3.ind[a])
		L3.atlas[a]<-l3.atlas[L3.ind.[a]]
		L3.type[a]<-l3.type[L3.ind.[a]]

	}else if(l4.img[i[a],j[a],k[a]]>0){
		atlas[a]<-"Cerebellum"
		L4.ind[a]<-l4.img[i[a],j[a],k[a]]
		L4.atlas[a]<-l4.atlas[L4.ind[a]]
		L4.type[a]<-l4.type[L4.ind[a]]
		L4.sub[a]<-l4.sub[L4.ind[a]]

	}else{
		atlas[a]<-"NA"
	}
}

sink()

l2d<-function(lab,len=length(atlas)){
	len.<-length(lab)
	if(len.<len) lab[(len.+1):len]<-NA
	lab
}

vox<-data.frame(
	ind=as.factor(ind),
	i,j,k,x,y,z,
	atlas,
	L1.ind<-l2d(L1.ind),
	L1.atlas<-l2d(L1.atlas),
	L1.type<-l2d(L1.type),
	L1.tract<-l2d(L1.tract),
	L1.sub<-l2d(L1.sub),
	L2.ind<-l2d(L2.ind),
	L2.atlas<-l2d(L2.atlas),
	L2.type<-l2d(L2.type),
	L3.ind<-l2d(L3.ind),
	L3.atlas<-l2d(L3.atlas),
	L3.type<-l2d(L3.type),
	L4.ind<-l2d(L4.ind),
	L4.atlas<-l2d(L4.atlas),
	L4.type<-l2d(L4.type),
	L4.sub<-l2d(L4.sub)
)

## saves file with voxel indices
save(vox,file=paste(paths$study,"vox.Rframe",sep="/"))

## identify tracts/rois for analysis
atlas.ind<-list()

## L1.type
atlas.ind$L1.type<-list()
uniq<-unique(L1.type)
for(u in 1:length(uniq)) atlas.ind$L1.type[[u]]<-which(L1.type==uniq[u])
names(atlas.ind$L1.type)<-uniq

## L1.tract
atlas.ind$L1.tract<-list()
uniq<-unique(L1.tract)
for(u in 1:length(uniq)) atlas.ind$L1.tract[[u]]<-which(L1.tract==uniq[u])
names(atlas.ind$L1.tract)<-uniq

## L1.sub
atlas.ind$L1.sub<-list()
uniq<-unique(L1.sub)
for(u in 1:length(uniq)) atlas.ind$L1.sub[[u]]<-which(L1.sub==uniq[u])
names(atlas.ind$L1.sub)<-uniq

## L2.atlas
atlas.ind$L2.atlas<-list()
uniq<-unique(L2.atlas)
for(u in 1:length(uniq)) atlas.ind$L2.atlas[[u]]<-which(L2.atlas==uniq[u])
names(atlas.ind$L2.atlas)<-uniq

## L2.type
atlas.ind$L2.type<-list()
uniq<-unique(L2.type)
for(u in 1:length(uniq)) atlas.ind$L2.type[[u]]<-which(L2.type==uniq[u])
names(atlas.ind$L2.type)<-uniq

## L3.type
atlas.ind$L3.type<-list()
uniq<-unique(L3.type)
for(u in 1:length(uniq)) atlas.ind$L3.type[[u]]<-which(L3.type==uniq[u])
names(atlas.ind$L3.type)<-uniq

## L4.type
atlas.ind$L4.type<-list()
uniq<-unique(L4.type)
for(u in 1:length(uniq)) atlas.ind$L4.type[[u]]<-which(L4.type==uniq[u])
names(atlas.ind$L4.type)<-uniq

## L4.sub
atlas.ind$L4.sub<-list()
uniq<-unique(L4.sub)
for(u in 1:length(uniq)) atlas.ind$L4.sub[[u]]<-which(L4.sub==uniq[u])
names(atlas.ind$L4.sub)<-uniq

save(atlas.ind,file=paste(paths$study,"atlas.ind.Rframe",sep="/"))
