## paths
paths<-list()
paths$data<-"/home/danisimmonds/Dani/dti_0511/tbss"
paths$scripts<-"/home/danisimmonds/Dropbox/scripts2"
paths$demo<-"/home/danisimmonds/Dani/dti_0511/demographics"
paths$analysis<-paste(paths$data,"analysis","cross","roi",sep="/")

## load data
## voxel-wise data frame
load(paste(paths$data,"data.vox.Rframe",sep="/"))
## demographics frame
load(paste(paths$demo,"demo.Rframe",sep="/"))
## outlier weights
load(paste(paths$analysis,"out_weights",sep="/"))

## subset
subset<-W[[1]]
w.<-W[[2]]

## calculate weighted average of "mature" subjects (age>18)
data.<-Y[subset,]
demo.<-demo[subset,]
mat.ind<-which(demo.$age>18)
data.<-t(t(data.)/as.numeric(colSums(data.[mat.ind,]*w.[mat.ind,])/colSums(w.[mat.ind,])))

#if(type=="long.slope"){
#	data.diff<-diff(data.)
#	demo.diff<-diff(demo.)
#}

## save data in analysis folder
save(data.,file=paste(paths$analysis,"data",sep="/"))
save(demo.,file=paste(paths$analysis,"demo",sep="/"))
