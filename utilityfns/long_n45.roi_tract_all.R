## load data/demographics
load("/home/danisimmonds/Dani/dti_0511/demographics/demo.Rframe")
load("/home/danisimmonds/Dani/dti_0511/tbss/data.vox.Rframe")
load("/home/danisimmonds/Dani/dti_0511/tbss/atlas.ind.Rframe")

## order (top to bottom): association/limbic, association, projection, cerebellar, callosal
y.ind<-c(15,21,19,20,7,3,5,4,6,2,9,16,14,17,11,13,22,10,18,8,12)
y.names<-names(atlas.ind[[2]])[y.ind]

## paths
path<-"/home/danisimmonds/Dani/dti_0511/tbss/analysis/long.n45"
path.<-paste(path,"roi_tract_all","ind",sep="/")
path.scripts<-"/home/danisimmonds/Dropbox/scripts2"

## outlier scripts
source(paste(path.scripts,"outlier.R",sep="/"))
source(paste(path.scripts,"out.sub.R",sep="/"))

## index of adults for maturation index calculation
ind<-demo$nscans.tot>=4
adult.ind<-intersect(which(demo$age>18),which(ind))

## convert to maturation index, then average across voxels in tract
Y<-sapply(1:dim(data.vox)[2], function(y) data.vox[,y]/mean(data.vox[adult.ind,y])) 
Y.<-sapply(1:length(y.ind), function(i) rowMeans(Y[,atlas.ind[[2]][[y.ind[i]]]]))

## outlier removal
W.<-outlier(Y.,sub=ind,within=TRUE,allY=TRUE,sd.thr=2,basename=paste(path.,"out",sep="/"))

## revised data/demographics
if(is.null(dim(Y.))) data.<-as.matrix(Y.[W.[[1]]]) else data.<-Y.[W.[[1]],]
save(data.,file=paste(path.,"data",sep="/"))
demo.<-demo[W.[[1]],]
save(demo.,file=paste(path.,"demo",sep="/"))

## load extra scripts (need to rewrite utilityfns.R)
source(paste(path.scripts,"utilityfns","model.setup.R",sep="/"))
source(paste(path.scripts,"utilityfns","pred.str.R",sep="/"))
source(paste(path.scripts,"utilityfns","deriv.check.R",sep="/"))
source(paste(path.scripts,"utilityfns","model.est.R",sep="/"))
source(paste(path.scripts,"utilityfns","fixef.set.R",sep="/"))
source(paste(path.scripts,"utilityfns","coef.est.R",sep="/"))
source(paste(path.scripts,"utilityfns","ll.test.R",sep="/"))
source(paste(path.scripts,"utilityfns","deriv.est.R",sep="/"))

## analysis setup and run analysis
source(paste(path.scripts,"analysis.setup_long.n45.roi_tract_all.ind.R",sep="/"))
source(paste(path.scripts,"analysis.R",sep="/"))

