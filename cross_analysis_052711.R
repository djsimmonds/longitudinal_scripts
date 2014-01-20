load("/home/danisimmonds/Dani/dti_0511/demographics/demo.Rframe")
load("/home/danisimmonds/Dani/dti_0511/tbss/data.vox.Rframe")
load("/home/danisimmonds/Dani/dti_0511/tbss/atlas.ind.Rframe")
#path<-"/home/danisimmonds/Dani/dti_0511/tbss/analysis/cross"
#path<-"/home/danisimmonds/Dani/dti_0511/tbss/analysis/long.n45"
path<-"/home/danisimmonds/Dani/dti_0511/tbss/analysis/long.all"
path.scripts<-"/home/danisimmonds/Dropbox/scripts2"

## outlier scripts
source(paste(path.scripts,"outlier.R",sep="/"))
source(paste(path.scripts,"out.sub.R",sep="/"))

## index of adults for maturation index calculation
#ind<-demo$nscans.ind==1
#ind<-demo$nscans.tot>=4
ind<-!logical(dim(demo)[[1]])
adult.ind<-intersect(which(demo$age>18),which(ind))

## convert to maturation index
Y<-sapply(1:dim(data.vox)[2], function(y) data.vox[,y]/mean(data.vox[adult.ind,y]))


## all1/ind
path.<-paste(path,"roi_all1","ind",sep="/")
## average across all voxels in skeleton
#Y.<-rowMeans(Y)

## all16/all ## switched to just 14
path.<-paste(path,"roi_all16","ind",sep="/")
#y.ind<-cbind(c(1,1,1,1,1,5,5,5,5,5,6,6,6,6,6,7),c(2,3,4,5,6,1,3,4,5,6,2,3,4,5,6,2))
y.ind<-cbind(c(5,5,5,5,5,6,6,6,1,1,1,1,1,7),c(5,4,6,3,2,4,5,6,2,6,5,4,3,2))
Y.<-sapply(1:dim(y.ind)[1], function(i) rowMeans(Y[,atlas.ind[[y.ind[i,1]]][[y.ind[i,2]]]]))
#ynames<-sapply(1:dim(y.ind)[1], function(i) paste(names(atlas.ind)[y.ind[i,1]],names(atlas.ind[[y.ind[i,1]]])[y.ind[i,2]],sep=","))
ynames=c("par","occ","sm","fron","temp","mt","bg","thal","proj","cal","assoc","assoc.l","cer.c","cer.p")  


## outliers
W.<-outlier(Y.,sub=ind,within=TRUE,allY=TRUE,sd.thr=2,
	basename=paste(path.,"out",sep="/"))
## revised data/demographics
if(is.null(dim(Y.))) data.<-as.matrix(Y.[W.[[1]]]) else data.<-Y.[W.[[1]],]
save(data.,file=paste(path.,"data",sep="/"))
demo.<-demo[W.[[1]],]
save(demo.,file=paste(path.,"demo",sep="/"))

source(paste(path.scripts,"analysis.setup.R",sep="/"))
#> setup
#$path
#[1] "/home/danisimmonds/Dani/dti_0511/tbss/analysis/cross/roi_all1/ind"
#
#$mixed
#[1] FALSE
#
#$main
#[1] TRUE
#
#$numeric
#[1] TRUE
#
#$numeric.weights
#[1] "/home/danisimmonds/Dani/dti_0511/demographics/beh.out/vgs.mRT_cross.out_weights"     
#[2] "/home/danisimmonds/Dani/dti_0511/demographics/beh.out/vgs.cv_cross.out_weights"      
#[3] "/home/danisimmonds/Dani/dti_0511/demographics/beh.out/anti.percErr_cross.out_weights"
#[4] "/home/danisimmonds/Dani/dti_0511/demographics/beh.out/vgs.mRT_cross.out_weights"     
#
#$fixef
#$fixef$main
#     [,1]  [,2]   [,3]                  
#[1,] "age" "age"  "age"                 
#[2,] "lin" "poly" "ns"                  
#[3,] ""    "2"    "k=c(13,15.5,18,20.5)"
#
#$fixef$numeric
#     [,1]      [,2]     [,3]           [,4]     
#[1,] "vgs.mRT" "vgs.cv" "anti.percErr" "anti.cv"
#[2,] "lin"     "lin"    "lin"          "lin"    
#[3,] ""        ""       ""             ""       
#
#$fixef$categorical
#$fixef$categorical[[1]]
#     [,1] 
#[1,] "sex"
#[2,] "lin"
#[3,] ""   
#
#
#$ranef
#     [,1] 
#[1,] "age"
#[2,] "lin"
#[3,] ""   
#[4,] "id" 
#
#$deriv
#$deriv$nsim
#[1] 1000
#
#$deriv$main.range.int
#[1] 0.1
#
#$deriv$main.small.int
#[1] 1e-04
#
#$deriv$num.range.int
#[1] 40
#
#$deriv$vars
#$deriv$vars$main
#[1] 2
#
#$deriv$vars$numeric
#[1] 1 2 3 4
#
#$deriv$vars$categorical
#$deriv$vars$categorical[[1]]
#[1] 1
#
#$deriv$p.thr
#[1] 1

source(paste(path.scripts,"utilityfns","model.setup.R",sep="/"))
source(paste(path.scripts,"utilityfns","pred.str.R",sep="/"))
source(paste(path.scripts,"utilityfns","deriv.check.R",sep="/"))
#models<-model.setup(setup)

## more scripts for analysis.R --> need to update utilityfns
source(paste(path.scripts,"utilityfns","model.est.R",sep="/"))
source(paste(path.scripts,"utilityfns","fixef.set.R",sep="/"))
source(paste(path.scripts,"utilityfns","coef.est.R",sep="/"))
source(paste(path.scripts,"utilityfns","ll.test.R",sep="/"))
source(paste(path.scripts,"utilityfns","deriv.est.R",sep="/"))

