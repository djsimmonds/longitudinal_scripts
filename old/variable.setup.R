## load demographics
#source("/home/danisimmonds/Dani/dti_0511/scripts/demographics.R")
#demo<-data.frame(
# nscans.ind=nscans.ind,
# nscans.tot=nscans.tot,
#	age=age,
#	id=id,
#	sex=sex,
#	tsr=tsr.2,
#	viq=viq,
#	piq=piq,
#	vgs.mRT=vgs.mRT,
#	vgs.sdRT=vgs.sdRT,
#	vgs.cv=vgs.cv,
#	vgs.mu=vgs.mu,
#	vgs.sigma=vgs.sigma,
#	vgs.tau=vgs.tau,
#	vgs.slow4=vgs.slow4,
#	anti.numErr=anti.numErr,
#	anti.percErr=anti.percErr,
#	anti.mRT=anti.mRT,
#	anti.sdRT=anti.sdRT,
#	anti.cv=anti.cv,
#	anti.mu=anti.mu,
#	anti.sigma=anti.sigma,
#	anti.tau=anti.tau,
#	anti.slow4corr=anti.slow4corr,
#	anti.slow4all=anti.slow4all
#)
#save(demo,file="/home/danisimmonds/Dani/dti_0511/demographics/demo.Rframe")
load("/home/danisimmonds/Dani/dti_0511/demographics/demo.Rframe")
## NEED TO RERUN TO PUT IN NSCANS.IND AND NSCANS.TOT

## code sample for getting data
#library(Rniftilib)
#paths<-list()
#paths$study<-"/home/danisimmonds/Dani/dti_0511"
#paths$data<-"/mnt/Schwarzenagger/Governator/DTI_STUDY/tbss/stats"
#mask<-nifti.image.read(paste(paths$data,"mean_FA_skeleton_mask_2mm_thr02_bin.nii.gz",sep="/"))
#d<-dim(mask)
#ind<-which(mask[,,]>0)
#i<-ind%%d[1]
#j<-(ind%%(d[1]*d[2]))%/%d[1]+1
#k<-ind%/%(d[1]*d[2])+1
#subjs<-nifti.image.read(paste(paths$data,"all_FA_skeletonised_2mm_div_mas.nii.gz",sep="/"))
#data.vox<-sapply(1:length(ind), function(n) subjs[i[n],j[n],k[n],])demo$age[ind]
#save(data.vox,file=paste(paths$study,"data.vox.Rframe",sep="/"))

## load data
load("/home/danisimmonds/Dani/dti_0511/tbss/data.vox.Rframe")
## atlas labels for ROI generation
load("/home/danisimmonds/Dani/dti_0511/tbss/atlas.ind.Rframe")

## outlier scripts
## new approach
	## no weighting
	## remove within-subject outliers at scan level (cleaning artifacts)
		## 1 - if >2 scans, fit linear model of age to each voxel, calculate residuals
		## 2 - average residuals across voxels for each each scan, calculate (res.avg-res.median)/(res.sd)
		## 3 - it appears that 0.5 or 1 are appropriate thresholds (2 scans @ 1, 3 scans @ 0.5)
	## will use regression analytics to remove outliers during analysis
library(lme4)
uid<-unique(demo$id)
out.ind<-array(NA,dim(data.vox)[1])
for(u in 1:length(uid)){
	ind<-which(demo$id==uid[u])
	while(length(ind)>2){
		cat(as.character(uid[u]),"\n")
		y<-as.numeric(data.vox[ind,])
		x1<-rep(demo$age[ind],dim(data.vox)[2])
		x2<-rep(1:dim(data.vox)[2],each=length(ind))
		m<-lmer(y~x1+(x1|x2))
		r<-m@resid
		r.med<-median(r)
		r.sd<-sd(r)
		r.avg<-rowMeans(matrix(r,length(ind),dim(data.vox)[2]))
		r.avg.n<-(r.avg-r.med)/r.sd
		i<-which.max(abs(r.avg))
		out.ind[ind[i]]<-r.avg.n[i]
		ind<-ind[-i]
	}
}
## plotting within-subject outliers to determine suitable cutoff
out<-which(abs(out.ind)>=0.5)
y<-rowMeans(data.vox)
plot(demo$age,y,pch=20,cex=0.5,xlab="age(y)",ylab="FA")
for(u in 1:length(uid)){
	ind<-which(demo$id==uid[u])
	lines(demo$age[ind],y[ind],lty=2)
}
for(i in 1:length(out.ind)) if(!is.na(out.ind[i])) if(abs(out.ind[i])>=1) points(demo$age[i],y[i],col="red") else if(abs(out.ind[i])>=0.5) points(demo$age[i],y[i],col="blue")
y.<-y[-out]
age.<-demo$age[-out]
id.<-demo$id[-out]
for(u in 1:length(uid)){
	ind<-which(id.==uid[u])
	lines(age.[ind],y.[ind])
}
## 3 removed: 30,43,74
## note - also used dixon outlier test on data averaged across all voxels
	## 10138, scan 3: p=0.008
	## 10152, scan 2: p=0.57
	## 10169, scan 3: p=0.056

## outliers removed
data.vox<-data.vox[-out,]
demo<-demo[-out,]
X<-demo
save(X,file="/home/danisimmonds/Dani/dti_0511/tbss/analysis/062011/X")
## L/R analyses
#X<-rbind(X,X)
#save(X,file="/home/danisimmonds/Dani/dti_0511/tbss/analysis/roiLR/X")

## normalized by mean of each voxel
vox.mean<-sapply(1:dim(data.vox)[2],function(v) mean(data.vox[,v]))
data.vox.n<-t(apply(data.vox,1,'/',vox.mean))

## maturation index (normalized by adult mean of each voxel)
adult.ind<-which(demo$age>18)
adult.mean<-sapply(1:dim(data.vox)[2],function(v) mean(data.vox[adult.ind,v]))
data.vox.mi<-t(apply(data.vox,1,'/',adult.mean))

## making various ROIs for analysis (will repeat all L/R)
ynames<-c("all","proj","cal","assoc","assoc.l","cer.c","cer.p","par","occ","sm","fron","temp","mt","bg","thal","ML","PCT","CST","CP","IC.A","IC.P","IC.R","CR.A","CR.S","CR.P","PTR","SFOF","SS","EC","UF","SLF","FOR.CB","FOR.C","CIN.CG","CIN.H","CAL.G","CAL.B","CAL.S","CAL.T","CER.I","CER.M","CER.S")
roi.ind<-list(
	1:dim(data.vox)[1],
	atlas.ind[[1]][[2]],
	atlas.ind[[1]][[5]],
	atlas.ind[[1]][[4]],
	atlas.ind[[1]][[6]],
	atlas.ind[[1]][[3]],
	atlas.ind[[7]][[2]],
	atlas.ind[[5]][[5]],
	atlas.ind[[5]][[4]],
	atlas.ind[[5]][[6]],
	atlas.ind[[5]][[3]],
	atlas.ind[[5]][[1]],
	atlas.ind[[6]][[4]],
	atlas.ind[[6]][[5]],
	atlas.ind[[6]][[6]],
	atlas.ind[[2]][[6]],
	atlas.ind[[2]][[4]],
	atlas.ind[[2]][[2]],
	atlas.ind[[2]][[9]],
	intersect(atlas.ind[[2]][[16]],atlas.ind[[3]][[5]]),
	intersect(atlas.ind[[2]][[16]],atlas.ind[[3]][[6]]),
	intersect(atlas.ind[[2]][[16]],atlas.ind[[3]][[7]]),
	intersect(atlas.ind[[2]][[14]],atlas.ind[[3]][[5]]),
	intersect(atlas.ind[[2]][[14]],atlas.ind[[3]][[10]]),
	intersect(atlas.ind[[2]][[14]],atlas.ind[[3]][[6]]),
	atlas.ind[[2]][[17]],
	atlas.ind[[2]][[22]],
	atlas.ind[[2]][[11]],
	atlas.ind[[2]][[13]],
	atlas.ind[[2]][[10]],
	atlas.ind[[2]][[18]],
	intersect(atlas.ind[[2]][[12]],atlas.ind[[3]][[8]]),
	intersect(atlas.ind[[2]][[12]],atlas.ind[[3]][[4]]),
	intersect(atlas.ind[[2]][[8]],atlas.ind[[3]][[9]]),
	intersect(atlas.ind[[2]][[8]],atlas.ind[[3]][[3]]),
	atlas.ind[[2]][[15]],
	atlas.ind[[2]][[21]],
	atlas.ind[[2]][[19]],
	atlas.ind[[2]][[20]],
	atlas.ind[[2]][[7]],
	atlas.ind[[2]][[3]],
	atlas.ind[[2]][[5]]
)

## ROI data
Y<-sapply(1:length(roi.ind),function(r) rowMeans(data.vox[,roi.ind[[r]]]))
## ROI data separated for L/R
#load("/home/danisimmonds/Dani/dti_0511/tbss/vox.Rframe")
#vox.l<-which(vox$x>0) ## radiological coords
#vox.r<-which(vox$x<0)
#Y<-sapply(1:length(roi.ind), function(r) rbind(rowMeans(data.vox[,intersect(roi.ind[[r]],vox.l)]),rowMeans(data.vox[,intersect(roi.ind[[r]],vox.r)])))
#save(Y,file="/home/danisimmonds/Dani/dti_0511/tbss/analysis/roiLR/Y")
save(Y,file="/home/danisimmonds/Dani/dti_0511/tbss/analysis/062011/Y")
Y.n<-sapply(1:length(roi.ind),function(r) rowMeans(data.vox.n[,roi.ind[[r]]]))
Y.mi<-sapply(1:length(roi.ind),function(r) rowMeans(data.vox.mi[,roi.ind[[r]]]))


## comparing curve fits for each data type
mean.Y<-apply(Y,2,mean)
mean.Y.n<-apply(Y.n,2,mean)
mean.Y.mi<-apply(Y.mi,2,mean)
adj.Y<-t(apply(Y,1,'/',mean.Y))
adj.Y.n<-t(apply(Y.n,1,'/',mean.Y.n))
adj.Y.mi<-t(apply(Y.mi,1,'/',mean.Y.mi))
age=demo$age
age.pred<-seq(8.5,28,0.5)
par(mfrow=c(2,2))
figs<-c(1,11,18,31)
for(i in 1:length(figs)){
	plot(range(age),range(c(adj.Y[,i],adj.Y.n[,i],adj.Y.mi[,i])),main=ynames[figs[i]],pch=30,xlab="age",ylab="FA")
	points(age,adj.Y[,i],pch=20,col=rgb(1,0,0,alpha=0.5))
	lines(age.pred,predict(loess(adj.Y[,i]~age),data.frame(age=age.pred)),col=rgb(1,0,0))
	points(age,adj.Y.n[,i],pch=20,col=rgb(0,1,0,alpha=0.5))
	lines(age.pred,predict(loess(adj.Y.n[,i]~age),data.frame(age=age.pred)),col=rgb(0,1,0),lty=2)
	points(age,adj.Y.mi[,i],pch=20,col=rgb(0,0,1,alpha=0.5))
	lines(age.pred,predict(loess(adj.Y.mi[,i]~age),data.frame(age=age.pred)),col=rgb(0,0,1),lty=3)
}

## model example for cook's distance calculations
library(influence.ME)
library(splines)
age.m<-round(mean(demo$age),2)
d<-data.frame(
	age=demo$age-age.m,
	id=demo$id,
	sex=demo$sex,
	y=Y[,1]
)
knots<-c(13-age.m,15.5-age.m,18-age.m)
m1<-lmer(y~ns(age,k=knots)+(age|id),d)
m1.est<-estex(m1,"id")
m1.cook<-ME.cook(m1.est)
m1.out<-which(m1.cook>(4/length(m1.cook)))
m2<-lmer(y~ns(age,k=knots)*sex+(age|id),d)
m2.est<-estex(m2,"id")
m2.cook<-ME.cook(m2.est)
m2.out<-which(m2.cook>(4/length(m2.cook)))
## plots
par(mfrow=c(1,2))
plot(range(age),range(y),main="age",pch=30,xlab="age",ylab="FA")
for(u in 1:length(uid)){
	ind<-which(id==uid[u])
	if(u%in%m1.out==TRUE) color<-"red" else color="black"
	if(length(ind)>1) lines(age[ind],Y[ind,1],col=color) else points(age[ind],Y[ind,1],pch=20,col=color)
}
plot(range(age),range(y),main="age*sex",pch=30,xlab="age",ylab="FA")
for(u in 1:length(uid)){
	ind<-which(id==uid[u])
	if(sex[min(ind)]=="m") { p.type=20;p.size=1;l.type=1 } else { p.type=1;p.size=0.5;l.type=2 }
	if(u%in%m2.out==TRUE) color<-"red" else color="black"
	if(length(ind)>1) lines(age[ind],Y[ind,1],lty=l.type,col=color) else points(age[ind],Y[ind,1],pch=p.type,cex=p.size,col=color)
}
