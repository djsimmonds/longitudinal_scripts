## paths
paths<-list()
paths$study<-"/home/danisimmonds/Dani/dti_0511"

## load data and demographics
load(paste(paths$study,"data.tract.Rframe",sep="/"))
load(paste(paths$study,"demo.Rframe",sep="/"))

tracts<-c(
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

## create plots directory if one does not already exists
if(length(which(dir(paths$study)=="plots"))==0) dir.create(paste(paths$study,"plots",sep="/"))
paths$plots<-paste(paths$study,"plots",sep="/")
if(length(which(dir(paths$plots)=="raw"))==0) dir.create(paste(paths$plots,"raw",sep="/"))
paths$plots<-paste(paths$plots,"raw",sep="/")
if(length(which(dir(paths$plots)=="tract"))==0) dir.create(paste(paths$plots,"tract",sep="/"))
paths$plots<-paste(paths$plots,"tract",sep="/")

## 1) all scans
## plots mean FA vs age for each tract, with points for subjects with single scans and lines connecting scans from the same subject
dir.create(paste(paths$plots,"all.scans",sep="/"))
for(i in 1:dim(data.tract)[2]){
	if(i<10) filename<-paste(0,i,".pdf",sep="") else filename<-paste(i,".pdf",sep="")
	pdf(paste(paths$plots,"all.scans",filename,sep="/"))
	plot(0,0,main=tracts[i],xlab="age(y)",xlim=c(min(demo$age),max(demo$age)),ylab="FA",ylim=c(min(data.tract[,i]),max(data.tract[,i])))
	for(j in 1:length(unique(demo$id))){
		id<-which(demo$id==unique(demo$id)[j])
		if(length(id)>1) lines(demo$age[id],data.tract[id,i]) else points(demo$age[id],data.tract[id,i],pch=20,cex=0.3)
	}
	dev.off()
}

## 2) excluding subjects with single scans
## only plotting first scan (or change/average of first two scans)
## series of four plots for each tract, all have dotted line for linear fit
	## top left --> mean FA vs age
	## top right --> average of mean FA at times 1 and 2 vs average of ages at time 1 and 2
	## bottom left --> slope of FA (change in FA/year) from time 1 to time 2 vs average of ages at time 1 and 2
	## bottom right --> %slope of FA (change in FA/year divided by FA at time 1 multiplied by 100 vs average of ages at time 1 and 2
dir.create(paste(paths$plots,"n2345.scans",sep="/"))
ind<-intersect(which(demo$nscans.tot>1),which(demo$nscans.ind==1))
age.m<-numeric(length(demo$age))
fa.m<-numeric(length(demo$age))
for(j in 1:length(ind)) age.m[ind[j]]<-mean(c(demo$age[ind[j]],demo$age[ind[j]+1]))
for(i in 1:dim(data.tract)[2]){
	for(j in 1:length(ind)) fa.m[ind[j]]<-mean(c(data.tract[ind[j],i],data.tract[ind[j]+1,i]))
	if(i<10) filename<-paste(0,i,".pdf",sep="") else filename<-paste(i,".pdf",sep="")
	pdf(paste(paths$plots,"n2345.scans",filename,sep="/"))
	par(mfrow=c(2,2))
	## top left
	plot(demo$age[ind],data.tract[ind,i],pch=20,cex=0.5,main=tracts[i],xlab="age at t1",ylab="FA at t1",ylim=c(min(data.tract[ind,i]),max(data.tract[ind,i])))
	abline(lm(data.tract[ind,i]~demo$age[ind]),lty=2)
	## top right
	plot(age.m[ind],fa.m[ind],pch=20,cex=0.5,xlab="mean(age at t1, age at t2)",ylab="mean(FA at t1, FA at t2)",ylim=c(min(data.tract[ind,i]),max(data.tract[ind,i])))
	abline(lm(fa.m[ind]~age.m[ind]),lty=2)
	## bottom left
	plot(age.m[ind],diff(data.tract[,i])[ind]/diff(demo$age)[ind],pch=20,cex=0.5,xlab="mean(age at t1, age at t2)",ylab="slope of FA from t1 to t2")
	abline(lm(diff(data.tract[,i])[ind]~age.m[ind]),lty=2)
	lines(c(min(age.m[ind]),max(age.m[ind])),c(0,0))
	## bottom right
	perc=100*diff(data.tract[,i])[ind]/data.tract[ind,i]/diff(demo$age)[ind]
	plot(age.m[ind],perc,pch=20,cex=0.5,xlab="mean(age at t1, age at t2)",ylab="slope of % change in FA from t1 to t2")
	abline(lm(perc~age.m[ind]),lty=2)
	lines(c(min(age.m[ind]),max(age.m[ind])),c(0,0))
	dev.off()
}

## 3) excluding subjects with <3 scans
## demonstrating outlier removal at a series of thresholds (2, 2.5, 3, 3.5, and 4 sd)
	## aaa: within subject residuals of fa~age, within tract, scan level
	## aab: within subject residuals of fa~age, within tract, subject level
	## aba: within subject residuals of fa~age, averaged across tracts, scan level
	## abb: within subject residuals of fa~age, averaged across tracts, subject level
	## baa: whole sample residuals of fa~age, within tract, scan level
	## bab: whole sample residuals of fa~age, within tract, subject level
	## bba: whole sample residuals of fa~age, averaged across tracts, scan level
	## bbb: whole sample residuals of fa~age, averaged across tracts, subject level
## series of four plots of fa~age for each tract
	## top left: original data
	## top right: outliers of within subject residuals marked (red=within tract, blue=across tract average)
	## bottom left: outliers of whole sample residuals marked (orange=within tract, green=across tract average)
	## bottom right: outliers removed
ind<-which(demo$nscans.tot>=3)
demo2=demo[ind,]
tract2=data.tract[ind,]
uniq<-unique(demo2$id)
thresh=c(2,2.5,3,3.5,4)
dir.create(paste(paths$plots,"n345.scans",sep="/"))
for(t in 1:length(thresh)){
	out<-list()
	## outliers based on within subject residuals of fa~age, applied at scan level
	out$aaa<-matrix(NA,dim(tract2)[1],dim(tract2)[2])
	out$aaa.thr<-numeric(dim(tract2)[2])
	## outliers based on within subject residuals of fa~age, applied at subject level
	out$aab<-matrix(NA,length(uniq),dim(tract2)[2])
	out$aab.thr<-numeric(dim(tract2)[2])
	## outliers based on whole sample residuals of fa~age, applied at scan level
	out$baa<-matrix(NA,dim(tract2)[1],dim(tract2)[2])
	out$baa.thr<-numeric(dim(tract2)[2])
	## outliers based on whole residuals of fa~age, applied at subject level
	out$bab<-matrix(NA,length(uniq),dim(tract2)[2])
	out$bab.thr<-numeric(dim(tract2)[2])
	## all 2sd outliers based on outlier measurements described above
	outliers<-array(TRUE,dim(tract2))
	outliers.subj<-array(TRUE,c(length(uniq),dim(tract2)[2]))
	## residuals within tracts
	for(i in 1:dim(tract2)[2]){
		out$baa[,i]<-resid(lm(tract2[,i]~demo2$age))^2
		for(j in 1:length(uniq)){
			ind<-which(demo2$id==uniq[j])
			out$aaa[ind,i]<-resid(lm(tract2[ind,i]~demo2$age[ind]))^2
			out$aab[j,i]<-mean(out$aaa[ind,i])
			out$bab[j,i]<-mean(out$baa[ind,i])		
		}
	}
	## residuals across tracts, and outlier threshold
	out$aba<-rowMeans(out$aaa)
	out$aba.thr<-mean(out$aba)+thresh[t]*sd(out$aba)
	out$abb<-rowMeans(out$aab)
	out$abb.thr<-mean(out$abb)+thresh[t]*sd(out$abb)
	out$bba<-rowMeans(out$baa)
	out$bba.thr<-mean(out$bba)+thresh[t]*sd(out$bba)
	out$bbb<-rowMeans(out$bab)
	out$bbb.thr<-mean(out$bbb)+thresh[t]*sd(out$bbb)
	## identify outliers across tracts
	for(j in 1:dim(tract2)[1]){
		if(out$aba[j]>out$aba.thr) outliers[j,]<-FALSE
		if(out$bba[j]>out$bba.thr) outliers[j,]<-FALSE
	}
	for(j in 1:length(uniq)){
		ind<-which(demo2$id==uniq[j])
		if(out$abb[j]>out$abb.thr){
			outliers[ind,]<-FALSE
			outliers.subj[j,]<-FALSE
		}
		if(out$bbb[j]>out$bbb.thr){
			outliers[ind,]<-FALSE
			outliers.subj[j,]<-FALSE
		}
	}
	## identify outliers within tracts, excluding those already identified across tracts
	for(i in 1:dim(tract2)[2]){
		## calculates 2sd threshold
		out$aaa.thr[i]<-mean(out$aaa[outliers[,i],i])+thresh[t]*sd(out$aaa[outliers[,i],i])
		out$aab.thr[i]<-mean(out$aab[outliers.subj[,i],i])+thresh[t]*sd(out$aab[outliers.subj[,i],i])
		out$baa.thr[i]<-mean(out$baa[outliers[,i],i])+thresh[t]*sd(out$baa[outliers[,i],i])
		out$bab.thr[i]<-mean(out$bab[outliers.subj[,i],i])+thresh[t]*sd(out$bab[outliers.subj[,i],i])
		## identify outliers
		for(j in 1:dim(tract2)[1]){
			if(out$aaa[j,i]>out$aaa.thr[i]) outliers[j,i]<-FALSE
			if(out$baa[j,i]>out$baa.thr[i]) outliers[j,i]<-FALSE
		}
		for(j in 1:length(uniq)){
			ind<-which(demo2$id==uniq[j])
			if(out$aab[j,i]>out$aab.thr[i]) outliers[ind,i]<-FALSE
			if(out$bab[j,i]>out$bab.thr[i]) outliers[ind,i]<-FALSE
		}
	}

	## plots
	dir.create(paste(paths$plots,"n345.scans",thresh[t],sep="/"))
	ind<-intersect(which(demo$nscans.tot>3),which(demo$nscans.ind==1))
	for(i in 1:dim(data.tract)[2]){
		if(i<10) filename<-paste(0,i,".pdf",sep="") else filename<-paste(i,".pdf",sep="")
		pdf(paste(paths$plots,"n345.scans",thresh[t],filename,sep="/"))
		par(mfrow=c(2,2),oma=c(1,1,3,1))
		## top left
		plot(0,0,main="original data",xlab="age(y)",xlim=c(min(demo2$age),max(demo2$age)),ylab="FA",ylim=c(min(tract2[,i]),max(tract2[,i])))
		for(j in 1:length(uniq)) lines(demo2$age[demo2$id==uniq[j]],tract2[demo2$id==uniq[j],i])
		## top right
		plot(0,0,main="within subject outliers",xlab="age(y)",xlim=c(min(demo2$age),max(demo2$age)),ylab="FA",ylim=c(min(tract2[,i]),max(tract2[,i])))
		for(j in 1:length(uniq)) lines(demo2$age[demo2$id==uniq[j]],tract2[demo2$id==uniq[j],i])
		for(j in 1:length(uniq)){
			if(out$aab[j,i]>out$aab.thr[i]) lines(demo2$age[demo2$id==uniq[j]],tract2[demo2$id==uniq[j],i],col="red",lwd=2)
			if(out$abb[j]>out$abb.thr) lines(demo2$age[demo2$id==uniq[j]],tract2[demo2$id==uniq[j],i],col="blue",lty=2,lwd=2)
		}
		for(j in 1:dim(tract2)[1]){
			if(out$aba[j]>out$aba.thr) points(demo2$age[j],tract2[j,i],pch=20,cex=2,col="blue")
			if(out$aaa[j,i]>out$aaa.thr[i]) points(demo2$age[j],tract2[j,i],pch=20,cex=1,col="red")
		}
		## bottom left
		plot(0,0,main="whole sample outliers",xlab="age(y)",xlim=c(min(demo2$age),max(demo2$age)),ylab="FA",ylim=c(min(tract2[,i]),max(tract2[,i])))
		for(j in 1:length(uniq)) lines(demo2$age[demo2$id==uniq[j]],tract2[demo2$id==uniq[j],i])
		for(j in 1:length(uniq)){
			if(out$bab[j,i]>out$bab.thr[i]) lines(demo2$age[demo2$id==uniq[j]],tract2[demo2$id==uniq[j],i],col="orange",lwd=2)
			if(out$bbb[j]>out$bbb.thr) lines(demo2$age[demo2$id==uniq[j]],tract2[demo2$id==uniq[j],i],col="green",lty=2,lwd=2)
		}
		for(j in 1:dim(tract2)[1]){
			if(out$bba[j]>out$bba.thr) points(demo2$age[j],tract2[j,i],pch=20,cex=2,col="green")
			if(out$baa[j,i]>out$baa.thr[i]) points(demo2$age[j],tract2[j,i],pch=20,cex=1,col="orange")
		}
		## bottom right
		demo.out<-demo2[outliers[,i],]
		tract.out<-tract2[outliers[,i],i]	
		plot(0,0,main="outliers removed",xlab="age(y)",xlim=c(min(demo2$age),max(demo2$age)),ylab="FA",ylim=c(min(tract2[,i]),max(tract2[,i])))
		for(j in 1:length(unique(demo.out$id))){
			id<-which(demo.out$id==unique(demo.out$id)[j])
			if(length(id)>1) lines(demo.out$age[id],tract.out[id]) else points(demo.out$age[id],tract.out[id],pch=20,cex=0.3)
		}
		mtext(paste(tracts[i],", thresh=",thresh[t],"sd"),outer=TRUE,line=1)
		dev.off()
	}
}

## 4) excluding subjects with <4 scans
## outliers as described previously (threshold 2sd)
## series of two plots of fa~age and two plots of slope of fa~age for each tract 
	## top left: fa~age, outliers marked
	## top right: fa~age, outliers removed
	## bottom left: slope of fa~age, outliers marked
	## bottom right: slope of fa~age, outliers removed
ind<-which(demo$nscans.tot>=4)
demo2=demo[ind,]
tract2=data.tract[ind,]
uniq<-unique(demo2$id)
dir.create(paste(paths$plots,"n45.scans",sep="/"))
out<-list()
outliers<-list()
outliers.subj<-list()

out$int<-list()
## outliers based on within subject residuals of fa~age, applied at scan level
out$int$aaa<-matrix(NA,dim(tract2)[1],dim(tract2)[2])
out$int$aaa.thr<-numeric(dim(tract2)[2])
## outliers based on within subject residuals of fa~age, applied at subject level
out$int$aab<-matrix(NA,length(uniq),dim(tract2)[2])
out$int$aab.thr<-numeric(dim(tract2)[2])
## outliers based on whole sample residuals of fa~age, applied at scan level
out$int$baa<-matrix(NA,dim(tract2)[1],dim(tract2)[2])
out$int$baa.thr<-numeric(dim(tract2)[2])
## outliers based on whole residuals of fa~age, applied at subject level
out$int$bab<-matrix(NA,length(uniq),dim(tract2)[2])
out$int$bab.thr<-numeric(dim(tract2)[2])
## all 2sd outliers based on outlier measurements described above
outliers$int<-array(TRUE,dim(tract2))
outliers.subj$int<-array(TRUE,c(length(uniq),dim(tract2)[2]))
## residuals within tracts
for(i in 1:dim(tract2)[2]){
	out$int$baa[,i]<-resid(lm(tract2[,i]~demo2$age))^2
	for(j in 1:length(uniq)){
		ind<-which(demo2$id==uniq[j])
		out$int$aaa[ind,i]<-resid(lm(tract2[ind,i]~demo2$age[ind]))^2
		out$int$aab[j,i]<-mean(out$int$aaa[ind,i])
		out$int$bab[j,i]<-mean(out$int$baa[ind,i])		
	}
}
## residuals across tracts, and outlier threshold
out$int$aba<-rowMeans(out$int$aaa)
out$int$aba.thr<-mean(out$int$aba)+2*sd(out$int$aba)
out$int$abb<-rowMeans(out$int$aab)
out$int$abb.thr<-mean(out$int$abb)+2*sd(out$int$abb)
out$int$bba<-rowMeans(out$int$baa)
out$int$bba.thr<-mean(out$int$bba)+2*sd(out$int$bba)
out$int$bbb<-rowMeans(out$int$bab)
out$int$bbb.thr<-mean(out$int$bbb)+2*sd(out$int$bbb)
## identify outliers across tracts
for(j in 1:dim(tract2)[1]){
	if(out$int$aba[j]>out$int$aba.thr) outliers$int[j,]<-FALSE
	if(out$int$bba[j]>out$int$bba.thr) outliers$int[j,]<-FALSE
}
for(j in 1:length(uniq)){
	ind<-which(demo2$id==uniq[j])
	if(out$int$abb[j]>out$int$abb.thr){
		outliers$int[ind,]<-FALSE
		outliers.subj$int[j,]<-FALSE
	}
	if(out$int$bbb[j]>out$int$bbb.thr){
		outliers$int[ind,]<-FALSE
		outliers.subj$int[j,]<-FALSE
	}
}
## identify outliers within tracts, excluding those already identified across tracts
for(i in 1:dim(tract2)[2]){
	## calculates 2sd threshold
	out$int$aaa.thr[i]<-mean(out$int$aaa[outliers$int[,i],i])+2*sd(out$int$aaa[outliers$int[,i],i])
	out$int$aab.thr[i]<-mean(out$int$aab[outliers.subj$int[,i],i])+2*sd(out$int$aab[outliers.subj$int[,i],i])
	out$int$baa.thr[i]<-mean(out$int$baa[outliers$int[,i],i])+2*sd(out$int$baa[outliers$int[,i],i])
	out$int$bab.thr[i]<-mean(out$int$bab[outliers.subj$int[,i],i])+2*sd(out$int$bab[outliers.subj$int[,i],i])
	## identify outliers
	for(j in 1:dim(tract2)[1]){
		if(out$int$aaa[j,i]>out$int$aaa.thr[i]) outliers$int[j,i]<-FALSE
		if(out$int$baa[j,i]>out$int$baa.thr[i]) outliers$int[j,i]<-FALSE
	}
	for(j in 1:length(uniq)){
		ind<-which(demo2$id==uniq[j])
		if(out$int$aab[j,i]>out$int$aab.thr[i]) outliers$int[ind,i]<-FALSE
		if(out$int$bab[j,i]>out$int$bab.thr[i]) outliers$int[ind,i]<-FALSE
	}
}

## indexes which differences to count (only within subject)
ind.diff<-which(demo2$nscans.ind!=1)-1
age.m<-demo2$age[ind.diff]+diff(demo2$age)[ind.diff]/2
tract2.slope<-matrix(NA,length(ind.diff),dim(tract2)[2])
for(i in 1:dim(tract2)[2]) tract2.slope[,i]<-diff(tract2[,i])[ind.diff]/diff(demo2$age)[ind.diff]
demo.diff<-demo2[ind.diff,]
uniq.diff<-unique(demo.diff$id)

out$slope<-list()
## outliers based on within subject residuals of fa~age, applied at scan level
out$slope$aaa<-matrix(NA,dim(tract2.slope)[1],dim(tract2.slope)[2])
out$slope$aaa.thr<-numeric(dim(tract2.slope)[2])
## outliers based on within subject residuals of fa~age, applied at subject level
out$slope$aab<-matrix(NA,length(uniq.diff),dim(tract2.slope)[2])
out$slope$aab.thr<-numeric(dim(tract2.slope)[2])
## outliers based on whole sample residuals of fa~age, applied at scan level
out$slope$baa<-matrix(NA,dim(tract2.slope)[1],dim(tract2.slope)[2])
out$slope$baa.thr<-numeric(dim(tract2.slope)[2])
## outliers based on whole residuals of fa~age, applied at subject level
out$slope$bab<-matrix(NA,length(uniq.diff),dim(tract2.slope)[2])
out$slope$bab.thr<-numeric(dim(tract2.slope)[2])
## all 2sd outliers based on outlier measurements described above
outliers$slope<-array(TRUE,dim(tract2.slope))
outliers.subj$slope<-array(TRUE,c(length(uniq.diff),dim(tract2.slope)[2]))
## residuals within tracts
for(i in 1:dim(tract2.slope)[2]){
	out$slope$baa[,i]<-resid(lm(tract2.slope[,i]~age.m))^2
	for(j in 1:length(uniq.diff)){
		ind<-which(demo.diff$id==uniq.diff[j])
		out$slope$aaa[ind,i]<-resid(lm(tract2.slope[ind,i]~age.m[ind]))^2
		out$slope$aab[j,i]<-mean(out$slope$aaa[ind,i])
		out$slope$bab[j,i]<-mean(out$slope$baa[ind,i])		
	}
}
## residuals across tracts, and outlier threshold
out$slope$aba<-rowMeans(out$slope$aaa)
out$slope$aba.thr<-mean(out$slope$aba)+2*sd(out$slope$aba)
out$slope$abb<-rowMeans(out$slope$aab)
out$slope$abb.thr<-mean(out$slope$abb)+2*sd(out$slope$abb)
out$slope$bba<-rowMeans(out$slope$baa)
out$slope$bba.thr<-mean(out$slope$bba)+2*sd(out$slope$bba)
out$slope$bbb<-rowMeans(out$slope$bab)
out$slope$bbb.thr<-mean(out$slope$bbb)+2*sd(out$slope$bbb)
## identify outliers across tracts
for(j in 1:dim(tract2.slope)[1]){
	if(out$slope$aba[j]>out$slope$aba.thr) outliers$slope[j,]<-FALSE
	if(out$slope$bba[j]>out$slope$bba.thr) outliers$slope[j,]<-FALSE
}
for(j in 1:length(uniq.diff)){
	ind<-which(demo.diff$id==uniq.diff[j])
	if(out$slope$abb[j]>out$slope$abb.thr){
		outliers$slope[ind,]<-FALSE
		outliers.subj$slope[j,]<-FALSE
	}
	if(out$slope$bbb[j]>out$slope$bbb.thr){
		outliers$slope[ind,]<-FALSE
		outliers.subj$slope[j,]<-FALSE
	}
}
## identify outliers within tracts, excluding those already identified across tracts
for(i in 1:dim(tract2.slope)[2]){
	## calculates 2sd threshold
	out$slope$aaa.thr[i]<-mean(out$slope$aaa[outliers$slope[,i],i])+2*sd(out$slope$aaa[outliers$slope[,i],i])
	out$slope$aab.thr[i]<-mean(out$slope$aab[outliers.subj$slope[,i],i])+2*sd(out$slope$aab[outliers.subj$slope[,i],i])
	out$slope$baa.thr[i]<-mean(out$slope$baa[outliers$slope[,i],i])+2*sd(out$slope$baa[outliers$slope[,i],i])
	out$slope$bab.thr[i]<-mean(out$slope$bab[outliers.subj$slope[,i],i])+2*sd(out$slope$bab[outliers.subj$slope[,i],i])
	## identify outliers
	for(j in 1:dim(tract2.slope)[1]){
		if(out$slope$aaa[j,i]>out$slope$aaa.thr[i]) outliers$slope[j,i]<-FALSE
		if(out$slope$baa[j,i]>out$slope$baa.thr[i]) outliers$slope[j,i]<-FALSE
	}
	for(j in 1:length(uniq.diff)){
		ind<-which(demo.diff$id==uniq.diff[j])
		if(out$slope$aab[j,i]>out$slope$aab.thr[i]) outliers$slope[ind,i]<-FALSE
		if(out$slope$bab[j,i]>out$slope$bab.thr[i]) outliers$slope[ind,i]<-FALSE
	}
}

## plots
for(i in 1:dim(data.tract)[2]){
	if(i<10) filename<-paste(0,i,".pdf",sep="") else filename<-paste(i,".pdf",sep="")
	pdf(paste(paths$plots,"n45.scans",filename,sep="/"))
	par(mfrow=c(2,2),oma=c(1,1,3,1))
	## top left
	plot(0,0,main="fa~age, outliers marked",xlab="age(y)",xlim=c(min(demo2$age),max(demo2$age)),ylab="FA",ylim=c(min(tract2[,i]),max(tract2[,i])))
	for(j in 1:length(uniq)) lines(demo2$age[demo2$id==uniq[j]],tract2[demo2$id==uniq[j],i])
	abline(lm(tract2[,i]~demo2$age),lty=2)
	for(j in 1:length(uniq)){
		if(out$int$aab[j,i]>out$int$aab.thr[i]) lines(demo2$age[demo2$id==uniq[j]],tract2[demo2$id==uniq[j],i],col="red",lwd=2)
		if(out$int$abb[j]>out$int$abb.thr) lines(demo2$age[demo2$id==uniq[j]],tract2[demo2$id==uniq[j],i],col="blue",lty=2,lwd=2)
		if(out$int$bab[j,i]>out$int$bab.thr[i]) lines(demo2$age[demo2$id==uniq[j]],tract2[demo2$id==uniq[j],i],col="orange",lty=3,lwd=2)
		if(out$int$bbb[j]>out$int$bbb.thr) lines(demo2$age[demo2$id==uniq[j]],tract2[demo2$id==uniq[j],i],col="green",lty=4,lwd=2)
	}
	for(j in 1:dim(tract2)[1]){
		if(out$int$aba[j]>out$int$aba.thr) points(demo2$age[j],tract2[j,i],pch=20,cex=2.5,col="blue")
		if(out$int$bba[j]>out$int$bba.thr) points(demo2$age[j],tract2[j,i],pch=20,cex=2,col="green")
		if(out$int$aaa[j,i]>out$int$aaa.thr[i]) points(demo2$age[j],tract2[j,i],pch=20,cex=1.5,col="red")
		if(out$int$baa[j,i]>out$int$baa.thr[i]) points(demo2$age[j],tract2[j,i],pch=20,cex=1,col="orange")
	}
	## top right
	demo.out<-demo2[outliers$int[,i],]
	tract.out<-tract2[outliers$int[,i],i]	
	plot(0,0,main="fa~age, outliers removed",xlab="age(y)",xlim=c(min(demo2$age),max(demo2$age)),ylab="FA",ylim=c(min(tract2[,i]),max(tract2[,i])))
	for(j in 1:length(unique(demo.out$id))){
		id<-which(demo.out$id==unique(demo.out$id)[j])
		if(length(id)>1) lines(demo.out$age[id],tract.out[id]) else points(demo.out$age[id],tract.out[id],pch=20,cex=0.3)
	}
	abline(lm(tract.out~demo.out$age),lty=2)
	## bottom left
	plot(0,0,main="slope of fa~age, outliers marked",xlab="age(y)",xlim=c(min(age.m),max(age.m)),ylab="FA",ylim=c(min(tract2.slope[,i]),max(tract2.slope[,i])))
	lines(c(min(age.m)-1,max(age.m)+1),c(0,0),lty=3)
	for(j in 1:length(uniq.diff)) lines(age.m[demo2$id[ind.diff]==uniq.diff[j]],tract2.slope[demo2$id[ind.diff]==uniq.diff[j],i])
	abline(lm(tract2.slope[,i]~age.m),lty=2)
	for(j in 1:length(uniq.diff)){
		if(out$slope$aab[j,i]>out$slope$aab.thr[i]) lines(age.m[demo2$id[ind.diff]==uniq.diff[j]],tract2.slope[demo2$id[ind.diff]==uniq.diff[j],i],col="red",lwd=2)
		if(out$slope$abb[j]>out$slope$abb.thr) lines(age.m[demo2$id[ind.diff]==uniq.diff[j]],tract2.slope[demo2$id[ind.diff]==uniq.diff[j],i],col="blue",lty=2,lwd=2)
		if(out$slope$bab[j,i]>out$slope$bab.thr[i]) lines(age.m[demo2$id[ind.diff]==uniq.diff[j]],tract2.slope[demo2$id[ind.diff]==uniq.diff[j],i],col="orange",lty=3,lwd=2)
		if(out$slope$bbb[j]>out$slope$bbb.thr) lines(age.m[demo2$id[ind.diff]==uniq.diff[j]],tract2.slope[demo2$id[ind.diff]==uniq.diff[j],i],col="green",lty=4,lwd=2)
	}
	for(j in 1:dim(tract2.slope)[1]){
		if(out$slope$aba[j]>out$slope$aba.thr) points(age.m[j],tract2.slope[j,i],pch=20,cex=2.5,col="blue")
		if(out$slope$bba[j]>out$slope$bba.thr) points(age.m[j],tract2.slope[j,i],pch=20,cex=2,col="green")
		if(out$slope$aaa[j,i]>out$slope$aaa.thr[i]) points(age.m[j],tract2.slope[j,i],pch=20,cex=1.5,col="red")
		if(out$slope$baa[j,i]>out$slope$baa.thr[i]) points(age.m[j],tract2.slope[j,i],pch=20,cex=1,col="orange")
	}
	## bottom right
	demo.out<-demo.diff[outliers$slope[,i],]
	tract.out<-tract2.slope[outliers$slope[,i],i]
	age.m.out<-age.m[outliers$slope[,i]]
	plot(0,0,main="slope of fa~age, outliers removed",xlab="age(y)",xlim=c(min(age.m),max(age.m)),ylab="FA",ylim=c(min(tract2.slope[,i]),max(tract2.slope[,i])))
	lines(c(min(age.m)-1,max(age.m)+1),c(0,0),lty=3)
	for(j in 1:length(unique(demo.out$id))){
		id<-which(demo.out$id==unique(demo.out$id)[j])
		if(length(id)>1) lines(age.m.out[id],tract.out[id]) else points(age.m.out[id],tract.out[id],pch=20,cex=0.3)
	}
	abline(lm(tract.out~age.m.out),lty=2)
	mtext(tracts[i],outer=TRUE,line=1)
	dev.off()
}

