## required libraries
library(lme4) ## for mixed models
library(influence.ME) ## for regression diagnostics
library(splines) ## for splines
#library(Rniftilib) ## for working with nifti images

## load utility scripts
source(paste(setup$path.scripts,"utilityfns.R",sep="/"))

## load data, demographics
load(paste(setup$path,"X",sep="/"))
load(paste(setup$path,"Y",sep="/"))
if(is.null(dim(Y))) Y <- as.matrix(Y)

## models
models<-model.setup(setup)
sink(paste(setup$path,"models.txt",sep="/"))
cat(date(),"\tMODEL SETUP\n\n")
print(models)
sink()

## logfile
log.txt=paste(setup$path,"log.txt",sep="/")
sink(log.txt)
cat(date(),"\tLONGITUDINAL ANALYSIS\n\n")

## set up structures for estimation
models$model<-list()
cat(date(),"\t\t\t\testimating null model\n")
models$model[[1]]<-model.est() ## null model
cat(date(),"\t\t\t\testimation completed\n")

## model comparisons
#n.comp<-sum(sapply(1:length(models$con), function(i) length(models$con[[i]])))
n.comp<-dim(models$X)[1] ## CHANGE
labels<-character(n.comp)
llr<-matrix(NA,n.comp,dim(Y)[2])
llr.df<-numeric(n.comp)
llr.p<-matrix(NA,n.comp,dim(Y)[2])
ind.comp<-1

n<-length(unique(X$id)) ## number of subjects for cook's distance calculations

## executes custom code from analysis setup
if(length(setup$custom.code)>0) for(i in 1:length(setup$custom.code)) eval(parse(text=setup$custom.code[i]))

## loop through all models
for(m in 2:dim(models$X)[1]){

	## print progress through m
	cat(date(),"\t\tm =",m,"/",dim(models$X)[1],", vars =",models$X.var[m,],"\n")

	## create directory for each m
	path<-paste(setup$path,m,sep="/")
	dir.create(path)

	## estimate model
	cat(date(),"\t\t\t\testimating model\n")
	models$model[[m]]<-model.est(fixef.set(models$X[m,]))
	cat(date(),"\t\t\t\testimation completed\n")

	## loop through contrasts
	#for(c in 1:length(models$con[[m]])){
	#	## print progress through c
	#	cat(date(),"\t\t\tc =",c,"/",length(models$con[[m]]),"\n")
	#	## calculate LL test
	temp_model<-model.est2()
	ll.list<-ll.test(m0.=temp_model)
	labels[ind.comp]<-paste(models$model[[m]]$formula,temp_model$formula,sep="_vs_")
	llr[ind.comp,]<-ll.list[,1]
	llr.df[ind.comp]<-ll.list[1,2]
	llr.p[ind.comp,]<-ll.list[,3]
	ind.comp<-ind.comp+1
	#}

	# derivatives analysis
	if(models$deriv.do[m]>0){
		cat(date(),"\t\t\t\tderivatives analysis...\n")
		deriv.est(m0.=temp_model)
		cat(date(),"\t\t\t\tderivatives analysis completed\n")
	}
}

## create files for log-likelihood tests

#source(paste(setup$path.scripts,"utilityfns","ll.tbl.R",sep="/"))

## llr
filename=paste(setup$path,"llr",sep="/")
sink(filename)
cat("",setup$ynames,"\n")
sink()
write.table(llr,col.names=FALSE,row.names=labels,file=filename,append=TRUE)

## df
filename=paste(setup$path,"llr.df",sep="/")
write.table(llr.df,col.names=FALSE,row.names=labels,file=filename,append=TRUE)

## pvals
filename=paste(setup$path,"llr.p",sep="/")
sink(filename)
cat("",setup$ynames,"\n")
sink()
write.table(llr.p,col.names=FALSE,row.names=labels,file=filename,append=TRUE)

## significance stars
filename=paste(setup$path,"llr.p.star",sep="/")
sink(filename)
cat("",setup$ynames,"\n")
sink()
write.table(
	ifelse(llr.p<.001,"***",ifelse(llr.p<.01,"**",ifelse(llr.p<.05,"*",""))),
	col.names=FALSE,row.names=labels,file=filename,append=TRUE
)

## save models and close logfile
save(models,file=paste(setup$path,"models.done",sep="/"))
sink()
