start.m<-232
sink(log.txt)

## loop through all models
for(m in start.m:dim(models$X)[1]){

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
	for(c in 1:length(models$con[[m]])){
		## print progress through c
		cat(date(),"\t\t\tc =",c,"/",length(models$con[[m]]),"\n")
		## calculate LL test
		ll.list<-ll.test()
		labels[ind.comp]<-paste(models$model[[m]]$formula,models$model[[models$con[[m]][c]]]$formula,sep="_vs_")
		llr[ind.comp,]<-ll.list[,1]
		llr.df[ind.comp]<-ll.list[1,2]
		llr.p[ind.comp,]<-ll.list[,3]
		ind.comp<-ind.comp+1
	}

	## estimate coefficients
	#if(models$coef.do[m]>0){
	#	cat(date(),"\t\t\t\tcoefficients analysis...\n")
	#	coef.est()
	#	cat(date(),"\t\t\t\tcoefficients analysis completed\n")
	#}

	# derivatives analysis
	if(models$deriv.do[m]>0){
		cat(date(),"\t\t\t\tderivatives analysis...\n")
		deriv.est()
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
