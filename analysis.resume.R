resume.m<-12

## logfile
log.txt=paste(setup$path,"log.txt",sep="/")
sink(log.txt)
cat(date(),"\tLONGITUDINAL ANALYSIS\n\n")

## loop through all models
for(m in resume.m:dim(models$X)[1]){

	## print progress through m
	cat(date(),"\t\tm =",m,"/",dim(models$X)[1],", vars =",models$X.var[m,],"\n")

	## create directory for each m
	path<-paste(setup$path,m,sep="/")
	dir.create(path)

	## estimate model
	models$model[[m]]<-model.est(fixef.set(models$X[m,]))

	## extra weights
	if(setup$numeric==TRUE & !is.na(models$X[m,2])) wts.<-Wts.[[models$X[m,2]]] else wts.<-NULL

	## loop through contrasts
	for(c in 1:length(models$con[[m]])){

		## print progress through c
		cat(date(),"\t\t\tc =",c,"/",length(models$con[[m]]),"\n")

		## calculate LL test
		ll.list<-ll.test()
		labels[ind.comp]<-paste(models$model[[m]]$formula,models$model[[models$con[[m]][c]]]$formula,sep="_vs_")
		llr[ind.comp,]<-ll.list[1,]
		llr.df[ind.comp]<-ll.list[2,1]
		llr.p[ind.comp,]<-ll.list[3,]
		ind.comp<-ind.comp+1
	}

	## estimate coefficients
	coef.est(data.,Wts,wts.,m1.=models$model[[m]],m0.=models$model[[1]],path=paste(setup$path,m,sep="/"),mixed=setup$mixed,ynames=setup$ynames)

	## derivatives analysis
	if(models$deriv.do[m]>0){
		cat(date(),"\t\t\t\tderivatives analysis...\n")
		deriv.est()
		cat(date(),"\t\t\t\tderivatives analysis completed\n")
	}

}

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
filename=paste(path,"llr.p.star",sep="/")
sink(filename)
cat("",ynames,"\n")
sink()
write.table(ifelse(llr.p<.001,"***",ifelse(llr.p<.01,"**",ifelse(llr.p<.05,"*",""))),
	col.names=FALSE,row.names=labels,file=filename,append=TRUE)

save(models,file=paste(setup$path,"models.done",sep="/"))

sink()
