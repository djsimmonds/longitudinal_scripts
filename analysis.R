## required libraries
library(lme4) ## for mixed models
library(lmtest) ## for log-likelihood tests in regular mixed models (lm)
library(splines) ## for splines
library(AICcmodavg) ## for predicting from mixed models

data<-get(setup$data)
demo.<-get(setup$demo)
demo.$age<-demo.$age-8 ## converting to avoid singularity for correlation of random intercepts/slopes

## models
models<-model.setup()
sink(paste(setup$path,"models.txt",sep="/"))
cat(date(),"\tMODEL SETUP\n\n")
print(models)
sink()

## logfile
log.txt=paste(setup$path,"log.txt",sep="/")
sink(log.txt)
cat(date(),"\tLONGITUDINAL ANALYSIS\n\n")

## subset of time points
#cat(date(),"\tcalculating subset...\n")
#sub.ind<-subset.index()
#demo.<-demo.[sub.ind,]
#data<-data[sub.ind,]
#insert age.info function here
#sink(paste(setup$path,"subset.txt",sep="/"))
#cat(date(),"\tSUBSET\n\nn =",length(sub.ind),"\nind =",sub.ind)
#sink()
#cat(date(),"\tcalculation complete\n\n")

## convert continuous variables to ordered factors
#cat(date(),"\tconverting continuous variables to ordered factors...\n")
#sink(paste(setup$path,"cont.to.fact.txt",sep="/"))
#cat(date(),"\tCONVERSION OF CONTINUOUS VARIABLES TO ORDERED FACTORS\n\n")
#demo.<-cont.to.fact()
#sink()
#cat(date(),"\tconversion complete\n\n")

## calculate outliers
#cat(date(),"\tcalculating outliers...\n")
#out<-outlier()
#sink(paste(setup$path,"outliers.txt",sep="/"))
#cat(date(),"\tOUTLIERS\n\n")
#print(out)
#sink()
#cat(date(),"\tcalculation complete\n\n")

## load outlier weights
load(paste(setup$paths,"weights.Rdata",sep="/"))

## set up results files
label.txt=paste(setup$path,"label.txt",sep="/")
label.W.txt=paste(setup$path,"label.W.txt",sep="/")
anova.txt=paste(setup$path,"anova.txt",sep="/")
anova.W.txt=paste(setup$path,"anova.W.txt",sep="/")

## loop through each y
for(y in 1:dim(get(setup$data))[2]){

	## create directory for each y
	path.<-paste(setup$path,y,sep="/")
	dir.create(path.)

	## data frame for each y
	data.<-data.frame(data[,y],demo.,w=out$w[,y])
	names(data.)[1]<-setup$resp

	## print progress through y
	cat(date(),"\ty =",y,"/",dim(get(setup$data))[2],"\n")

	## set up structures for estimation and results
	models$model<-list()
	models$ll.test<-list()
	models$ll.test.W<-list()
	models$deriv<-list()

	## null model
	if(setup$mixed==TRUE){
		models$model[[1]]<-lmer.est()
	}else{
		models$model[[1]]<-lm.est()
	}

	## loop through all models
	for(m in 2:dim(models$X)[1]){

		## print progress through m
		cat(date(),"\t\tm =",m,"/",dim(models$X)[1],", vars =",models$X.var[m,],"\n")

		## create directory for each m
		path..<-paste(path.,m,sep="/")
		dir.create(path..)

		## estimate model
		if(setup$mixed==TRUE){
			models$model[[m]]<-lmer.est(fixef.set(models$X[m,]))
		}else{
			models$model[[m]]<-lm.est(fixef.set(models$X[m,]))
		}

		## log-likelihood tests for each contrast and derivatives analysis
		models$ll.test[[m]]<-list()
		for(c in 1:length(models$con[[m]])){

			## print progress through c
			cat(date(),"\t\t\tc =",c,"/",length(models$con[[m]]),"\n")

			## record contrast labels in label.txt for lining up with anova.txt
			if(y==1){
				sink(label.txt,append=TRUE)
				cat(models$model[[m]]$formula,"vs",models$model[[models$con[[m]][c]]]$formula,"\n")
				sink()
			}

			## calculate LL test and write p-value to "anova.txt" file
			comb.na<-numeric(0)
			if(!length(models$model[[m]]$na)==0 | !length(models$model[[models$con[[m]][c]]]$na)==0) comb.na<-unique(union(models$model[[m]]$na,models$model[[models$con[[m]][c]]]$na))
			models$ll.test[[m]][[c]]<-ll.test()

			## derivatives analysis
			if(c==1 & deriv.check(models$X[m,],models$ll.test[[m]][[c]]$p)){
				cat(date(),"\t\t\t\tderivatives analysis...\n")
				if(length(comb.na)==0){
					deriv.est()
				}else{
					if(setup$mixed==TRUE){
						m1<-lmer.est(models$model[[m]]$fixef,data=data.[-comb.na,])
						m0<-lmer.est(models$model[[models$con[[m]][c]]]$fixef,data=data.[-comb.na,])
					}else{
						m1<-lm.est(models$model[[m]]$fixef,data=data.[-comb.na,])
						m0<-lm.est(models$model[[models$con[[m]][c]]]$fixef,data=data.[-comb.na,])
					}
					deriv.est(m1,m0,data=data.[-comb.na,])
				}
				cat(date(),"\t\t\t\tderivatives analysis completed\n")
			}
		}

		## if within is true, calculate these LL tests also (in separate file)
		if(length(models$con.W[[m]])>0){
			models$ll.test.W[[m]]<-list()
			for(c in 1:length(models$con.W[[m]])){

				## print progress through c
				cat(date(),"\t\t\tc.W =",c,"/",length(models$con.W[[m]]),"\n")

				## record contrast labels in label.W.txt for lining up with anova.W.txt
				if(y==1){
					sink(label.W.txt,append=TRUE)
					cat(models$model[[m]]$formula,"vs",models$model[[models$con.W[[m]][c]]]$formula,"\n")
					sink()
				}

				## calculate LL test and write p-value to "anova.W.txt" file
				models$ll.test.W[[m]][[c]]<-ll.test(m0=models$model[[models$con.W[[m]][c]]],anova.=anova.W.txt)
			}
		}
	}

	## adjust results file for next iteration and save model
	sink(anova.txt,append=TRUE)
	cat("\n")
	sink()
	sink(anova.W.txt,append=TRUE)
	cat("\n")
	sink()
	cat("\nSaving model file...\n\n")
	save(models,file=paste(path.,"models.R",sep="/"))
}
sink()
