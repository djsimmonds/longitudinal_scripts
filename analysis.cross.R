## paths
paths<-list()
paths$study<-"/home/danisimmonds/Dani/dti_0511"

## load utility functions
source(paste(paths$study,"scripts","utilityfns.R",sep="/"))

## load data and demographics
load(paste(paths$study,"data.tract.Rframe",sep="/"))
load(paste(paths$study,"demo.Rframe",sep="/"))
#tracts<-c("Middle cerebellar peduncle","Pontine crossing tract (a part of MCP)","Genu of corpus callosum","Body of corpus callosum","Splenium of corpus callosum","Fornix (column and body of fornix)","Corticospinal tract R","Corticospinal tract L","Medial lemniscus R","Medial lemniscus L","Inferior cerebellar peduncle R","Inferior cerebellar peduncle L","Superior cerebellar peduncle R","Superior cerebellar peduncle L","Cerebral peduncle R","Cerebral peduncle L","Anterior limb of internal capsule R","Anterior limb of internal capsule L","Posterior limb of internal capsule R","Posterior limb of internal capsule L","Retrolenticular part of internal capsule R","Retrolenticular part of internal capsule L","Anterior corona radiata R","Anterior corona radiata L","Superior corona radiata R","Superior corona radiata L","Posterior corona radiata R","Posterior corona radiata L","Posterior thalamic radiation (include optic radiation) R","Posterior thalamic radiation (include optic radiation) L","Sagittal stratum (include inferior longitidinal fasciculus and inferior fronto-occipital fasciculus) R","Sagittal stratum (include inferior longitidinal fasciculus and inferior fronto-occipital fasciculus) L","External capsule R","External capsule L","Cingulum (cingulate gyrus) R","Cingulum (cingulate gyrus) L","Cingulum (hippocampus) R","Cingulum (hippocampus) L","Fornix (cres) / Stria terminalis (can not be resolved with current resolution) R","Fornix (cres) / Stria terminalis (can not be resolved with current resolution) L","Superior longitudinal fasciculus R","Superior longitudinal fasciculus L","Superior fronto-occipital fasciculus (could be a part of anterior internal capsule) R","Superior fronto-occipital fasciculus (could be a part of anterior internal capsule) L","Uncinate fasciculus R","Uncinate fasciculus L","Tapetum R","Tapetum L")

## setup
setup<-list(
	path="/home/danisimmonds/Dani/dti_0511/analysis/tract.cross", ## make sure it exists!
	## data frame containing Y variables
	data="data.tract",
	## data frame containing X variables
	demo="demo",
	## convert continuous variables to ordered factors
	cont.to.fact=list(
		vars=c("vgs.mRT","vgs.cv","anti.percErr","anti.cv"),
		cut=list(c(1/3,2/3),c(1/3,2/3),c(1/3,2/3),c(1/3,2/3))
	),
	## response variable
	resp="fa",
	## mixed or regular regression
	mixed=FALSE,
	## subset of time points for analysis (variable in demo, logical operation, values)
	subset=list(
		#list("nscans.tot","g",3)
		list("nscans.ind","e",1)
	),
	## outliers
		## note: Y, X and id are included here, but changing them won't do anything 
	out=list(
		Y="",
		X="age", 
		id="id",
		within=FALSE,
		allY=TRUE,
		sd.thr=2
	),
	## fixed effects (var,type,param)
	fixef=list(
		## main models
		cbind(
			c("age","lin",""),
			c("age","ns","df=3")
		),
		## interactions with main models
		cbind(c("sex","lin","")),
		## interactions with main, int(1) and main*int(1) models
		cbind(
			#c("viq","lin",""),
			#c("piq","lin",""),
			c("vgs.mRT","lin",""),
			c("vgs.cv","lin",""),
			c("anti.percErr","lin",""),
			c("anti.cv","lin","")
		)
		## more nested interactions if necessary go here
	),
	## random effects (var,type,param,grp), if longitudinal
	ranef=cbind(c("age","lin","","id")),
	## compare main model effects to each other?
	within=TRUE,
	## derivatives analysis (only with main variables and their interactions)
	deriv=list(
		nsim=1000,
		range=seq(8.5,24.5,0.1), ## values of main at which derivatives are examined
		int=1e-4, ## tiny interval for derivative calculations
		vars=list(2,1,3), ## corresponds to models structure above, add more levels as necessary
		p.thr=1
	)
)
save(setup,file=paste(setup$path,"model.setup",sep="/"))

##############################################################################################

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
cat(date(),"\tcalculating subset...\n")
sub.ind<-subset.index()
demo.<-demo.[sub.ind,]
data<-data[sub.ind,]
#insert age.info function here
sink(paste(setup$path,"subset.txt",sep="/"))
cat(date(),"\tSUBSET\n\nn =",length(sub.ind),"\nind =",sub.ind)
sink()
cat(date(),"\tcalculation complete\n\n")

## convert continuous variables to ordered factors
cat(date(),"\tconverting continuous variables to ordered factors...\n")
sink(paste(setup$path,"cont.to.fact.txt",sep="/"))
cat(date(),"\tCONVERSION OF CONTINUOUS VARIABLES TO ORDERED FACTORS\n\n")
demo.<-cont.to.fact()
sink()
cat(date(),"\tconversion complete\n\n")

## calculate outliers
cat(date(),"\tcalculating outliers...\n")
out<-outlier()
sink(paste(setup$path,"outliers.txt",sep="/"))
cat(date(),"\tOUTLIERS\n\n")
print(out)
sink()
cat(date(),"\tcalculation complete\n\n")

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
