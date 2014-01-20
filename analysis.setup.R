## paths
paths<-list()
paths$study<-"/home/danisimmonds/Dani/dti_0511"
paths$scripts<-"/home/danisimmonds/Dropbox/scripts"
paths$analysis<-paste(paths$study,"tbss","analysis","cross","vox",sep="/")

## load utility functions
source(paste(paths$scripts,"utilityfns.R",sep="/"))

## load data, demographics and weights
load(paste(paths$analysis,"data.Rdata",sep="/"))
load(paste(paths$analysis,"demo.Rdata",sep="/"))
load(paste(paths$analysis,"weights.Rdata",sep="/"))

## setup
setup<-list(
	path=paste(paths$analysis,"age.only",sep="/"), ## make sure it exists!
	## data frame containing Y variables
	data="data.",
	## data frame containing X variables
	demo="demo",
	## convert continuous variables to ordered factors
	cont.to.fact=list(
		#vars=c("vgs.mRT","vgs.cv","anti.percErr","anti.cv"),
		#cut=list(c(1/3,2/3),c(1/3,2/3),c(1/3,2/3),c(1/3,2/3))
	),
	## response variable
	resp="fa",
	## mixed or regular regression
	mixed=TRUE,
	## subset of time points for analysis (variable in demo, logical operation, values)
	subset=list(
		#list("nscans.tot","g",3)
		#list("nscans.ind","e",1)
	),
	## outliers
		## note: Y, X and id are included here, but changing them won't do anything 
	#out=list(
	#	Y="",
	#	X="age", 
	#	id="id",
	#	within=TRUE,
	#	allY=TRUE,
	#	sd.thr=2
	#),
	## fixed effects (var,type,param)
	fixef=list(
		## main models
		cbind(
			c("age","lin",""),
			c("age","ns","2")
		)#,
		## interactions with main models
		#cbind(c("sex","lin","")),
		## interactions with main, int(1) and main*int(1) models
		#cbind(
		#	#c("viq","lin",""),
		#	#c("piq","lin",""),
		#	c("vgs.mRT","lin",""),
		#	c("vgs.cv","lin",""),
		#	c("anti.percErr","lin",""),
		#	c("anti.cv","lin","")
		#)
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
		vars=list(2), ## corresponds to models structure above, add more levels as necessary
		p.thr=1
	)
)
save(setup,file=paste(setup$path,"model.setup",sep="/"))
