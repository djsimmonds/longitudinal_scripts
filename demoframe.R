#####################################
# longitudinal DTI study            #
#   ~	~	~	~	~	~	~	~                 #
# written by:       when:           #
#   Dani Simmonds     12/10 - 5/11  #
#####################################

## paths
paths<-list()
paths$study<-"/home/danisimmonds/Dani/dti_0511"

## run demographics script
source(paste(paths$study,"scripts","demographics.R",sep="/"))

## pick variables for inclusion 
demo<-data.frame(
	include=include,
	nscans.ind=nscans.ind,
	nscans.tot=nscans.tot,
	id=id,
	age=age,
	sex=sex,
	tsr=tsr.3,
	viq=viq.3,
	piq=piq.3,
	vgs.mRT=vgs.mRT,
	vgs.cv=vgs.cv,
	anti.percErr=anti.percErr,
	anti.cv=anti.cv
)

save(demo,file=paste(paths$study,"demo.Rframe",sep="/"))
