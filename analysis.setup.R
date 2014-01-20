## setup
setup<-list(
	path="/Volumes/Governator/ANTISTATELONG/ROIs/Data/dani_scripts/analysis20120508",
	path.scripts="/Volumes/Governator/ANTISTATELONG/ROIs/Data/dani_scripts/scripts0512",
	ynames=c("dlPFC_L", "dlPFC_R", "vlPFC_L", "vlPFC_R", "insula_L", "insula_R", "SEF", "preSMA", "FEF_L", "FEF_R", "putamen_L", "putamen_R", "PPC_L", "PPC_R", "V1_bilat",	"cerebellum_L", "cerebellum_R", "dACC"),
	mixed=TRUE, ## longitudinal?
	main=TRUE, ## main numeric vector to be included for derivatives analyses?
	numeric=TRUE, ## are there any numeric vectors to be included for derivatives analyses?
	## fixed effects (var,type,param)
	fixef=list(
		## main
		main=cbind(
			c("age","ns","k=knots")
		),
		## numeric
		numeric=cbind(
			c("ASpErrCorr","lin",""),
			c("AS.lat.corr.AVG","lin",""),
			c("AS.lat.errCorr.AVG","lin",""),
			c("VGS.lat.corr.AVG","lin","")
		),		
		## categorical (or numeric which will be converted to categorical for prediction/derivative analyses)
		categorical=list()
		#	cbind(
		#		c("sex","lin","")
		#	),
		#	cbind(
		#		c("tsr","lin","")
		#	)
		## more nested interactions if necessary go here
		)
	),
	## random effects (var,type,param,grp), if longitudinal
	ranef=cbind(c("age","lin","","id")),
	## coefficients analysis
	coef=list(
		nsim=1000,
		vars=list(main=0,numeric=0,categorical=list()) ## corresponds to  models structure above, add more levels as necessary
	),
	## derivatives analysis (only with main variables and their interactions)
	deriv=list(
		nsim=1000,
		main.range.int=0.1, ## interval width for predicting main
		main.small.int=1e-4, ## tiny interval for derivative calculations
		num.range.int=20, ## how many total intervals to divide numeric variables into for predictions
		vars=list(main=1,numeric=1:4,categorical=list()) ## corresponds to models structure above, add more levels as necessary
	),
  exclude="none",
	custom.code=c(
		"knots<-c(13,15.5,18)",
		"knots<-knots-mean(X$age)",
		"X$age<-X$age-mean(X$age)" ## converting to avoid singularity for correlation of random intercepts/slopes
	)
)
save(setup,file=paste(setup$path,"model.setup",sep="/"))
