## setup
setup<-list(
	path="/home/danisimmonds/Dani/dti_0511/tbss/analysis/vox",
	path.scripts="/home/danisimmonds/Dropbox/scripts0611",
	ynames=1:24880, ## number of vox
	mixed=TRUE, ## longitudinal?
	main=TRUE, ## main numeric vector to be included for derivatives analyses?
	numeric=FALSE, ## are there any numeric vectors to be included for derivatives analyses?
	## fixed effects (var,type,param)
	fixef=list(
		## main
		main=cbind(
			c("age","ns","k=knots") ## roughly divides into childhood (<13), early adolescence (13-15.5), late adolescence (15.5-18) and adulthood (>18)
		),
		## numeric
		numeric=list(),		
		## categorical (or numeric which will be converted to categorical for prediction/derivative analyses)
		categorical=list()
	),
	## random effects (var,type,param,grp), if longitudinal
	ranef=cbind(c("age","lin","","id")),
	## coefficients analysis
	coef=list(
		nsim=1000,
		vars=list(main=1,numeric=NA,categorical=NA) ## corresponds to  models structure above, add more levels as necessary
	),
	## derivatives analysis (only with main variables and their interactions)
	deriv=list(
		nsim=1000,
		main.range.int=0.1, ## interval width for predicting main
		main.small.int=1e-4, ## tiny interval for derivative calculations
		num.range.int=20, ## how many total intervals to divide numeric variables into for predictions
		vars=list(main=1,numeric=NA,categorical=NA) ## corresponds to models structure above, add more levels as necessary
	),
	custom.code=c(
		"knots<-c(13,15.5,18)",
		"knots<-knots-mean(X$age)",
		"X$age<-X$age-mean(X$age)" ## converting to avoid singularity for correlation of random intercepts/slopes
	)
)
save(setup,file=paste(setup$path,"model.setup",sep="/"))
