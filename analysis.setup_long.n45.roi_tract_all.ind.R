## setup
setup<-list(
	path="/home/danisimmonds/Dani/dti_0511/tbss/analysis/long.n45/roi_tract_all/ind",
	path.scripts="/home/danisimmonds/Dropbox/scripts2",
	ynames=c("CCG","CCB","CCS","CCT","ICP","MCP","SCP","PCT","ML","CST","CP","IC","CR","PTR","SS","EC","SFOF","UF","SLF","CG","FOR"),

	mixed=TRUE, ## longitudinal?
	main=TRUE, ## main numeric vector to be included for derivatives analyses?
	numeric=FALSE, ## are there any numeric vectors to be included for derivatives analyses?
	#numeric.weights=c( ## weights files for numeric variables
	#	"/home/danisimmonds/Dani/dti_0511/demographics/beh.out/vgs.mRT_n45.out_weights",
	#	"/home/danisimmonds/Dani/dti_0511/demographics/beh.out/vgs.cv_n45.out_weights",
	#	"/home/danisimmonds/Dani/dti_0511/demographics/beh.out/anti.percErr_n45.out_weights",
	#	"/home/danisimmonds/Dani/dti_0511/demographics/beh.out/anti.cv_n45.out_weights"
	#),
	## fixed effects (var,type,param)
	fixef=list(
		## main
		main=cbind(
			c("age","lin",""), ## linear fit
			#c("age","poly","2"), ## quadratic fit
			c("age","ns","k=c(5,7.5,10)") ## roughly divides into childhood (<13), early adolescence (13-15.5), late adolescence (15.5-18), and adulthood (>18)
		)#,
		## numeric
		#numeric=cbind(
		#	c("vgs.mRT","lin",""),
		#	c("vgs.cv","lin",""),
		#	c("anti.percErr","lin",""),
		#	c("anti.cv","lin","")
		#),		
		## categorical (or numeric which will be converted to categorical for prediction/derivative analyses)
		#categorical=list(
		#	cbind(
		#		c("sex","lin","")
		#	)
		## more nested interactions if necessary go here
		#)
	),
	## random effects (var,type,param,grp), if longitudinal
	ranef=cbind(c("age","lin","","id")),
	## derivatives analysis (only with main variables and their interactions)
	deriv=list(
		nsim=1000,
		main.range.int=0.1, ## interval width for predicting main
		main.small.int=1e-4, ## tiny interval for derivative calculations
		num.range.int=20, ## how many total intervals to divide numeric variables into for predictions
		vars=list(main=2),#,numeric=1:4,categorical=list(1)), ## corresponds to models structure above, add more levels as necessary
		p.thr=1 ## no threshold, always run derivative analysis
	)
)
save(setup,file=paste(setup$path,"model.setup",sep="/"))
