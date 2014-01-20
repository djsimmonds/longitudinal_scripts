## setup
setup<-list(
	path="/home/danisimmonds/Dani/dti_0511/tbss/analysis/062011",
	path.scripts="/home/danisimmonds/Dropbox/scripts0611",
	ynames=c("all","proj","assoc","assoc.l","cal","cer.c","cer.p","par","occ","sm","fron","temp","mt","bg","thal","ML","PCT","CST","CP","IC.A","IC.P","IC.R","CR.A","CR.S","CR.P","PTR","SFOF","SS","EC","UF","SLF","FOR.CB","FOR.C","CIN.CG","CIN.H","CAL.G","CAL.B","CAL.S","CAL.T","CER.I","CER.M","CER.S"),
	mixed=TRUE, ## longitudinal?
	main=TRUE, ## main numeric vector to be included for derivatives analyses?
	numeric=TRUE, ## are there any numeric vectors to be included for derivatives analyses?
	## fixed effects (var,type,param)
	fixef=list(
		## main
		main=cbind(
			c("age","lin",""), ## linear fit
			c("age","ns","k=knots") ## roughly divides into childhood (<13), early adolescence (13-15.5), late adolescence (15.5-18) and adulthood (>18)
		),
		## numeric
		numeric=cbind(
			c("viq","lin",""),
			c("piq","lin",""),
			c("vgs.mRT","lin",""),
			c("vgs.sdRT","lin",""),
			c("vgs.cv","lin",""),
			c("vgs.mu","lin",""),
			c("vgs.sigma","lin",""),
			c("vgs.tau","lin",""),
			c("vgs.slow4","lin",""),
			c("anti.percErr","lin",""),
			c("anti.mRT","lin",""),
			c("anti.sdRT","lin",""),
			c("anti.cv","lin",""),
			c("anti.mu","lin",""),
			c("anti.sigma","lin",""),
			c("anti.tau","lin",""),
			c("anti.slow4corr","lin",""),
			c("anti.slow4all","lin","")
		),		
		## categorical (or numeric which will be converted to categorical for prediction/derivative analyses)
		categorical=list(
			cbind(
				c("sex","lin","")
			),
			cbind(
				c("tsr","lin","")
			)
		## more nested interactions if necessary go here
		)
	),
	## random effects (var,type,param,grp), if longitudinal
	ranef=cbind(c("age","lin","","id")),
	## coefficients analysis
	coef=list(
		nsim=1000,
		vars=list(main=1,numeric=1:18,categorical=list(1,1)) ## corresponds to  models structure above, add more levels as necessary
	),
	## derivatives analysis (only with main variables and their interactions)
	deriv=list(
		nsim=1000,
		main.range.int=0.1, ## interval width for predicting main
		main.small.int=1e-4, ## tiny interval for derivative calculations
		num.range.int=20, ## how many total intervals to divide numeric variables into for predictions
		vars=list(main=2,numeric=1:18,categorical=list(1,1)) ## corresponds to models structure above, add more levels as necessary
	),
	custom.code=c(
		"knots<-c(13,15.5,18)",
		"knots<-knots-mean(X$age)",
		"X$age<-X$age-mean(X$age)" ## converting to avoid singularity for correlation of random intercepts/slopes
	)
)
save(setup,file=paste(setup$path,"model.setup",sep="/"))
