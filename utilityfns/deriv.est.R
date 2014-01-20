deriv.est<-function(
	y=Y,
	model=models$X[m,],
	m1.=models$model[[m]],
	m0.=models$model[[models$con[[m]][1]]],
	setup.=setup,
	deriv.do=models$deriv.do[m]
){

	## path
	path<-paste(setup.$path,m,"deriv",sep="/")
	dir.create(path)

	## log for deriv analysis
	sink(paste(path,"deriv.log",sep="/"))
	cat(date(),"DERIVATIVES ANALYSIS\n\n",date(),"setting up...\n\n")

	## gets indices of rows which are excluded in m1; if any excluded, adjust data matrix
	if(length(m1.$all.exc)>0) Y<-as.matrix(y[-m1.$all.exc,]) else Y<-y
	## data frame for m1 (NAs already removed)
	d<-m1.$data
	## updates models with appropriate data frame
	m1<-update(m1.$fit,data=d)
	m0<-update(m0.$fit,data=d)

	## number of simulations
	nsim<-setup.$deriv$nsim
	## small interval for calculation of derivatives
	int<-setup.$deriv$main.small.int
	## interval of predicted values of main variable
	int.r<-setup.$deriv$main.range.int
	## range for prediction of main variable
	range.main<-seq(ceiling(min(d[,1])/int.r)*int.r,floor(max(d[,1])/int.r)*int.r,int.r)
	## range for derivative prediction
	range.d<-rep(range.main,each=2)
	ind<-which((1:length(range.d))%%2==0)
	range.d[ind]<-range.d[ind]+int	
	## range for predictions
	pred=list(range.main)
	names(pred)[1]<-setup.$fixef$main[1,model[1]]
	
	## range for numeric variable
	if(deriv.do==2 | deriv.do==4){
		int.n<-setup.$deriv$num.range.int
		range.num<-range(d[,2])
		range.num<-seq(range.num[1],range.num[2],length.out=int.n)
		pred[[2]]=range.num
		names(pred)[2]<-setup.$fixef$numeric[1,model[2]]
	}

	i.off1<-(setup.$numeric & is.na(model[2]))
	i.off2<-length(pred)

	## other variables
	if(length(setup.$fixef$categorical)>0){
		for(i in (1+i.off2):length(model)){
			if(is.na(model[i+i.off1])) next
			col<-which(names(d)==setup.$fixef$categorical[[i-i.off2]][1,model[i+i.off1]])
			if(class(d[,col])=="numeric") d[,col]<-cont.to.fact(d[,col])
			pred[[i]]<-unique(d[,col])
			names(pred)[i]<-names(d)[col]
		}
	}

	## prediction frame
	if(length(pred)>1) for(i in length(pred):2) if(is.null(pred[[i]])) pred<-pred[-i]
	pred.grid<-expand.grid(pred)
	pred[[1]]<-range.d
	pred.grid.<-expand.grid(pred)
	ind.grid<-which((1:dim(pred.grid)[1])%%2==0)
	ind.grid.<-which((1:dim(pred.grid.)[1])%%2==0)

	## longitudinal
	if(setup.$mixed==TRUE){
		## coefficients
		cat(date(),"estimating models...\n")
		coef<-sapply(
			1:dim(Y)[2],
			function(i){
				cat(i,"")
				refit(m1,Y[,i])@fixef
			}
		)
		cat(date(),"estimation completed...\n\n")

		## predictions
		pred.<-as.matrix(X.(m1,pred.grid.)%*%coef)
		pred<-as.matrix(pred.[ind.grid.-1,])
		pred.dif<-avg.low(pred,pred.grid,range.main)
		pred.d<-as.matrix(diff(pred.)[ind.grid.-1,]/int)
		pred.d.dif<-avg.low(pred.d,pred.grid,range.main)

		## simulate coefficients
		cat(date(),"simulating null coefficients...\n")
		pred.list<-sapply(
			1:dim(Y)[2],
			function(i){
				cat(i,"")
				coef.=coef[,i]
				S=simulate(refit(m0,Y[,i]),nsim)
				s.coef<-sapply(
					1:dim(S)[2],
					function(s) refit(m1,S[,s])@fixef
				)

				## simulation predictions
				pred.sim.<-X.(m1,pred.grid.)%*%s.coef
				pred.sim.d<-diff(pred.sim.)[ind.grid.-1,]/int
				pred.sim<-pred.sim.[ind.grid.-1,]
				pred.sim.m<-apply(pred.sim,1,mean)
				pred.sim.sd<-apply(pred.sim,1,sd)
				pred.p<-1-pnorm(pred[,i],pred.sim.m,pred.sim.sd)
				pred.pboot<-sapply(1:dim(pred.sim)[1], function(c) sum(pred.sim[c,]>pred[c,i])/nsim)
				pred.sim.d.m<-apply(pred.sim.d,1,mean)
				pred.sim.d.sd<-apply(pred.sim.d,1,sd)
				pred.d.p<-1-pnorm(pred.d[,i],pred.sim.d.m,pred.sim.d.sd)
				pred.d.pboot<-sapply(1:dim(pred.sim.d)[1], function(c) sum(pred.sim.d[c,]>pred.d[c,i])/nsim)
				if(deriv.do==1){
					list(pred.sim.m,pred.sim.sd,pred.p,pred.pboot,pred.sim.d.m,pred.sim.d.sd,pred.d.p,pred.d.pboot)
				}else{
					## this whole portion is for calculating the difference from the lower level model, based on x-mean(x) of the last column in the data matrix
					pred.sim.dif<-avg.low(pred.sim,pred.grid,range.main)
					pred.sim.dif.m<-apply(pred.sim.dif,1,mean)
					pred.sim.dif.sd<-apply(pred.sim.dif,1,sd)
					pred.dif.p<-1-pnorm(pred.dif[,i],pred.sim.dif.m,pred.sim.dif.sd)
					pred.dif.pboot<-sapply(1:dim(pred.sim.dif)[1], function(c) sum(pred.sim.dif[c,]>pred.dif[c,i])/nsim)
					pred.sim.d.dif<-avg.low(pred.sim.d,pred.grid,range.main)
					pred.sim.d.dif.m<-apply(pred.sim.d.dif,1,mean)
					pred.sim.d.dif.sd<-apply(pred.sim.d.dif,1,sd)
					pred.d.dif.p<-1-pnorm(pred.d.dif[,i],pred.sim.d.dif.m,pred.sim.d.dif.sd)
					pred.d.dif.pboot<-sapply(1:dim(pred.sim.d.dif)[1], function(c) sum(pred.sim.d.dif[c,]>pred.d.dif[c,i])/nsim)
					list(pred.sim.m,pred.sim.sd,pred.p,pred.pboot,pred.sim.d.m,pred.sim.d.sd,pred.d.p,pred.d.pboot,pred.sim.dif.m,pred.sim.dif.sd,pred.dif.p,pred.dif.pboot,pred.sim.d.dif.m,pred.sim.d.dif.sd,pred.d.dif.p,pred.d.dif.pboot)
				}
			}
		)
		cat(date(),"simulation completed...\n\n")

		
	## cross-sectional
	}else{
		## coefficients
		cat(date(),"estimating models...\n")
		X<-model.matrix(m1)
		coef<-sapply(1:dim(Y)[2],
			function(i){
				cat(i,"")
				lm.fit(X,Y[,i])$coef
			}
		)
		cat(date(),"estimation completed...\n\n")

		## predictions
		pred.<-as.matrix(X.(m1,pred.grid.)%*%coef)
		pred<-as.matrix(pred.[ind.grid.-1,])
		pred.dif<-avg.low(pred,pred.grid,range.main)
		pred.d<-as.matrix(diff(pred.)[ind.grid.-1,]/int)
		pred.d.dif<-avg.low(pred.d,pred.grid,range.main)

		## simulated coefficients
		cat(date(),"simulating null coefficients...\n")
		pred.list<-sapply(
			1:dim(Y)[2],
			function(i){
				cat(i,"")
				coef.=coef[,i]
				S=simulate(update(m0),nsim)
				s.coef<-sapply(1:dim(S)[2], function(s) lm.fit(X,S[,s])$coef)
				## simulation predictions
				pred.sim.<-X.(m1,pred.grid.)%*%s.coef
				pred.sim.d<-diff(pred.sim.)[ind.grid.-1,]/int
				pred.sim<-pred.sim.[ind.grid.-1,]
				pred.sim.m<-apply(pred.sim,1,mean)
				pred.sim.sd<-apply(pred.sim,1,sd)
				pred.p<-1-pnorm(pred[,i],pred.sim.m,pred.sim.sd)
				pred.pboot<-sapply(1:dim(pred.sim)[1], function(c) sum(pred.sim[c,]>pred[c,i])/nsim)
				pred.sim.d.m<-apply(pred.sim.d,1,mean)
				pred.sim.d.sd<-apply(pred.sim.d,1,sd)
				pred.d.p<-1-pnorm(pred.d[,i],pred.sim.d.m,pred.sim.d.sd)
				pred.d.pboot<-sapply(1:dim(pred.sim.d)[1], function(c) sum(pred.sim.d[c,]>pred.d[c,i])/nsim)
				if(deriv.do==1){
					list(pred.sim.m,pred.sim.sd,pred.p,pred.pboot,pred.sim.d.m,pred.sim.d.sd,pred.d.p,pred.d.pboot)
				}else{
					## this whole portion is for calculating the difference from the lower level model, based on x-mean(x) of the last column in the data matrix
					pred.sim.dif<-avg.low(pred.sim,pred.grid,range.main)
					pred.sim.dif.m<-apply(pred.sim.dif,1,mean)
					pred.sim.dif.sd<-apply(pred.sim.dif,1,sd)
					pred.dif.p<-1-pnorm(pred.dif[,i],pred.sim.dif.m,pred.sim.dif.sd)
					pred.dif.pboot<-sapply(1:dim(pred.sim.dif)[1], function(c) sum(pred.sim.dif[c,]>pred.dif[c,i])/nsim)
					pred.sim.d.dif<-avg.low(pred.sim.d,pred.grid,range.main)
					pred.sim.d.dif.m<-apply(pred.sim.d.dif,1,mean)
					pred.sim.d.dif.sd<-apply(pred.sim.d.dif,1,sd)
					pred.d.dif.p<-1-pnorm(pred.d.dif[,i],pred.sim.d.dif.m,pred.sim.d.dif.sd)
					pred.d.dif.pboot<-sapply(1:dim(pred.sim.d.dif)[1], function(c) sum(pred.sim.d.dif[c,]>pred.d.dif[c,i])/nsim)
					list(pred.sim.m,pred.sim.sd,pred.p,pred.pboot,pred.sim.d.m,pred.sim.d.sd,pred.d.p,pred.d.pboot,pred.sim.dif.m,pred.sim.dif.sd,pred.dif.p,pred.dif.pboot,pred.sim.d.dif.m,pred.sim.d.dif.sd,pred.d.dif.p,pred.d.dif.pboot)
				}
			}
		)
		cat(date(),"simulation completed...\n\n")

	}

	if(deriv.do==1) pred.vars<-c("pred.sim.m","pred.sim.sd","pred.p","pred.pboot","pred.sim.d.m","pred.sim.d.sd","pred.d.p","pred.d.pboot") else pred.vars<-c("pred.sim.m","pred.sim.sd","pred.p","pred.pboot","pred.sim.d.m","pred.sim.d.sd","pred.d.p","pred.d.pboot","pred.sim.dif.m","pred.sim.dif.sd","pred.dif.p","pred.dif.pboot","pred.sim.d.dif.m","pred.sim.d.dif.sd","pred.d.dif.p","pred.d.dif.pboot")
	for(i in 1:length(pred.vars)) assign(pred.vars[i],matrix(unlist(pred.list[i,]),length(pred.list[[i,1]]),dim(pred.list)[2]))

	## p-vals and direction of difference
	p<-1-abs(pred.p-0.5)*2
	p.sign<-sign(pred.p-0.5)	
	p.d<-1-abs(pred.d.p-0.5)*2
	p.d.sign<-sign(pred.d.p-0.5)	
	pboot<-1-abs(pred.pboot-0.5)*2
	pboot.sign<-sign(pred.pboot-0.5)	
	pboot.d<-1-abs(pred.d.pboot-0.5)*2
	pboot.d.sign<-sign(pred.d.pboot-0.5)
	if(deriv.do>1){
		dif.p<-1-abs(pred.dif.p-0.5)*2
		dif.p.sign<-sign(pred.dif.p-0.5)	
		dif.p.d<-1-abs(pred.d.dif.p-0.5)*2
		dif.p.d.sign<-sign(pred.d.dif.p-0.5)	
		dif.pboot<-1-abs(pred.dif.pboot-0.5)*2
		dif.pboot.sign<-sign(pred.dif.pboot-0.5)	
		dif.pboot.d<-1-abs(pred.d.dif.pboot-0.5)*2
		dif.pboot.d.sign<-sign(pred.d.dif.pboot-0.5)
	}
	
	## TABLES
	path.tables<-paste(path,"tables",sep="/")
	dir.create(path.tables)
#	source(paste(setup$path.scripts,"utilityfns","deriv.tbl.R",sep="/"))

	## pred.grid
	filename=paste(path.tables,"pred.grid",sep="/")
	sink(filename)
	cat("",names(pred.grid),"\n")
	sink()
	write.table(pred.grid,col.names=FALSE,file=filename,append=TRUE)

	## predicted values
	filename=paste(path.tables,"pred",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pred,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## mean of simulated predictions
	filename=paste(path.tables,"sim.pred.mean",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pred.sim.m,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## sd of simulated predictions
	filename=paste(path.tables,"sim.pred.sd",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pred.sim.sd,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## predicted derivative values
	filename=paste(path.tables,"pred.d",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pred.d,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## mean of simulated derivative predictions
	filename=paste(path.tables,"sim.pred.d.mean",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pred.sim.d.m,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## sd of simulated derivative predictions
	filename=paste(path.tables,"sim.pred.d.sd",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pred.sim.d.sd,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## pred.p
	filename=paste(path.tables,"pred.p",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(p,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## pred.p significance stars
	filename=paste(path.tables,"pred.p.star",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*",""))),
		col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## pred.p sign (-1 for low/1 for high)
	filename=paste(path.tables,"pred.p.sign",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(p.sign,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## pred.d.p
	filename=paste(path.tables,"pred.d.p",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(p.d,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## pred.d.p significance stars
	filename=paste(path.tables,"pred.d.p.star",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(ifelse(p.d<.001,"***",ifelse(p.d<.01,"**",ifelse(p.d<.05,"*",""))),
		col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## pred.d.p sign (-1 for low/1 for high)
	filename=paste(path.tables,"pred.d.p.sign",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(p.d.sign,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## pred.pboot
	filename=paste(path.tables,"pred.pboot",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pboot,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## pred.pboot significance stars
	filename=paste(path.tables,"pred.pboot.star",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(ifelse(pboot<.001,"***",ifelse(pboot<.01,"**",ifelse(pboot<.05,"*",""))),
		col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## pred.pboot sign (-1 for low/1 for high)
	filename=paste(path.tables,"pred.pboot.sign",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pboot.sign,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## pred.d.pboot
	filename=paste(path.tables,"pred.d.pboot",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pboot.d,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## pred.d.pboot significance stars
	filename=paste(path.tables,"pred.d.pboot.star",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(ifelse(pboot.d<.001,"***",ifelse(pboot.d<.01,"**",ifelse(pboot.d<.05,"*",""))),
		col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	## pred.d.p sign (-1 for low/1 for high)
	filename=paste(path.tables,"pred.d.pboot.sign",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pboot.d.sign,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

	if(deriv.do>1){

		## predicted.dif values
		filename=paste(path.tables,"pred.dif",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(pred.dif,col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## mean of simulated predictions.dif
		filename=paste(path.tables,"sim.pred.dif.mean",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(pred.sim.dif.m,col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## sd of simulated predictions.dif
		filename=paste(path.tables,"sim.pred.dif.sd",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(pred.sim.dif.sd,col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## predicted.dif derivative values
		filename=paste(path.tables,"pred.d.dif",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(pred.d.dif,col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## mean of simulated derivative predictions.dif
		filename=paste(path.tables,"sim.pred.d.dif.mean",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(pred.sim.d.dif.m,col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## sd of simulated derivative predictions.dif
		filename=paste(path.tables,"sim.pred.d.dif.sd",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(pred.sim.d.dif.sd,col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## pred.dif.p
		filename=paste(path.tables,"pred.dif.p",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(dif.p,col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## pred.dif.p significance stars
		filename=paste(path.tables,"pred.dif.p.star",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(ifelse(dif.p<.001,"***",ifelse(dif.p<.01,"**",ifelse(dif.p<.05,"*",""))),
			col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## pred.dif.p sign (-1 for low/1 for high)
		filename=paste(path.tables,"pred.dif.p.sign",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(dif.p.sign,col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## pred.d.dif.p
		filename=paste(path.tables,"pred.d.dif.p",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(dif.p.d,col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## pred.d.dif.p significance stars
		filename=paste(path.tables,"pred.d.dif.p.star",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(ifelse(dif.p.d<.001,"***",ifelse(dif.p.d<.01,"**",ifelse(dif.p.d<.05,"*",""))),
			col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## pred.d.dif.p sign (-1 for low/1 for high)
		filename=paste(path.tables,"pred.d.dif.p.sign",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(dif.p.d.sign,col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## pred.dif.pboot
		filename=paste(path.tables,"pred.dif.pboot",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(dif.pboot,col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## pred.dif.pboot significance stars
		filename=paste(path.tables,"pred.dif.pboot.star",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(ifelse(dif.pboot<.001,"***",ifelse(dif.pboot<.01,"**",ifelse(dif.pboot<.05,"*",""))),
			col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## pred.dif.pboot sign (-1 for low/1 for high)
		filename=paste(path.tables,"pred.dif.pboot.sign",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(dif.pboot.sign,col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## pred.d.dif.pboot
		filename=paste(path.tables,"pred.d.dif.pboot",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(dif.pboot.d,col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## pred.d.dif.pboot significance stars
		filename=paste(path.tables,"pred.d.dif.pboot.star",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(ifelse(dif.pboot.d<.001,"***",ifelse(dif.pboot.d<.01,"**",ifelse(dif.pboot.d<.05,"*",""))),col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

		## pred.d.dif.p sign (-1 for low/1 for high)
		filename=paste(path.tables,"pred.d.dif.pboot.sign",sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(dif.pboot.d.sign,col.names=FALSE,row.names=range.main,file=filename,append=TRUE)

	}
	
	cat(date(),"derivatives analysis completed...\n\n")
	sink()
	
	"deriv.est() completed"
}


## EXTRA FUNCTIONS

## function for getting design matrix for prediction
X.<-function(object,newdata){
	tt<-terms(object)
	Terms<-delete.response(tt)
	m<-model.frame(Terms,newdata)
	X<-model.matrix(Terms,m)
	X
}

avg.low<-function(x,grid=pred.grid,r=range.main){
	r<-length(r)
	u<-unique(grid[,dim(grid)[2]])
	ind<-dim(x)[1]
	ind.low<-ind/length(u)
	x<-as.matrix(x)
	x.low<-sapply(1:ind.low, function(i){
		if(i==ind.low) i.<-0 else i.<-i
		colMeans(as.matrix(x[which(1:ind%%ind.low==i.),]))
	})

	if(!is.matrix(x.low)) x.low<-as.matrix(x.low) else x.low<-t(x.low)
	Xx<-x.low
	for(i in 2:length(u)) Xx<-rbind(Xx,x.low)
	dif<-x-Xx
	dif.low<-sapply(1:r, function(i){
		if(i==r) i.<-0 else i.<-i
		colSums(as.matrix(dif[which(1:ind%%r==i.),]^2))
	})
	if(!is.matrix(dif.low)) as.matrix(dif.low) else t(dif.low)
}
