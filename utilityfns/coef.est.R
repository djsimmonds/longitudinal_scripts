coef.est<-function(
	y=Y,
	m1.=models$model[[m]],
	m0.=models$model[[1]],
	path=paste(setup$path,m,sep="/"),
	mixed=setup$mixed,
	ynames=setup$ynames,
	nsim=setup$coef$nsim
){

	path<-paste(path,"coef",sep="/")
	dir.create(path)

	## gets indices of rows which are excluded in m1; if any excluded, adjust data matrix
	if(length(m1.$all.exc)>0) Y<-as.matrix(y[-m1.$all.exc,]) else Y<-y
	## data frame for m1 (NAs already removed)
	d<-m1.$data
	## updates models with appropriate data frame
	m0<-update(m0.$fit,data=d)
	m1<-update(m1.$fit,data=d)

	## estimates

	## longitudinal
	if(mixed==TRUE){
		xnames<-names(m1@fixef)
		coef<-sapply(
			1:dim(Y)[2],
			function(i) refit(m1,Y[,i])@fixef
		)
		coef.list<-sapply(
			1:dim(Y)[2],
			function(i){
				coef.=coef[,i]
				S=simulate(refit(m0,Y[,i]),nsim)
				s.coef<-sapply(1:dim(S)[2], function(s) refit(m1,S[,s])@fixef)
				s.coef.m<-apply(s.coef,1,mean)
				s.coef.sd<-apply(s.coef,1,sd)
				coef.p<-pnorm(coef.,s.coef.m,s.coef.sd)
				coef.pboot<-sapply(1:length(coef.), function(c) sum(s.coef[c,]>coef.[c])/nsim)
				list(s.coef.m,s.coef.sd,coef.p,coef.pboot)		
			}
		)
		
	## cross-sectional
	}else{
		xnames<-names(m1$coef)
		X<-model.matrix(m1)
		coef<-sapply(1:dim(Y)[2], function(i) lm.fit(X,Y[,i])$coef)
		coef.list<-sapply(
			1:dim(Y)[2],
			function(i){
				coef.<-coef[,i]
				S=simulate(m0,nsim)
				s.coef<-sapply(1:dim(S)[2], function(s) lm.fit(X,S[,s])$coef)
				s.coef.m<-apply(s.coef,1,mean)
				s.coef.sd<-apply(s.coef,1,sd)
				coef.p<-pnorm(coef.,s.coef.m,s.coef.sd)
				coef.pboot<-sapply(1:length(coef.), function(c) sum(s.coef[c,]>coef.[c])/nsim)
				list(s.coef.m,s.coef.sd,coef.p,coef.pboot)	
			}
		)
	}
	s.coef.m<-sapply(1:dim(coef.list)[2], function(i) coef.list[[1,i]])
	s.coef.sd<-sapply(1:dim(coef.list)[2], function(i) coef.list[[2,i]])
	coef.p<-sapply(1:dim(coef.list)[2], function(i) coef.list[[3,i]])
	coef.pboot<-sapply(1:dim(coef.list)[2], function(i) coef.list[[4,i]])

	## p-vals and direction of difference
	p<-1-abs(coef.p-0.5)*2
	p.sign<-sign(coef.p-0.5)	
	pboot<-1-abs(coef.pboot-0.5)*2
	pboot.sign<-sign(coef.pboot-0.5)	
	
	## coefficients
	filename=paste(path,"coef",sep="/")
	sink(filename)
	cat("",ynames,"\n")
	sink()
	write.table(coef,col.names=FALSE,row.names=xnames,file=filename,append=TRUE)

	## mean of simulated coefficients
	filename=paste(path,"sim.mean",sep="/")
	sink(filename)
	cat("",ynames,"\n")
	sink()
	write.table(s.coef.m,col.names=FALSE,row.names=xnames,file=filename,append=TRUE)

	## sd of simulated coefficients
	filename=paste(path,"sim.sd",sep="/")
	sink(filename)
	cat("",ynames,"\n")
	sink()
	write.table(s.coef.sd,col.names=FALSE,row.names=xnames,file=filename,append=TRUE)

	## p-values (pnorm)
	filename=paste(path,"p",sep="/")
	sink(filename)
	cat("",ynames,"\n")
	sink()
	write.table(p,col.names=FALSE,row.names=xnames,file=filename,append=TRUE)

	## significance stars (pnorm)
	filename=paste(path,"p.star",sep="/")
	sink(filename)
	cat("",ynames,"\n")
	sink()
	write.table(ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*",""))),
		col.names=FALSE,row.names=xnames,file=filename,append=TRUE)

	## sign of difference (-1 of significantly lower, 1 if significantly higher) (pnorm)
	filename=paste(path,"p.sign",sep="/")
	sink(filename)
	cat("",ynames,"\n")
	sink()
	write.table(p.sign,col.names=FALSE,row.names=xnames,file=filename,append=TRUE)
	
	## p-values (boot)
	filename=paste(path,"pboot",sep="/")
	sink(filename)
	cat("",ynames,"\n")
	sink()
	write.table(pboot,col.names=FALSE,row.names=xnames,file=filename,append=TRUE)

	## significance stars (boot)
	filename=paste(path,"pboot.star",sep="/")
	sink(filename)
	cat("",ynames,"\n")
	sink()
	write.table(ifelse(pboot<.001,"***",ifelse(pboot<.01,"**",ifelse(pboot<.05,"*",""))),
		col.names=FALSE,row.names=xnames,file=filename,append=TRUE)

	## sign of difference (-1 of significantly lower, 1 if significantly higher) (boot)
	filename=paste(path,"pboot.sign",sep="/")
	sink(filename)
	cat("",ynames,"\n")
	sink()
	write.table(pboot.sign,col.names=FALSE,row.names=xnames,file=filename,append=TRUE)
	
	"coef.est() completed"
}
