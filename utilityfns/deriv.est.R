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
	dir.create(path,recursive=TRUE)

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

	## tables to print
	tbl.names<-list(
		all=c("pred.sim.m","pred.sim.sd","pred.p","pred.p.sign","pred.p.star","pred.pboot","pred.pboot.sign","pred.pboot.star","pred.sim.d.m","pred.sim.d.sd","pred.d.p","pred.d.p.sign","pred.d.p.star","pred.d.pboot","pred.d.pboot.sign","pred.d.pboot.star"),
		cor=c("pred.sim.cor.m","pred.sim.cor.sd","pred.cor.p","pred.cor.p.sign","pred.cor.p.star","pred.cor.pboot","pred.cor.pboot.sign","pred.cor.pboot.star","pred.sim.d.cor.m","pred.sim.d.cor.sd","pred.d.cor.p","pred.d.cor.p.sign","pred.d.cor.p.star","pred.d.cor.pboot","pred.d.cor.pboot.sign","pred.d.cor.pboot.star"),
		dif=c("pred.sim.dif.m","pred.sim.dif.sd","pred.dif.p","pred.dif.p.sign","pred.dif.p.star","pred.dif.pboot","pred.dif.pboot.sign","pred.dif.pboot.star","pred.sim.d.dif.m","pred.sim.d.dif.sd","pred.d.dif.p","pred.d.dif.p.sign","pred.d.dif.p.star","pred.d.dif.pboot","pred.d.dif.pboot.sign","pred.d.dif.pboot.star")
	)

	tbl.listStr<-paste("list",tbl.names$all[1],sep="(")
	tbl.names.<-tbl.names$all
	for(i in 2:length(tbl.names$all)) tbl.listStr<-paste(tbl.listStr,tbl.names$all[i],sep=",")
	if(deriv.do==2|deriv.do==4){
		for(i in 1:length(tbl.names$cor)) tbl.listStr<-paste(tbl.listStr,tbl.names$cor[i],sep=",")
		tbl.names.<-c(tbl.names.,tbl.names$cor)
	}
	if(deriv.do>2){
		for(i in 1:length(tbl.names$dif)) tbl.listStr<-paste(tbl.listStr,tbl.names$dif[i],sep=",")
		tbl.names.<-c(tbl.names.,tbl.names$dif)		
	}
	tbl.listStr=paste(tbl.listStr,")",sep="")

	## prediction frame
	for(i in length(pred):1) if(is.null(pred[[i]])) pred<-pred[-i]
	pred.grid<-expand.grid(pred)
	pred[[1]]<-range.d
	pred.grid.<-expand.grid(pred)
	ind.grid<-which((1:dim(pred.grid)[1])%%2==0)
	ind.grid.<-which((1:dim(pred.grid.)[1])%%2==0)

	## longitudinal
	if(setup.$mixed==TRUE){
		## coefficients
		coef<-sapply(
			1:dim(Y)[2],
			function(i) refit(m1,Y[,i])@fixef
		)

		## predictions
		pred.<-as.matrix(X.(m1,pred.grid.)%*%coef)
		pred<-as.matrix(pred.[ind.grid.-1,])
		pred.d<-as.matrix(diff(pred.)[ind.grid.-1,]/int)
		if(deriv.do==2 | deriv.do==4){
			u.var<-unique(pred.grid[,2])
			ind.v1<-which(pred.grid[,2]==u.var[1])
			ind.v2<-which(pred.grid[,2]==u.var[2])
			pred.cor<-as.matrix((pred[ind.v2,]-pred[ind.v1,])/(u.var[2]-u.var[1]))
			pred.d.cor<-as.matrix((pred.d[ind.v2,]-pred.d[ind.v1,])/(u.var[2]-u.var[1]))
			if(deriv.do==4){
				pred.dif<-avg.low(pred.cor,pred.grid[ind.v1,-2],range.main)
				pred.d.dif<-avg.low(pred.d.cor,pred.grid[ind.v1,-2],range.main)
			}
		}else if(deriv.do==3 | deriv.do==5){
			pred.dif<-avg.low(pred,pred.grid,range.main)
			pred.d.dif<-avg.low(pred.d,pred.grid,range.main)
		}

		## simulate coefficients
		pred.list<-sapply(
			1:dim(Y)[2],
			function(i){
				coef.=coef[,i]
				S=simulate(refit(m0,Y[,i]),nsim)
				s.coef<-sapply(
					1:dim(S)[2],
					function(s) refit(m1,S[,s])@fixef
				)

				## simulation predictions and p-values
				pred.sim.<-X.(m1,pred.grid.)%*%s.coef
				pred.sim.d<-diff(pred.sim.)[ind.grid.-1,]/int
				pred.sim<-pred.sim.[ind.grid.-1,]
				pred.sim.m<-apply(pred.sim,1,mean)
				pred.sim.sd<-apply(pred.sim,1,sd)
				pred.p<-2*(pnorm(pred[,i],pred.sim.m,pred.sim.sd)-0.5)
				pred.p.sign<-sign(pred.p)
				pred.p<-1-abs(pred.p)
				pred.p.star<-sig.star(pred.p)
				pred.pboot<-sapply(1:dim(pred.sim)[1], function(c) 2*(0.5-sum(pred.sim[c,]>pred[c,i])/nsim))
				pred.pboot.sign<-sign(pred.pboot)
				pred.pboot<-1-abs(pred.pboot)
				pred.pboot.star<-sig.star(pred.pboot)
				pred.sim.d.m<-apply(pred.sim.d,1,mean)
				pred.sim.d.sd<-apply(pred.sim.d,1,sd)
				pred.d.p<-2*(pnorm(pred.d[,i],pred.sim.d.m,pred.sim.d.sd)-0.5)
				pred.d.p.sign<-sign(pred.d.p)
				pred.d.p<-1-abs(pred.d.p)
				pred.d.p.star<-sig.star(pred.d.p)
				pred.d.pboot<-sapply(1:dim(pred.sim.d)[1], function(c) 2*(0.5-sum(pred.sim.d[c,]>pred.d[c,i])/nsim))
				pred.d.pboot.sign<-sign(pred.d.pboot)
				pred.d.pboot<-1-abs(pred.d.pboot)
				pred.d.pboot.star<-sig.star(pred.d.pboot)

				if(deriv.do==2 | deriv.do==4){
					pred.sim.cor<-sapply(1:nsim, function(s) (pred.sim[ind.v2,s]-pred.sim[ind.v1,s])/(u.var[2]-u.var[1]))
					pred.sim.cor.m<-apply(pred.sim.cor,1,mean)
					pred.sim.cor.sd<-apply(pred.sim.cor,1,sd)
					pred.cor.p<-2*(pnorm(pred.cor[,i],pred.sim.cor.m,pred.sim.cor.sd)-0.5)
					pred.cor.p.sign<-sign(pred.cor.p)
					pred.cor.p<-1-abs(pred.cor.p)
					pred.cor.p.star<-sig.star(pred.cor.p)
					pred.cor.pboot<-sapply(1:dim(pred.sim.cor)[1], function(c) 2*(0.5-sum(pred.sim.cor[c,]>pred.cor[c,i])/nsim))
					pred.cor.pboot.sign<-sign(pred.cor.pboot)
					pred.cor.pboot<-1-abs(pred.cor.pboot)
					pred.cor.pboot.star<-sig.star(pred.cor.pboot)
					pred.sim.d.cor<-sapply(1:nsim, function(s) (pred.sim.d[ind.v2,s]-pred.sim.d[ind.v1,s])/(u.var[2]-u.var[1]))
					pred.sim.d.cor.m<-apply(pred.sim.d.cor,1,mean)
					pred.sim.d.cor.sd<-apply(pred.sim.d.cor,1,sd)
					pred.d.cor.p<-2*(pnorm(pred.d.cor[,i],pred.sim.d.cor.m,pred.sim.d.cor.sd)-0.5)
					pred.d.cor.p.sign<-sign(pred.d.cor.p)
					pred.d.cor.p<-1-abs(pred.d.cor.p)
					pred.d.cor.p.star<-sig.star(pred.d.cor.p)
					pred.d.cor.pboot<-sapply(1:dim(pred.sim.d.cor)[1], function(c) 2*(0.5-sum(pred.sim.d.cor[c,]>pred.d.cor[c,i])/nsim))
					pred.d.cor.pboot.sign<-sign(pred.d.cor.pboot)
					pred.d.cor.pboot<-1-abs(pred.d.cor.pboot)
					pred.d.cor.pboot.star<-sig.star(pred.d.cor.pboot)
				}
				
				if(deriv.do>2){
					if(deriv.do==4){
						pred.sim.dif<-sapply(1:nsim, function(s) avg.low(as.matrix(pred.sim.cor[,s]),pred.grid[ind.v1,-2],range.main))
						pred.sim.d.dif<-sapply(1:nsim, function(s) avg.low(as.matrix(pred.sim.d.cor[,s]),pred.grid[ind.v1,-2],range.main))
					}else{
						pred.sim.dif<-avg.low(pred.sim,pred.grid,range.main)
						pred.sim.d.dif<-avg.low(pred.sim.d,pred.grid,range.main)
					}
					pred.sim.dif.m<-apply(pred.sim.dif,1,mean)
					pred.sim.dif.sd<-apply(pred.sim.dif,1,sd)
					pred.dif.p<-2*(pnorm(pred.dif[,i],pred.sim.dif.m,pred.sim.dif.sd)-0.5)
					pred.dif.p.sign<-sign(pred.dif.p)
					pred.dif.p<-1-abs(pred.dif.p)
					pred.dif.p.star<-sig.star(pred.dif.p)
					pred.dif.pboot<-sapply(1:dim(pred.sim.dif)[1], function(c) 2*(0.5-sum(pred.sim.dif[c,]>pred.dif[c,i])/nsim))
					pred.dif.pboot.sign<-sign(pred.dif.pboot)
					pred.dif.pboot<-1-abs(pred.dif.pboot)
					pred.dif.pboot.star<-sig.star(pred.dif.pboot)
					pred.sim.d.dif.m<-apply(pred.sim.d.dif,1,mean)
					pred.sim.d.dif.sd<-apply(pred.sim.d.dif,1,sd)
					pred.d.dif.p<-2*(pnorm(pred.d.dif[,i],pred.sim.d.dif.m,pred.sim.d.dif.sd)-0.5)
					pred.d.dif.p.sign<-sign(pred.d.dif.p)
					pred.d.dif.p<-1-abs(pred.d.dif.p)
					pred.d.dif.p.star<-sig.star(pred.d.dif.p)
					pred.d.dif.pboot<-sapply(1:dim(pred.sim.d.dif)[1], function(c) 2*(0.5-sum(pred.sim.d.dif[c,]>pred.d.dif[c,i])/nsim))
					pred.d.dif.pboot.sign<-sign(pred.d.dif.pboot)
					pred.d.dif.pboot<-1-abs(pred.d.dif.pboot)
					pred.d.dif.pboot.star<-sig.star(pred.d.dif.pboot)
				}

				eval(parse(text=tbl.listStr))
			}
		)
		
	## cross-sectional
	}else{
		## coefficients
		X<-model.matrix(m1)
		coef<-sapply(1:dim(Y)[2], function(i) lm.fit(X,Y[,i])$coef)

		## predictions
		pred.<-as.matrix(X.(m1,pred.grid.)%*%coef)
		pred<-as.matrix(pred.[ind.grid.-1,])
		pred.dif<-avg.low(pred,pred.grid,range.main)
		pred.d<-as.matrix(diff(pred.)[ind.grid.-1,]/int)
		pred.d.dif<-avg.low(pred.d,pred.grid,range.main)

		## simulated coefficients
		pred.list<-sapply(
			1:dim(Y)[2],
			function(i){
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

## FIX CODE HERE AS WITH LONGITUDINAL ABOVE

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
	}

	for(i in 1:length(tbl.names.)) assign(tbl.names.[i],matrix(unlist(pred.list[i,]),length(pred.list[[i,1]]),dim(pred.list)[2]))
	tbl.names$all<-c("pred","pred.d",tbl.names$all)
	tbl.names$cor<-c("pred.cor","pred.d.cor",tbl.names$cor)
	tbl.names$dif<-c("pred.dif","pred.d.dif",tbl.names$dif)

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

	## write all tables
	for(i in 1:length(tbl.names$all)){
		filename=paste(path.tables,tbl.names$all[i],sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(get(tbl.names$all[i]),col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)
	}
	if(deriv.do==2|deriv.do==4) for(i in 1:length(tbl.names$cor)){
		filename=paste(path.tables,tbl.names$cor[i],sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(get(tbl.names$cor[i]),col.names=FALSE,row.names=apply(as.matrix(pred.grid[ind.v1,-2]),1,paste,collapse=","),file=filename,append=TRUE)
	}
	if(deriv.do>2) for(i in 1:length(tbl.names$dif)){
		filename=paste(path.tables,tbl.names$dif[i],sep="/")
		sink(filename)
		cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
		sink()
		write.table(get(tbl.names$dif[i]),col.names=FALSE,row.names=range.main,file=filename,append=TRUE)
	}

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

sig.star<-function(x) ifelse(x<.001,"***",ifelse(x<.01,"**",ifelse(x<.05,"*","")))

