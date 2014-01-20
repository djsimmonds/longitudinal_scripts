## utility functions for formatting variables, formula strings and estimating models
fixef.set<-function(X=models$X[m,]),fixef=setup$fixef){
	vars<-which(!is.na(X))
	f<-matrix("",3,length(vars))
	for(i in 1:length(vars)) f[,i]<-fixef[[vars[i]]][,X[vars[i]]]
	f
}

## x = c(var,type,param)
pred.str<-function(x){
	if(x[2]=="null") "1" else
	if(x[2]=="lin") x[1] else
	paste(x[2],"(",x[1],",",x[3],")",sep="")
}

lmer.est<-function(f=cbind(c("","null","")),r=setup$ranef,y=setup$resp,data=get(setup$data)){
	na<-numeric(0)
	## fixef
	f.<-pred.str(f[,1])
	if(nchar(f[1,1])>0) na<-c(na,which(is.na(data[,which(names(data)==f[1,1])])))
	if(dim(f)[2]>1){
		for(i in 2:dim(f)[2]){
			f.<-paste(f.,"*",pred.str(f[,i]),sep="")
			na<-c(na,which(is.na(data[,which(names(data)==f[1,i])])))
		}
	}
	## ranef
	r.<-paste("(",pred.str(r[1:3,1]),"|",r[4,1],")",sep="")
	if(dim(r)[2]>1) for(i in 2:dim(r)[2]) r.<-paste(r.,"+","(",pred.str(r[1:3,i]),"|",r[4,i],")",sep="")
	## formula
	formula<-paste(y,"~",f.,"+",r.,sep="")
	if(length(na)>0) data=data[-unique(na),]
	list(fixef=f,formula=formula,na=unique(na),fit=lmer(formula,data,weights=data$wt))
}

## check against setup file to see if derivative analysis should be performed
deriv.check<-function(vars,p,p.thr=setup$deriv$p.thr,inc=setup$deriv$vars){
	## need to have main variable for analysis
	if(is.na(vars[1])) return(FALSE)
	## only perform analysis if log-likelihood test is significant
	if(p>p.thr) return(FALSE)
	## all variables should be included in derivatives list in setup file
	for(i in 1:length(vars)) if(!is.na(vars[i]) & length(which(inc[[i]]==vars[i]))==0) return(FALSE)
	## if all conditions are met, do analysis
	TRUE
}

## convert continuous variable to factor
cont.to.fact<-function(var,data,cut=c(1/3,2/3)){
	col<-which(names(data)==var)
	var.new<-numeric(length(data[,col]))+1
	for(i in 1:length(cut)) var.new<-var.new+ifelse(data[,col]>=quantile(data[,col],cut[i],na.rm=TRUE),1,0)
	var.new<-ifelse(is.na(data[,col]),NA,var.new)
	ord<-character(length(cut)+1)
	ord[1]<-paste(0,round(cut[1]*100,0),sep="-")
	ord[length(ord)]<-paste(round(cut[length(cut)]*100,0),100,sep="-")	
	if(length(ord)>2) for(i in 1:(length(cut)-1)) ord[i+1]<-paste(round(cut[i]*100,0),round(cut[i+1]*100,0),sep="-")
	ord<-ordered(ord)
	data[,col]<-ord[var.new]
	data
}

## significance stars for p values
p.stars<-function(p) if(p<=0.001|p>=0.999) "***" else if(p<=0.01|p>=0.99) "**" else if(p<=0.05|p>=0.95) "*" else if(p<=0.1|p>=0.9) "+" else ""

## derivatives analysis
deriv.est<-function(model,model.sim,path,deriv=setup$deriv,data=get(setup$data),resp=setup$resp){
	## simulate from lower level model
	sim<-simulate(model.sim$fit,deriv$nsim)
	## create new data structure for predictions of bootstrapped models
	pred<-list(names=model$fixef[1,1],values=list(deriv$range))
	if(dim(model$fixef)[2]>1){
		for(i in 2:dim(model$fixef)[2]){
			col<-which(names(data)==model$fixef[1,i])
			## if there is a continuous variable, convert to factor before refitting models (default divides from 0-33, 33-66, 66-100)
			if(class(data[,col])=="numeric") data<-cont.to.fact(model$fixef[1,i],data)
			pred$names[i]<-model$fixef[1,i]
			pred$values[[i]]<-unique(factor(data[,col]))
		}
	}
	## refits model with new data frame replacing numeric variables with factors
	model$fit<-refit(model$fit,data[,which(names(data)==resp)])
	names(pred$values)<-pred$names
	pred.grid<-expand.grid(pred$values)
	for(i in 1:length(pred$names)) assign(pred$names[i],pred.grid[,i],envir=.GlobalEnv) ## necessary for predictSE.mer, has to do with the use of functions like ns() within formulas within functions, without having the variable in the function outside the ns() structure --> i don't fully understand, but somehow by declaring the variables in the main environment, it gets around the problem
	## predictions and calculations of derivatives from model
	pred.0<-predictSE.mer(model$fit,pred.grid,se.fit=FALSE)
	pred.0.sim<-matrix(NA,length(pred.0),deriv$nsim)
	cat("pred.0.sim: ")
	for(i in 1:deriv$nsim){
		cat(i,"")
		pred.0.sim[,i]<-predictSE.mer(refit(model$fit,sim[,i]),pred.grid,se.fit=FALSE)
	}
	assign(pred$names[1],get(pred$names[1])+deriv$int,envir=.GlobalEnv) ## see note above
	pred.1<-predictSE.mer(model$fit,pred.grid,se.fit=FALSE)
	pred.1.sim<-matrix(NA,length(pred.1),deriv$nsim)
	cat("\npred.1.sim: ")
	for(i in 1:deriv$nsim){
		cat(i,"")
		pred.1.sim[,i]<-predictSE.mer(refit(model$fit,sim[,i]),pred.grid,se.fit=FALSE)
	}
	assign(pred$names[1],get(pred$names[1])+deriv$int,envir=.GlobalEnv) ## see note above
	pred.2<-predictSE.mer(model$fit,pred.grid,se.fit=FALSE)
	pred.2.sim<-matrix(NA,length(pred.2),deriv$nsim)
	cat("\npred.2.sim: ")
	for(i in 1:deriv$nsim){
		cat(i,"")
		pred.2.sim[,i]<-predictSE.mer(refit(model$fit,sim[,i]),pred.grid,se.fit=FALSE)
	}	
	## first derivative
	pred.d1.0<-(pred.1-pred.0)/deriv$int
	pred.d1.0.sim<-(pred.1.sim-pred.0.sim)/deriv$int
	pred.d1.1<-(pred.2-pred.1)/deriv$int
	pred.d1.1.sim<-(pred.2.sim-pred.1.sim)/deriv$int
	## second derivative
	pred.d2<-(pred.d1.0-pred.d1.1)/deriv$int
	pred.d2.sim<-(pred.d1.0.sim-pred.d1.1.sim)/deriv$int
	## mean and sd calculations
	## if L>1, average across last dimension for plotting and distance calculations
	if(dim(model$fixef)[2]==1){
		last<-numeric(0)
		mid<-numeric(0)
		## mean of simulations
		pred.sim.mean<-rowSums(pred.0.sim)/dim(pred.0.sim)[2]
		pred.d1.sim.mean<-rowSums(pred.d1.0.sim)/dim(pred.d1.0.sim)[2]
		pred.d2.sim.mean<-rowSums(pred.d2.sim)/dim(pred.d2.sim)[2]
		## sd of simulations
		pred.sim.sd<-numeric(length(pred.sim.mean))
		pred.d1.sim.sd<-numeric(length(pred.sim.mean))
		pred.d2.sim.sd<-numeric(length(pred.sim.mean))
		for(i in 1:dim(pred.grid)[1]){
			pred.sim.sd[i]<-sd(pred.0.sim[i,])
			pred.d1.sim.sd[i]<-sd(pred.d1.0.sim[i,])
			pred.d2.sim.sd[i]<-sd(pred.d2.sim[i,])
		}
	}else{
		last<-pred$values[[length(pred$values)]]
		if(dim(model$fixef)[2]==2){
			mid<-numeric(0)
		}else{
			mid<-pred$values
			mid[[1]]<-NULL
			mid[[length(mid)]]<-NULL
			mid<-expand.grid(mid)
			mid<-do.call("paste",mid)
		} 
		pred.grid.low<-pred.grid[1:(dim(pred.grid)[1]/length(last)),1:(dim(pred.grid)[2]-1)]
		pred.sim.low<-matrix(0,dim(pred.grid.low)[1],deriv$nsim)
		pred.d1.sim.low<-matrix(0,dim(pred.grid.low)[1],deriv$nsim)
		pred.d2.sim.low<-matrix(0,dim(pred.grid.low)[1],deriv$nsim)
		for(i in 1:length(last)){
			pred.sim.low<-(pred.sim.low+pred.0.sim[1:dim(pred.grid.low)[1]+(i-1)*dim(pred.grid.low)[1],])/length(last)
			pred.d1.sim.low<-(pred.d1.sim.low+pred.d1.0.sim[1:dim(pred.grid.low)[1]+(i-1)*dim(pred.grid.low)[1],])/length(last)
			pred.d2.sim.low<-(pred.d2.sim.low+pred.d2.sim[1:dim(pred.grid.low)[1]+(i-1)*dim(pred.grid.low)[1],])/length(last)
		}
		## mean, averaged across each level of last dimension
		pred.sim.mean<-rowSums(pred.sim.low)/dim(pred.sim.low)[2]
		pred.d1.sim.mean<-rowSums(pred.d1.sim.low)/dim(pred.d1.sim.low)[2]
		pred.d2.sim.mean<-rowSums(pred.d2.sim.low)/dim(pred.d2.sim.low)[2]
		## standard deviation, averaged across each level of last dimension
		pred.sim.sd<-numeric(dim(pred.grid.low)[1])
		pred.d1.sim.sd<-numeric(dim(pred.grid.low)[1])
		pred.d2.sim.sd<-numeric(dim(pred.grid.low)[1])
		for(i in 1:dim(pred.grid.low)[1]){
			pred.sim.sd[i]<-sd(pred.sim.low[i,])
			pred.d1.sim.sd[i]<-sd(pred.d1.sim.low[i,])
			pred.d2.sim.sd[i]<-sd(pred.d2.sim.low[i,])
		}
	}
	## difference measures
	pred.dif<-pred.0-pred.sim.mean
	pred.sim.dif<-pred.0.sim-pred.sim.mean
	pred.d1.dif<-pred.d1.0-pred.d1.sim.mean
	pred.d1.sim.dif<-pred.d1.0.sim-pred.d1.sim.mean
	pred.d2.dif<-pred.d2-pred.d2.sim.mean
	pred.d2.sim.dif<-pred.d2.sim-pred.d2.sim.mean

	## results
	p<-numeric(dim(pred.grid)[1])
	sig<-character(dim(pred.grid)[1])
	d1.p<-numeric(dim(pred.grid)[1])
	d1.sig<-character(dim(pred.grid)[1])
	d2.p<-numeric(dim(pred.grid)[1])
	d2.sig<-character(dim(pred.grid)[1])
	for(i in 1:length(deriv$range)){
		p[i]<-2*abs((length(which(pred.sim.dif[i,]>pred.dif[i]))/deriv$nsim)-0.5)
		sig[i]<-p.stars(p[i])
		d1.p[i]<-2*abs((length(which(pred.d1.sim.dif[i,]>pred.d1.dif[i]))/deriv$nsim)-0.5)
		d1.sig[i]<-p.stars(d1.p[i])
		d2.p[i]<-2*abs((length(which(pred.d2.sim.dif[i,]>pred.d2.dif[i]))/deriv$nsim)-0.5)
		d2.sig[i]<-p.stars(d2.p[i])
	}
	
	## data plots
	## response
	pdf(paste(path,"plot.pred.pdf",sep="/"))
		## set up before adding data
		plot(0,0,pch=30,xlab=pred$names[1],ylab=resp,xlim=c(min(deriv$range),max(deriv$range)),ylim=c(min(c(pred.0,pred.sim.mean-pred.sim.sd)),max(c(pred.0,pred.sim.mean+pred.sim.sd))))
		## mean/sd of lower level fits
		if(length(mid)==0){
			color<-col2rgb("black")
			polygon(c(deriv$range,rev(deriv$range)),c(pred.sim.mean+pred.sim.sd,rev(pred.sim.mean-pred.sim.sd)),col=rgb(color[1],color[2],color[3],75,maxColorValue=255),border=NA)
		}else{
			color<-col2rgb(rainbow(length(mid)))
			for(i in 1:length(mid)) polygon(c(deriv$range,rev(deriv$range)),c(pred.sim.mean[1:length(deriv$range)+(i-1)*length(deriv$range)]+pred.sim.sd[1:length(deriv$range)+(i-1)*length(deriv$range)],rev(pred.sim.mean[1:length(deriv$range)+(i-1)*length(deriv$range)]-pred.sim.sd[1:length(deriv$range)+(i-1)*length(deriv$range)])),col=rgb(color[1,i],color[2,i],color[3,i],75,maxColorValue=255),border=NA)
		}
		## predicted lines
		if(length(last)==0){
			lines(deriv$range,pred.0,lwd=2)
		}else if(length(mid)==0){
			for(i in 1:length(last)) lines(deriv$range,pred.0[1:length(deriv$range)+(i-1)*length(deriv$range)],col=rainbow(length(last))[i],lwd=2)
		}else{
			for(i in 1:length(mid)) for(j in 1:length(last)) lines(deriv$range,pred.0[1:length(deriv$range)+(j-1)*length(deriv$range)+(i-1)*length(deriv$range)*length(last)],col=rainbow(length(mid))[i],lwd=2,lty=j)
		}
	dev.off()
	## d1
	pdf(paste(path,"plot.pred.d1.pdf",sep="/"))
		## set up before adding data
		plot(0,0,pch=30,xlab=pred$names[1],ylab=paste(resp,"'",sep=""),xlim=c(min(deriv$range),max(deriv$range)),ylim=c(min(c(pred.d1.0,pred.d1.sim.mean-pred.d1.sim.sd)),max(c(pred.d1.0,pred.d1.sim.mean+pred.d1.sim.sd))))
		## mean/sd of lower level fits
		if(length(mid)==0){
			color<-col2rgb("black")
			polygon(c(deriv$range,rev(deriv$range)),c(pred.d1.sim.mean+pred.d1.sim.sd,rev(pred.d1.sim.mean-pred.d1.sim.sd)),col=rgb(color[1],color[2],color[3],75,maxColorValue=255),border=NA)
		}else{
			color<-col2rgb(rainbow(length(mid)))
			for(i in 1:length(mid)) polygon(c(deriv$range,rev(deriv$range)),c(pred.d1.sim.mean[1:length(deriv$range)+(i-1)*length(deriv$range)]+pred.d1.sim.sd[1:length(deriv$range)+(i-1)*length(deriv$range)],rev(pred.d1.sim.mean[1:length(deriv$range)+(i-1)*length(deriv$range)]-pred.d1.sim.sd[1:length(deriv$range)+(i-1)*length(deriv$range)])),col=rgb(color[1,i],color[2,i],color[3,i],75,maxColorValue=255),border=NA)
		}
		## predicted lines
		if(length(last)==0){
			lines(deriv$range,pred.d1.0,lwd=2)
		}else if(length(mid)==0){
			for(i in 1:length(last)) lines(deriv$range,pred.d1.0[1:length(deriv$range)+(i-1)*length(deriv$range)],col=rainbow(length(last))[i],lwd=2)
		}else{
			for(i in 1:length(mid)) for(j in 1:length(last)) lines(deriv$range,pred.d1.0[1:length(deriv$range)+(j-1)*length(deriv$range)+(i-1)*length(deriv$range)*length(last)],col=rainbow(length(mid))[i],lwd=2,lty=j)
		}		
	dev.off()	
	## d2
	pdf(paste(path,"plot.pred.d2.pdf",sep="/"))
		## set up before adding data
		plot(0,0,pch=30,xlab=pred$names[1],ylab=paste(resp,"''",sep=""),xlim=c(min(deriv$range),max(deriv$range)),ylim=c(min(c(pred.d2,pred.d2.sim.mean-pred.d2.sim.sd)),max(c(pred.d2,pred.d2.sim.mean+pred.d2.sim.sd))))
		## mean/sd of lower level fits
		if(length(mid)==0){
			color<-col2rgb("black")
			polygon(c(deriv$range,rev(deriv$range)),c(pred.d2.sim.mean+pred.d2.sim.sd,rev(pred.d2.sim.mean-pred.d2.sim.sd)),col=rgb(color[1],color[2],color[3],75,maxColorValue=255),border=NA)
		}else{
			color<-col2rgb(rainbow(length(mid)))
			for(i in 1:length(mid)) polygon(c(deriv$range,rev(deriv$range)),c(pred.d2.sim.mean[1:length(deriv$range)+(i-1)*length(deriv$range)]+pred.d2.sim.sd[1:length(deriv$range)+(i-1)*length(deriv$range)],rev(pred.d2.sim.mean[1:length(deriv$range)+(i-1)*length(deriv$range)]-pred.d2.sim.sd[1:length(deriv$range)+(i-1)*length(deriv$range)])),col=rgb(color[1,i],color[2,i],color[3,i],75,maxColorValue=255),border=NA)
		}
		## predicted lines
		if(length(last)==0){
			lines(deriv$range,pred.d2,lwd=2)
		}else if(length(mid)==0){
			for(i in 1:length(last)) lines(deriv$range,pred.d2[1:length(deriv$range)+(i-1)*length(deriv$range)],col=rainbow(length(last))[i],lwd=2)
		}else{
			for(i in 1:length(mid)) for(j in 1:length(last)) lines(deriv$range,pred.d2[1:length(deriv$range)+(j-1)*length(deriv$range)+(i-1)*length(deriv$range)*length(last)],col=rainbow(length(mid))[i],lwd=2,lty=j)
		}		
	dev.off()	
	
	## p-value plots
	## response
	pdf(paste(path,"plot.p.pdf",sep="/"))
		## set up before adding data
		plot(0,0,pch=30,xlab=pred$names[1],ylab=resp,xlim=c(min(deriv$range),max(deriv$range)),ylim=c(0.9,1),yaxt='n')
		axis(2,c(0.9,0.95,0.99,0.999),c(".1",".05",".01",".001"))
		lines(c(min(deriv$range)-1,max(deriv$range)+1),c(0.9,0.9),lty=2)
		lines(c(min(deriv$range)-1,max(deriv$range)+1),c(0.95,0.95),lty=2)
		lines(c(min(deriv$range)-1,max(deriv$range)+1),c(0.99,0.99),lty=2)
		lines(c(min(deriv$range)-1,max(deriv$range)+1),c(0.999,0.999),lty=2)
		## add data lines
		if(length(last)==0){
			lines(deriv$range,p,lwd=2)
		}else if(length(mid)==0){
			for(i in 1:length(last)) lines(deriv$range,p[1:length(deriv$range)+(i-1)*length(deriv$range)],col=rainbow(length(last))[i],lwd=2)
		}else{
			for(i in 1:length(mid)) for(j in 1:length(last)) lines(deriv$range,p[1:length(deriv$range)+(j-1)*length(deriv$range)+(i-1)*length(deriv$range)*length(last)],col=rainbow(length(mid))[i],lwd=2,lty=j)
		}
	dev.off()
	## d1
	pdf(paste(path,"plot.d1.p.pdf",sep="/"))
		## set up before adding data
		plot(0,0,pch=30,xlab=pred$names[1],ylab=resp,xlim=c(min(deriv$range),max(deriv$range)),ylim=c(0.9,1),yaxt='n')
		axis(2,c(0.9,0.95,0.99,0.999),c(".1",".05",".01",".001"))
		lines(c(min(deriv$range)-1,max(deriv$range)+1),c(0.9,0.9),lty=2)
		lines(c(min(deriv$range)-1,max(deriv$range)+1),c(0.95,0.95),lty=2)
		lines(c(min(deriv$range)-1,max(deriv$range)+1),c(0.99,0.99),lty=2)
		lines(c(min(deriv$range)-1,max(deriv$range)+1),c(0.999,0.999),lty=2)
		## add data lines
		if(length(last)==0){
			lines(deriv$range,d1.p,lwd=2)
		}else if(length(mid)==0){
			for(i in 1:length(last)) lines(deriv$range,d1.p[1:length(deriv$range)+(i-1)*length(deriv$range)],col=rainbow(length(last))[i],lwd=2)
		}else{
			for(i in 1:length(mid)) for(j in 1:length(last)) lines(deriv$range,d1.p[1:length(deriv$range)+(j-1)*length(deriv$range)+(i-1)*length(deriv$range)*length(last)],col=rainbow(length(mid))[i],lwd=2,lty=j)
		}
	dev.off()
	## d2
	pdf(paste(path,"plot.d2.p.pdf",sep="/"))
		## set up before adding data
		plot(0,0,pch=30,xlab=pred$names[1],ylab=resp,xlim=c(min(deriv$range),max(deriv$range)),ylim=c(0.9,1),yaxt='n')
		axis(2,c(0.9,0.95,0.99,0.999),c(".1",".05",".01",".001"))
		lines(c(min(deriv$range)-1,max(deriv$range)+1),c(0.9,0.9),lty=2)
		lines(c(min(deriv$range)-1,max(deriv$range)+1),c(0.95,0.95),lty=2)
		lines(c(min(deriv$range)-1,max(deriv$range)+1),c(0.99,0.99),lty=2)
		lines(c(min(deriv$range)-1,max(deriv$range)+1),c(0.999,0.999),lty=2)
		## add data lines
		if(length(last)==0){
			lines(deriv$range,d2.p,lwd=2)
		}else if(length(mid)==0){
			for(i in 1:length(last)) lines(deriv$range,d2.p[1:length(deriv$range)+(i-1)*length(deriv$range)],col=rainbow(length(last))[i],lwd=2)
		}else{
			for(i in 1:length(mid)) for(j in 1:length(last)) lines(deriv$range,d2.p[1:length(deriv$range)+(j-1)*length(deriv$range)+(i-1)*length(deriv$range)*length(last)],col=rainbow(length(mid))[i],lwd=2,lty=j)
		}
	dev.off()

	## legends for data plots (only necessary if length(last)>0)
	if(length(last)>0){
		pdf(paste(path,"plot.legend.pdf",sep="/"))
			## set up before adding data
			plot(c(0,1),c(0,1),pch=30,axes=FALSE,ann=FALSE)
			## add legend
			if(length(mid)==0){
				y<-seq(1,0,length.out=(length(last)+2))[2:(length(last)+1)]
				for(i in 1:length(last)){
					lines(c(1/5,2/5),c(y[i],y[i]),col=rainbow(length(last))[i],lwd=2)
					text(3/5,y[i],last[i])
				}
			}else{
				y<-seq(1,0,length.out=(length(last)+length(mid)+2))[2:(length(last)+length(mid)+1)]
				for(i in 1:length(mid)){
					lines(c(1/5,2/5),c(y[i],y[i]),col=rainbow(length(mid))[i],lwd=2)
					text(3/5,y[i],mid[i])
				}
				for(i in 1:length(last)){
					lines(c(1/5,2/5),c(y[i+length(mid)],y[i+length(mid)]),lwd=2,lty=i)
					text(3/5,y[i+length(mid)],last[i])
				}
			}
		dev.off()
	}
	
	## return results
	data.frame(pred.grid,pred.dif,p,sig,pred.d1.dif,d1.p,d1.sig,pred.d2.dif,d2.p,d2.sig)
}

## ANALYSIS

## data
load(paste(setup$path,"data.R",sep="/"))

## tract.data
load(paste(setup$path,"tract.data.R",sep="/"))

## tract names
tract.names<-c("skel_Middle_cerebellar_peduncle","skel_Pontine_crossing_tract","skel_Genu_of_corpus_callosum","skel_Body_of_corpus_callosum","skel_Splenium_of_corpus_callosum","skel_Fornix_column_body","skel_Corticospinal_tract_R","skel_Corticospinal_tract_L","skel_Medial_lemniscus_R","skel_Medial_lemniscus_L","skel_Inferior_cerebellar_peduncle_R","skel_Inferior_cerebellar_peduncle_L","skel_Superior_cerebellar_peduncle_R","skel_Superior_cerebellar_peduncle_L","skel_Cerebral_peduncle_R","skel_Cerebral_peduncle_L","skel_Anterior_limb_of_internal_capsule_R","skel_Anterior_limb_of_internal_capsule_L","skel_Posterior_limb_of_internal_capsule_R","skel_Posterior_limb_of_internal_capsule_L","skel_Retrolenticular_part_of_internal_capsule_R","skel_Retrolenticular_part_of_internal_capsule_L","skel_Anterior_corona_radiata_R","skel_Anterior_corona_radiata_L","skel_Superior_corona_radiata_R","skel_Superior_corona_radiata_L","skel_Posterior_corona_radiata_R","skel_Posterior_corona_radiata_L","skel_Posterior_thalamic_radiation_R","skel_Posterior_thalamic_radiation_L","skel_Sagittal_stratum_R","skel_Sagittal_stratum_L","skel_External_capsule_R","skel_External_capsule_L","skel_Cingulum_cingulate_gyrus_R","skel_Cingulum_cingulate_gyrus_L","skel_Cingulum_hippocampus_R","skel_Cingulum_hippocampus_L","skel_Fornix_cres-Stria_terminalis_R","skel_Fornix_cres-Stria_terminalis_L","skel_Superior_longitudinal_fasciculus_R","skel_Superior_longitudinal_fasciculus_L","skel_Superior_fronto-occipital_fasciculus_R","skel_Superior_fronto-occipital_fasciculus_L","skel_Uncinate_fasciculus_R","skel_Uncinate_fasciculus_L","skel_Tapetum_R.txt","skel_Tapetum_L.txt")

## models
models<-modelsetup()
sink(paste(setup$path,"models.txt",sep="/"))
print(models)
sink()

## logfile
log.txt=paste(setup$path,"log.txt",sep="/")
sink(log.txt)
cat("LONGITUDINAL ANALYSIS\t",date(),"\n\nNOTE: there are occasionally warnings from the graphs concerning pch=30. you can ignore these.\n\n")

## temporary
data$fa<-0
cat("vgs.mRT cutpoints:",quantile(data$vgs.mRT,1/3,na.rm=TRUE),quantile(data$vgs.mRT,2/3,na.rm=TRUE),"\n")
cat("vgs.cv cutpoints:",quantile(data$vgs.cv,1/3,na.rm=TRUE),quantile(data$vgs.cv,2/3,na.rm=TRUE),"\n")
cat("anti.percErr cutpoints:",quantile(data$anti.percErr,1/3,na.rm=TRUE),quantile(data$anti.percErr,2/3,na.rm=TRUE),"\n\n")
data<-cont.to.fact("vgs.mRT",data,cut=c(1/3,2/3))
data<-cont.to.fact("vgs.cv",data,cut=c(1/3,2/3))
data<-cont.to.fact("anti.percErr",data,cut=c(1/3,2/3))

for(t in 1:length(tract.names)){
	## temporary code
	data.<-data
	data.$fa<-tract.data[,t]
	data.na<-unique(c(which(is.na(data.$fa)),which(is.na(data.$age)),which(is.na(data.$id))))
	if(length(data.na)>0) data.<-data.[-data.na,]
	data.$wt<-outlier.wt(data.$fa,data.$age,data.$id)[,3]
	
	cat(tract.names[t],"\n")
	path<-paste(setup$path,t,sep="/")
	dir.create(path)
	## set up results directories and files
	anova.txt=paste(path,"anova.txt",sep="/")
	anova.path=paste(path,"anova",sep="/")
	dir.create(anova.path)
	deriv.path=paste(path,"deriv",sep="/")
	dir.create(deriv.path)
	## set up structures for estimation and results
	models$model<-list()
	models$ll.test<-list()
	models$ll.test.W<-list()
	models$deriv<-list()
	## null model
	models$model[[1]]<-lmer.est()
	## loop through all models
	for(m in 2:dim(models$X)[1]){
		cat("\n",models$X.var[m,],"\n")
		models$model[[m]]<-lmer.est(fixef.set(models$X[m,]))
		models$ll.test[[m]]<-list()
		for(c in 1:length(models$con[[m]])){
			cat(models$model[[m]]$formula,"vs",models$model[[models$con[[m]][c]]]$formula,"\n")
			comb.na<-numeric(0)
			if(length(models$model[[m]]$na)==0 & length(models$model[[models$con[[m]][c]]]$na)==0){
				models$ll.test[[m]][[c]]<-anova(models$model[[m]]$fit,(models$model[[models$con[[m]][c]]]$fit))
			}else{
				comb.na<-unique(union(models$model[[m]]$na,models$model[[models$con[[m]][c]]]$na))
				models$ll.test[[m]][[c]]<-anova(lmer.est(models$model[[m]]$fixef,data=data.[-comb.na,])$fit,(lmer.est(models$model[[models$con[[m]][c]]]$fixef,data=data.[-comb.na,])$fit))
			}
			sink(paste(anova.path,"/anova_",m,"-",models$con[[m]][c],".txt",sep=""))
			print(models$ll.test[[m]][[c]])
			cat("NA:",comb.na)
			sink()
			sink(anova.txt,append=TRUE)
			cat(models$model[[m]]$formula,"vs",models$model[[models$con[[m]][c]]]$formula,"\tp =",models$ll.test[[m]][[c]][[7]][2],"\n")
			sink()
			if(c==1 & deriv.check(models$X[m,],models$ll.test[[m]][[c]][[7]][2])){
				cat("derivatives analysis\n")
				cat("START:",date(),"\n")
				dir.create(paste(deriv.path,m,sep="/"))
				if(length(comb.na)==0){
					models$deriv[[m]]<-deriv.est(models$model[[m]],models$model[[models$con[[m]][c]]],path=paste(deriv.path,m,sep="/"))
				}else{
					models$deriv[[m]]<-deriv.est(lmer.est(models$model[[m]]$fixef,data=data.[-comb.na,]),lmer.est(models$model[[models$con[[m]][c]]]$fixef,data=data.[-comb.na,]),path=paste(deriv.path,m,sep="/"),data=data.[-comb.na,])
				}
				cat("END:",date(),"\n")
			}
		}
		if(length(models$con.W[[m]])>0){
			models$ll.test.W[[m]]<-list()
			for(c in 1:length(models$con.W[[m]])){
				cat(models$model[[m]]$formula,"vs",models$model[[models$con.W[[m]][c]]]$formula,"\n")
				if(length(models$model[[m]]$na)==0 & length(models$model[[models$con.W[[m]][c]]]$na)==0){
					models$ll.test.W[[m]][[c]]<-anova(models$model[[m]]$fit,(models$model[[models$con.W[[m]][c]]]$fit))
				}else{
					comb.na<-unique(union(models$model[[m]]$na,models$model[[models$con.W[[m]][c]]]$na))
					models$ll.test.W[[m]][[c]]<-anova(lmer.est(models$model[[m]]$fixef,data=data.[-comb.na,])$fit,(lmer.est(models$model[[models$con.W[[m]][c]]]$fixef,data=data.[-comb.na,])$fit))
				}
				sink(paste(anova.path,"/anova_",m,"-",models$con.W[[m]][c],".txt",sep=""))
				print(models$ll.test.W[[m]][[c]])
				cat("NA:",comb.na)
				sink()
				sink(anova.txt,append=TRUE)
				cat(models$model[[m]]$formula,"vs",models$model[[models$con.W[[m]][c]]]$formula,"\tp =",models$ll.test.W[[m]][[c]][[7]][2],"\n")
				sink()
			}
		}
	}
	save(models,file=paste(path,paste("models",t,".R",sep=""),sep="/"))
}
sink()
