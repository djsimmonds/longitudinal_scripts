deriv.est<-function(m1=models$model[[m]],m0=models$model[[models$con[[m]][c]]],path=path..,deriv=setup$deriv,mixed=setup$mixed,data=data.){
## derivatives analysis
## details:
	## need to fill in

	## simulate from lower level model
	cat(date(),"\t\t\t\t\tsimulating from null, n =",deriv$nsim,"\n")
	sim<-simulate(m0$fit,deriv$nsim)

	## create new data structure for predictions of bootstrapped models
	pred<-list(names=m1$fixef[1,1],values=list(deriv$range-8))
	if(dim(m1$fixef)[2]>1){
		for(i in 2:dim(m1$fixef)[2]){
			col<-which(names(data)==m1$fixef[1,i])
			pred$names[i]<-m1$fixef[1,i]
			pred$values[[i]]<-unique(factor(data[,col]))
		}
	}
	names(pred$values)<-pred$names
	pred.grid<-expand.grid(pred$values)

	## expanding grid for short interval calculations
	pred.grid.ind<-rep(1:dim(pred.grid)[1],each=2)
	pred.grid.int<-rbind(pred.grid,pred.grid)
	for(i in 1:length(pred.grid.ind)) pred.grid.int[i,]<-pred.grid[pred.grid.ind[i],]
	pred.grid.int[which(1:length(pred.grid.ind)%%2==0),1]<-pred.grid.int[which(1:length(pred.grid.ind)%%2==0),1]+deriv$int
	for(i in 1:length(pred$names)) assign(pred$names[i],pred.grid.int[,i],envir=.GlobalEnv) ## necessary for predictSE.mer, has to do with the use of functions like ns() within formulas within functions, without having the variable in the function outside the ns() structure --> i don't fully understand, but somehow by declaring the variables in the main environment, it gets around the problem

	## predicting values for true and simulated fits
	cat(date(),"\t\t\t\t\tfitting simulated data:")
	if(mixed==TRUE){
		pred.<-predictSE.mer(m1$fit,pred.grid.int,se.fit=FALSE)
		pred.sim<-matrix(NA,length(pred.),deriv$nsim)
		for(i in 1:deriv$nsim){
			cat(i,"")
			pred.sim[,i]<-predictSE.mer(refit(m1$fit,sim[,i]),pred.grid.int,se.fit=FALSE)
		}
	}else{
		pred.<-predict(m1$fit,pred.grid.int)
		pred.sim<-matrix(NA,length(pred.),deriv$nsim)
		for(i in 1:deriv$nsim){
			cat(i,"")
			data[,1]<-sim[,1]
			m.sim<-update(m1$fit,.~.)
			pred.sim[,i]<-predict(m.sim,pred.grid.int)
		}
	}
	cat("\n",date(),"\t\t\t\t\tfitting complete\n")

	## first derivative
	pred.d1<-(pred.[which(1:length(pred.grid.ind)%%2==0)]-pred.[which(1:length(pred.grid.ind)%%2==1)])/deriv$int
	pred.d1.sim<-(pred.sim[which(1:length(pred.grid.ind)%%2==0),]-pred.sim[which(1:length(pred.grid.ind)%%2==1),])/deriv$int
	## reorganizing initial predictions to remove small interval
	pred.<-pred.[which(1:length(pred.grid.ind)%%2==1)]
	pred.sim<-pred.sim[which(1:length(pred.grid.ind)%%2==1),]

	## mean and sd calculations
	## if L>1, average across last dimension for plotting and distance calculations
	if(dim(m1$fixef)[2]==1){
		last<-numeric(0)
		mid<-numeric(0)
		## mean of simulations
		pred.sim.mean<-rowSums(pred.sim)/dim(pred.sim)[2]
		pred.d1.sim.mean<-rowSums(pred.d1.sim)/dim(pred.d1.sim)[2]
		## sd of simulations
		pred.sim.sd<-numeric(length(pred.sim.mean))
		pred.d1.sim.sd<-numeric(length(pred.d1.sim.mean))
		for(i in 1:dim(pred.grid)[1]){
			pred.sim.sd[i]<-sd(pred.sim[i,])
			pred.d1.sim.sd[i]<-sd(pred.d1.sim[i,])
		}
	}else{
		last<-pred$values[[length(pred$values)]]
		if(dim(m1$fixef)[2]==2){
			mid<-numeric(0)
		}else{
			mid<-pred$values
			mid[[1]]<-NULL
			mid[[length(mid)]]<-NULL
			mid<-expand.grid(mid)
			mid<-do.call("paste",c(mid,sep="_"))
		} 
		pred.grid.low<-pred.grid[1:(dim(pred.grid)[1]/length(last)),1:(dim(pred.grid)[2]-1)]
		dim.low<-c(min(length(pred.grid.low),dim(pred.grid.low)[1]),max(1,dim(pred.grid.low)[2]))
		pred.sim.low<-matrix(0,dim.low[1],deriv$nsim)
		pred.d1.sim.low<-matrix(0,dim.low[1],deriv$nsim)
		for(i in 1:length(last)){
			pred.sim.low<-(pred.sim.low+pred.sim[1:dim.low[1]+(i-1)*dim.low[1],])/length(last)
			pred.d1.sim.low<-(pred.d1.sim.low+pred.d1.sim[1:dim.low[1]+(i-1)*dim.low[1],])/length(last)
		}
		## mean, averaged across each level of last dimension
		pred.sim.mean<-rowSums(pred.sim.low)/dim.low[2]
		pred.d1.sim.mean<-rowSums(pred.d1.sim.low)/dim.low[2]
		## standard deviation, averaged across each level of last dimension
		pred.sim.sd<-numeric(dim.low[1])
		pred.d1.sim.sd<-numeric(dim.low[1])
		for(i in 1:dim.low[1]){
			pred.sim.sd[i]<-sd(pred.sim.low[i,])
			pred.d1.sim.sd[i]<-sd(pred.d1.sim.low[i,])
		}
	}

	## difference measures
	pred.dif<-pred.-pred.sim.mean
	mean.dif<-numeric(length(deriv$range))
	pred.sim.dif<-pred.sim-pred.sim.mean
	mean.sim.dif<-matrix(0,length(deriv$range),deriv$nsim)
	pred.d1.dif<-pred.d1-pred.d1.sim.mean
	mean.d1.dif<-numeric(length(deriv$range))
	pred.d1.sim.dif<-pred.d1.sim-pred.d1.sim.mean
	mean.d1.sim.dif<-matrix(0,length(deriv$range),deriv$nsim)
	for(i in 1:length(deriv$range)){
		mean.dif[i]<-sqrt(sum(pred.dif[which((1:length(pred.dif))%%length(deriv$range)==i)]^2))
		mean.d1.dif[i]<-sqrt(sum(pred.d1.dif[which((1:length(pred.d1.dif))%%length(deriv$range)==i)]^2))
		for(j in 1:length(deriv$nsim)){
			mean.sim.dif[i,j]<-sqrt(sum(pred.sim.dif[which((1:length(pred.dif))%%length(deriv$range)==i),j]^2))
			mean.d1.sim.dif[i,j]<-sqrt(sum(pred.d1.sim.dif[which((1:length(pred.d1.dif))%%length(deriv$range)==i),j]^2))
		}
	}

	## p-values
	p<-numeric(length(pred.dif))
	p.abs<-numeric(length(pred.dif))
	p.d1<-numeric(length(pred.d1.dif))
	p.d1.abs<-numeric(length(pred.d1.dif))
	for(i in 1:length(pred.dif)){
		p[i]<-length(which(pred.sim.dif[i,]>pred.dif[i]))/deriv$nsim
		p.abs[i]<-2*abs(p[i]-0.5)
		p.d1[i]<-length(which(pred.d1.sim.dif[i,]>pred.d1.dif[i]))/deriv$nsim
		p.d1.abs[i]<-2*abs(p.d1[i]-0.5)
	}
	p.dif<-numeric(length(deriv$range))
	p.dif.abs<-numeric(length(deriv$range))
	p.d1.dif<-numeric(length(deriv$range))
	p.d1.dif.abs<-numeric(length(deriv$range))
	for(i in 1:length(deriv$range)){
		p.dif[i]<-length(which(mean.sim.dif[i,]>mean.dif[i]))/deriv$nsim
		p.dif.abs[i]<-2*abs(p.dif[i]-0.5)
		p.d1.dif[i]<-length(which(mean.d1.sim.dif[i,]>mean.d1.dif[i]))/deriv$nsim
		p.d1.dif.abs[i]<-2*abs(p.d1.dif[i]-0.5)
	}

	## results
	sink(paste(path,"grid.txt",sep="/"))
	print(pred.grid)
	sink()
	sink(paste(path,"pred.txt",sep="/"))
	cat(pred.)
	sink()
	sink(paste(path,"pred.d1.txt",sep="/"))
	cat(pred.d1)
	sink()
	sink(paste(path,"sim.mean.txt",sep="/"))
	cat(pred.sim.mean)
	sink()
	sink(paste(path,"sim.sd.txt",sep="/"))
	cat(pred.sim.sd)
	sink()
	sink(paste(path,"sim.d1.mean.txt",sep="/"))
	cat(pred.d1.sim.mean)
	sink()
	sink(paste(path,"sim.d1.sd.txt",sep="/"))
	cat(pred.d1.sim.sd)
	sink()
	sink(paste(path,"p.txt",sep="/"))
	cat(p)
	sink()
	sink(paste(path,"p.abs.txt",sep="/"))
	cat(p.abs)
	sink()
	sink(paste(path,"p.dif.txt",sep="/"))
	cat(p.dif)
	sink()
	sink(paste(path,"p.dif.abs.txt",sep="/"))
	cat(p.dif.abs)
	sink()
	sink(paste(path,"p.d1.txt",sep="/"))
	cat(p.d1)
	sink()
	sink(paste(path,"p.d1.abs.txt",sep="/"))
	cat(p.d1.abs)
	sink()
	sink(paste(path,"p.d1.dif.txt",sep="/"))
	cat(p.d1.dif)
	sink()
	sink(paste(path,"p.d1.dif.abs.txt",sep="/"))
	cat(p.d1.dif.abs)
	sink()

	return(NULL)
}
