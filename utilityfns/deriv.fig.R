deriv.fig<-function(analysis.path,model.number,type,ynames=NULL,figs=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),main.offset=8){
## path: within model derivatives directory
	path<-paste(analysis.path,model.number,"deriv",sep="/")
## type: 1 to 5
## figs: logical of whether to create figures
	## predicted values
		## 1) predicted values
		## 2) pvals (from normal distribution)
		## 3) pvals (from monte carlo)
	## predicted derivatives
		## 4) predicted values
		## 5) pvals (from normal distribution)
		## 6) pvals (from monte carlo)

	## create directories for images and load tables
	tables<-"pred.grid"
	if(figs[1]){
		dir.create(paste(path,"pred",sep="/"))
		tables<-c(tables,c("pred","sim.pred.mean","sim.pred.sd"))
		if(type>1) tables<-c(tables,c("pred.dif","sim.pred.dif.mean","sim.pred.dif.sd"))
	}
	if(figs[2]){
		dir.create(paste(path,"pred.p",sep="/"))
		tables<-c(tables,"pred.p")
		if(type>1) tables<-c(tables,"pred.dif.p")
	}
	if(figs[3]){
		dir.create(paste(path,"pred.pboot",sep="/"))
		tables<-c(tables,"pred.pboot")
		if(type>1) tables<-c(tables,"pred.dif.pboot")
	}
	if(figs[4]){
		dir.create(paste(path,"pred.d",sep="/"))
		tables<-c(tables,c("pred.d","sim.pred.d.mean","sim.pred.d.sd"))
		if(type>1) tables<-c(tables,c("pred.d.dif","sim.pred.d.dif.mean","sim.pred.d.dif.sd"))
	}
	if(figs[5]){
		dir.create(paste(path,"pred.d.p",sep="/"))
		tables<-c(tables,"pred.d.p")
		if(type>1) tables<-c(tables,"pred.d.dif.p")
	}
	if(figs[6]){
		dir.create(paste(path,"pred.d.pboot",sep="/"))
		tables<-c(tables,"pred.d.pboot")
		if(type>1) tables<-c(tables,"pred.d.dif.pboot")
	}
	## load tables
	path.tables<-paste(path,"tables",sep="/")
	for(t in 1:length(tables)){
		filename<-paste(path.tables,tables[t],sep="/")
		if(t==1) assign(tables[t],read.table(filename)) else assign(tables[t],read.table(filename,header=TRUE,row.names=1))
		if(!is.matrix(get(tables[t]))) assign(tables[t],as.matrix(get(tables[t])))
	}
	## range for main variables
	range<-as.numeric(unique(pred.grid[,1]))


	## different plots for different types

	
	#######
	## 1 ##
	#######

	if(type==1){

		## for each y
		for(i in 1:dim(pred)[2]){

			## predicted fits
			if(figs[1]){
				## for setting y limits on graph
				pred.<-as.numeric(pred[,i])
				sim.pred.mean.<-as.numeric(sim.pred.mean[,i])
				sim.pred.sd.<-as.numeric(sim.pred.sd[,i])

				pdf(paste(path,"pred",paste(i,"pdf",sep="."),sep="/"))
					## prepare plot
					yL<-min(c(min(pred.),min(sim.pred.mean.-sim.pred.sd.)))
					yH<-max(c(max(pred.),max(sim.pred.mean.+sim.pred.sd.)))
					plot(0,0,pch=30,main=ynames[i],xlim=range(range)+main.offset,xlab="age(y)",ylim=c(yL,yH),ylab="FA maturation index")
					## add null distribution
					polygon(c(range+main.offset,rev(range+main.offset)),c(sim.pred.mean.+sim.pred.sd.,rev(sim.pred.mean.-sim.pred.sd.)),col=rgb(0,0,0,75,maxColorValue=255),border=NA)
					## add predicted line
					lines(range+main.offset,pred.,lwd=2)
				dev.off()
			}
			
			## pvals (normal distribution)
			if(figs[2]){
				pdf(paste(path,"pred.p",paste(i,"pdf",sep="."),sep="/"))
					plot(0,0,pch=30,main=ynames[i],xlim=range(range)+main.offset,xlab="age(y)",ylim=c(0,0.1),ylab="pval")
					lines(range(range)+main.offset,c(0.05,0.05),lty=3)
					lines(range(range)+main.offset,c(0.01,0.01),lty=3)
					lines(range(range)+main.offset,c(0.001,0.001),lty=3)
					lines(range+main.offset,pred.p[,i],lwd=2)
				dev.off()
			}
			
			## pvals (monte carlo)
			if(figs[3]){
				pdf(paste(path,"pred.pboot",paste(i,"pdf",sep="."),sep="/"))
					plot(0,0,pch=30,main=ynames[i],xlim=range(range)+main.offset,xlab="age(y)",ylim=c(0,0.1),ylab="pval")
					lines(range(range)+main.offset,c(0.05,0.05),lty=3)
					lines(range(range)+main.offset,c(0.01,0.01),lty=3)
					lines(range(range)+main.offset,c(0.001,0.001),lty=3)
					lines(range+main.offset,pred.pboot[,i],lwd=2)
				dev.off()
			}
			
			## predicted fit derivatives
			if(figs[4]){
				## for setting y limits on graph
				pred.d.<-as.numeric(pred.d[,i])
				sim.pred.d.mean.<-as.numeric(sim.pred.d.mean[,i])
				sim.pred.d.sd.<-as.numeric(sim.pred.d.sd[,i])

				## prediction - derivatives
				pdf(paste(path,"pred.d",paste(i,"pdf",sep="."),sep="/"))
					## prepare plot
					yL<-min(c(min(pred.d.),min(sim.pred.d.mean.-sim.pred.d.sd.)))
					yH<-max(c(max(pred.d.),max(sim.pred.d.mean.+sim.pred.d.sd.)))
					plot(0,0,pch=30,main=ynames[i],xlim=range(range)+main.offset,xlab="age(y)",ylim=c(yL,yH),ylab="derivative of maturation index")
					## add null distribution
					polygon(c(range+main.offset,rev(range+main.offset)),c(sim.pred.d.mean.+sim.pred.d.sd.,rev(sim.pred.d.mean.-sim.pred.d.sd.)),col=rgb(0,0,0,75,maxColorValue=255),border=NA)
					## add predicted line
					lines(range+main.offset,pred.d.,lwd=2)
				dev.off()
			}

			## derivative pvals (normal distribution)
			if(figs[5]){
				pdf(paste(path,"pred.d.p",paste(i,"pdf",sep="."),sep="/"))
					plot(0,0,pch=30,main=ynames[i],xlim=range(range)+main.offset,xlab="age(y)",ylim=c(0,0.1),ylab="pval")
					lines(range(range)+main.offset,c(0.05,0.05),lty=3)
					lines(range(range)+main.offset,c(0.01,0.01),lty=3)
					lines(range(range)+main.offset,c(0.001,0.001),lty=3)
					lines(range+main.offset,pred.d.p[,i],lwd=2)
				dev.off()
			}

			## derivative pvals (monte carlo)
			if(figs[6]){
				pdf(paste(path,"pred.d.pboot",paste(i,"pdf",sep="."),sep="/"))
					plot(0,0,pch=30,main=ynames[i],xlim=range(range)+main.offset,xlab="age(y)",ylim=c(0,0.1),ylab="pval")
					lines(range(range)+main.offset,c(0.05,0.05),lty=3)
					lines(range(range)+main.offset,c(0.01,0.01),lty=3)
					lines(range(range)+main.offset,c(0.001,0.001),lty=3)
					lines(range+main.offset,pred.d.pboot[,i],lwd=2)
				dev.off()
			}
		}

		## all y in 1 graph
		if(dim(pred)[2]>1){

			for(t in 2:length(tables)) assign(tables[t],as.matrix(get(tables[t])))
			
			if(figs[1]){
				pdf(paste(path,"pred.pdf",sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(min(pred),max(pred),length.out=20)
					image(matrix(range.,1,length(range.)),col=topo.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,3))
					title(ylab="maturation index")
					image(range+main.offset,1:dim(pred)[2],pred,col=topo.colors(length(range.)),xlim=range(range+main.offset),xlab="age(y)",yaxt="n",ylab="Y")
					axis(2,1:length(ynames),ynames,las=1)	
				dev.off()
			}

			if(figs[2]){
				pdf(paste(path,"pred.p.pdf",sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(0,0.1,length.out=20)
					image(matrix(range.,1,length(range.)),col=heat.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,2))
					title(ylab="p")
					image(range+main.offset,1:dim(pred)[2],pred.p,col=heat.colors(length(range.)),zlim=c(0,0.1),xlim=range(range+main.offset),xlab="age(y)",yaxt="n",ylab="Y")
					axis(2,1:length(ynames),ynames,las=1)	
				dev.off()
			}

			if(figs[3]){
				pdf(paste(path,"pred.pboot.pdf",sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(0,0.1,length.out=20)
					image(matrix(range.,1,length(range.)),col=heat.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,2))
					title(ylab="p")
					image(range+main.offset,1:dim(pred)[2],pred.pboot,col=heat.colors(length(range.)),zlim=c(0,0.1),xlim=range(range+main.offset),xlab="age(y)",yaxt="n",ylab="Y")
					axis(2,1:length(ynames),ynames,las=1)	
				dev.off()
			}

			if(figs[4]){
				pdf(paste(path,"pred.d.pdf",sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(min(pred.d),max(pred.d),length.out=20)
					image(matrix(range.,1,length(range.)),col=topo.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,3))
					title(ylab="maturation index derivative")
					image(range+main.offset,1:dim(pred)[2],pred.d,col=topo.colors(length(range.)),xlim=range(range+main.offset),xlab="age(y)",yaxt="n",ylab="Y")
					axis(2,1:length(ynames),ynames,las=1)	
				dev.off()
			}

			if(figs[5]){
				pdf(paste(path,"pred.d.p.pdf",sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(0,0.1,length.out=20)
					image(matrix(range.,1,length(range.)),col=heat.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,2))
					title(ylab="p")
					image(range+main.offset,1:dim(pred)[2],pred.d.p,col=heat.colors(length(range.)),zlim=c(0,0.1),xlim=range(range+main.offset),xlab="age(y)",yaxt="n",ylab="Y")
					axis(2,1:length(ynames),ynames,las=1)	
				dev.off()
			}

			if(figs[6]){
				pdf(paste(path,"pred.d.pboot.pdf",sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(0,0.1,length.out=20)
					image(matrix(range.,1,length(range.)),col=heat.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,2))
					title(ylab="p")
					image(range+main.offset,1:dim(pred)[2],pred.d.pboot,col=heat.colors(length(range.)),zlim=c(0,0.1),xlim=range(range+main.offset),xlab="age(y)",yaxt="n",ylab="Y")
					axis(2,1:length(ynames),ynames,las=1)
				dev.off()
			}
		}
	}


	#######
	## 2 ##
	#######

	if(type==2){

		num.levels<-unique(pred.grid[,2])

		## convert each y to a matrix format
		pred.list<-list()
		pred.p.list<-list()
		pred.pboot.list<-list()
		pred.d.list<-list()
		pred.d.p.list<-list()
		pred.d.pboot.list<-list()
		for(i in 1:dim(pred)[2]){
			pred.list[[i]]<-sapply(1:length(num.levels), function(j) pred[which(pred.grid[,2]==num.levels[j]),i])
			pred.p.list[[i]]<-sapply(1:length(num.levels), function(j) pred.p[which(pred.grid[,2]==num.levels[j]),i])
			pred.pboot.list[[i]]<-sapply(1:length(num.levels), function(j) pred.pboot[which(pred.grid[,2]==num.levels[j]),i])
			pred.d.list[[i]]<-sapply(1:length(num.levels), function(j) pred.d[which(pred.grid[,2]==num.levels[j]),i])
			pred.d.p.list[[i]]<-sapply(1:length(num.levels), function(j) pred.d.p[which(pred.grid[,2]==num.levels[j]),i])
			pred.d.pboot.list[[i]]<-sapply(1:length(num.levels), function(j) pred.d.pboot[which(pred.grid[,2]==num.levels[j]),i])
		}

		## for each y
		for(i in 1:dim(pred)[2]){

			## predicted fits
			if(figs[1]){

				## for setting y limits on graph
				pred.<-as.numeric(pred[,i])
				sim.pred.mean.<-as.numeric(sim.pred.mean[,i])
				sim.pred.sd.<-as.numeric(sim.pred.sd[,i])

				pdf(paste(path,"pred",paste(i,"pdf",sep="."),sep="/"))
					yL<-min(c(min(pred.),min(sim.pred.mean.-sim.pred.sd.)))
					yH<-max(c(max(pred.),max(sim.pred.mean.+sim.pred.sd.)))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(yL,yH,length.out=20)
					image(matrix(range.,1,length(range.)),col=topo.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,3))
					title(ylab="maturation index")
					image(range+main.offset,num.levels,pred.list[[i]],col=topo.colors(length(range.)),xlim=range(range+main.offset),xlab="age(y)",ylim=range(num.levels),ylab=names(pred.grid[2]))
					title(main=ynames[i])
				dev.off()
			}

			## pvals (normal distribution)
			if(figs[2]){
				pdf(paste(path,"pred.p",paste(i,"pdf",sep="."),sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(0,0.1,length.out=20)
					image(matrix(range.,1,length(range.)),col=heat.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,2))
					title(ylab="p")
					image(range+main.offset,num.levels,pred.p.list[[i]],col=heat.colors(length(range.)),zlim=c(0,0.1),xlim=range(range+main.offset),xlab="age(y)",ylim=range(num.levels),ylab=names(pred.grid[2]))
					title(main=ynames[i])
				dev.off()
			}
			
			## pvals (monte carlo)
			if(figs[3]){
				pdf(paste(path,"pred.pboot",paste(i,"pdf",sep="."),sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(0,0.1,length.out=20)
					image(matrix(range.,1,length(range.)),col=heat.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,2))
					title(ylab="p")
					image(range+main.offset,num.levels,pred.pboot.list[[i]],col=heat.colors(length(range.)),zlim=c(0,0.1),xlim=range(range+main.offset),xlab="age(y)",ylim=range(num.levels),ylab=names(pred.grid[2]))
					title(main=ynames[i])
				dev.off()
			}
			
			## predicted fit derivatives
			if(figs[4]){
				## for setting y limits on graph
				pred.d.<-as.numeric(pred.d[,i])
				sim.pred.d.mean.<-as.numeric(sim.pred.d.mean[,i])
				sim.pred.d.sd.<-as.numeric(sim.pred.d.sd[,i])

				## prediction - derivatives
				pdf(paste(path,"pred.d",paste(i,"pdf",sep="."),sep="/"))
					## prepare plot
					yL<-min(c(min(pred.d.),min(sim.pred.d.mean.-sim.pred.d.sd.)))
					yH<-max(c(max(pred.d.),max(sim.pred.d.mean.+sim.pred.d.sd.)))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(yL,yH,length.out=20)
					image(matrix(range.,1,length(range.)),col=topo.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,3))
					title(ylab="maturation index")
					image(range+main.offset,num.levels,pred.d.list[[i]],col=topo.colors(length(range.)),xlim=range(range+main.offset),xlab="age(y)",ylim=range(num.levels),ylab=names(pred.grid[2]))
					title(main=ynames[i])
				dev.off()
			}

			## derivative pvals (normal distribution)
			if(figs[5]){
				pdf(paste(path,"pred.d.p",paste(i,"pdf",sep="."),sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(0,0.1,length.out=20)
					image(matrix(range.,1,length(range.)),col=heat.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,2))
					title(ylab="p")
					image(range+main.offset,num.levels,pred.d.p.list[[i]],col=heat.colors(length(range.)),zlim=c(0,0.1),xlim=range(range+main.offset),xlab="age(y)",ylim=range(num.levels),ylab=names(pred.grid[2]))
					title(main=ynames[i])
				dev.off()
			}
			
			## derivative pvals (monte carlo)
			if(figs[6]){
				pdf(paste(path,"pred.d.pboot",paste(i,"pdf",sep="."),sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(0,0.1,length.out=20)
					image(matrix(range.,1,length(range.)),col=heat.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,2))
					title(ylab="p")
					image(range+main.offset,num.levels,pred.d.pboot.list[[i]],col=heat.colors(length(range.)),zlim=c(0,0.1),xlim=range(range+main.offset),xlab="age(y)",ylim=range(num.levels),ylab=names(pred.grid[2]))
					title(main=ynames[i])
				dev.off()
			}
		}
		## all y in 1 graph
		## need some sort of difference measure here
	}


	#######
	## 3 ##
	#######

	if(type==3){

		## for each y
		for(i in 1:dim(pred)[2]){

			## predicted fits
			if(figs[1]){
				## for setting y limits on graph
				pred.<-as.numeric(pred[,i])
				sim.pred.mean.<-as.numeric(sim.pred.mean[,i])
				sim.pred.sd.<-as.numeric(sim.pred.sd[,i])
				## calculate this in derivatives file?
				sim.pred.mean..<-sapply(1:length(range), function(x){
					if(x==length(range)){
						mean(sim.pred.mean.[which(1:length(sim.pred.mean.)%%length(range)==0)])
					}else{
						mean(sim.pred.mean.[which(1:length(sim.pred.mean.)%%length(range)==x)])
					}
				})
				sim.pred.sd..<-sapply(1:length(range), function(x){
					if(x==length(range)){
						mean(sim.pred.sd.[which(1:length(sim.pred.sd.)%%length(range)==0)])
					}else{
						mean(sim.pred.sd.[which(1:length(sim.pred.sd.)%%length(range)==x)])
					}
				})
				uniq<-unique(pred.grid[,2])

				pdf(paste(path,"pred",paste(i,"pdf",sep="."),sep="/"))
					## prepare plot
					yL<-min(c(min(pred.),min(sim.pred.mean.-sim.pred.sd.)))
					yH<-max(c(max(pred.),max(sim.pred.mean.+sim.pred.sd.)))
					plot(0,0,pch=30,main=ynames[i],xlim=range(range)+main.offset,xlab="age(y)",ylim=c(yL,yH),ylab="FA maturation index")
					## add null distribution
					polygon(c(range+main.offset,rev(range+main.offset)),c(sim.pred.mean..+sim.pred.sd..,rev(sim.pred.mean..-sim.pred.sd..)),col=rgb(0,0,0,75,maxColorValue=255),border=NA)
					## add predicted lines
					legend(min(range)+main.offset,yH,uniq,border=NULL,bty="n",fill=rainbow(length(uniq)))
					for(j in 1:length(uniq)) lines(range+main.offset,pred.[which(pred.grid[,2]==uniq[j])],col=rainbow(length(uniq))[j],lwd=2)
				dev.off()
			}
			
			## pvals (normal distribution)
			if(figs[2]){
				pdf(paste(path,"pred.p",paste(i,"pdf",sep="."),sep="/"))
					plot(0,0,pch=30,main=ynames[i],xlim=range(range)+main.offset,xlab="age(y)",ylim=c(0,0.1),ylab="pval")
					lines(range(range)+main.offset,c(0.05,0.05),lty=3)
					lines(range(range)+main.offset,c(0.01,0.01),lty=3)
					lines(range(range)+main.offset,c(0.001,0.001),lty=3)
					legend(min(range)+main.offset,0.1,uniq,border=NULL,bty="n",fill=rainbow(length(uniq)))
					for(j in 1:length(uniq)) lines(range+main.offset,pred.p[which(pred.grid[,2]==uniq[j]),i],col=rainbow(length(uniq))[j],lwd=2)
				dev.off()
			}
			
			## pvals (monte carlo)
			if(figs[3]){
				pdf(paste(path,"pred.pboot",paste(i,"pdf",sep="."),sep="/"))
					plot(0,0,pch=30,main=ynames[i],xlim=range(range)+main.offset,xlab="age(y)",ylim=c(0,0.1),ylab="pval")
					lines(range(range)+main.offset,c(0.05,0.05),lty=3)
					lines(range(range)+main.offset,c(0.01,0.01),lty=3)
					lines(range(range)+main.offset,c(0.001,0.001),lty=3)
					legend(min(range)+main.offset,0.1,uniq,border=NULL,bty="n",fill=rainbow(length(uniq)))
					for(j in 1:length(uniq)) lines(range+main.offset,pred.pboot[which(pred.grid[,2]==uniq[j]),i],col=rainbow(length(uniq))[j],lwd=2)
				dev.off()
			}
			
			## predicted fit derivatives
			if(figs[4]){
				## for setting y limits on graph
				pred.d.<-as.numeric(pred.d[,i])
				sim.pred.d.mean.<-as.numeric(sim.pred.d.mean[,i])
				sim.pred.d.sd.<-as.numeric(sim.pred.d.sd[,i])
				## calculate this in derivatives file?
				sim.pred.d.mean..<-sapply(1:length(range), function(x){
					if(x==length(range)){
						mean(sim.pred.d.mean.[which(1:length(sim.pred.d.mean.)%%length(range)==0)])
					}else{
						mean(sim.pred.d.mean.[which(1:length(sim.pred.d.mean.)%%length(range)==x)])
					}
				})
				sim.pred.d.sd..<-sapply(1:length(range), function(x){
					if(x==length(range)){
						mean(sim.pred.d.sd.[which(1:length(sim.pred.d.sd.)%%length(range)==0)])
					}else{
						mean(sim.pred.d.sd.[which(1:length(sim.pred.d.sd.)%%length(range)==x)])
					}
				})

				## prediction - derivatives
				pdf(paste(path,"pred.d",paste(i,"pdf",sep="."),sep="/"))
					## prepare plot
					yL<-min(c(min(pred.d.),min(sim.pred.d.mean.-sim.pred.d.sd.)))
					yH<-max(c(max(pred.d.),max(sim.pred.d.mean.+sim.pred.d.sd.)))
					plot(0,0,pch=30,main=ynames[i],xlim=range(range)+main.offset,xlab="age(y)",ylim=c(yL,yH),ylab="derivative of maturation index")
					## add null distribution
					polygon(c(range+main.offset,rev(range+main.offset)),c(sim.pred.d.mean..+sim.pred.d.sd..,rev(sim.pred.d.mean..-sim.pred.d.sd..)),col=rgb(0,0,0,75,maxColorValue=255),border=NA)
					## add predicted line
					legend(min(range)+main.offset,yH,uniq,border=NULL,bty="n",fill=rainbow(length(uniq)))
					for(j in 1:length(uniq)) lines(range+main.offset,pred.d.[which(pred.grid[,2]==uniq[j])],col=rainbow(length(uniq))[j],lwd=2)
				dev.off()
			}

			## derivative pvals (normal distribution)
			if(figs[5]){
				pdf(paste(path,"pred.d.p",paste(i,"pdf",sep="."),sep="/"))
					plot(0,0,pch=30,main=ynames[i],xlim=range(range)+main.offset,xlab="age(y)",ylim=c(0,0.1),ylab="pval")
					lines(range(range)+main.offset,c(0.05,0.05),lty=3)
					lines(range(range)+main.offset,c(0.01,0.01),lty=3)
					lines(range(range)+main.offset,c(0.001,0.001),lty=3)
					legend(min(range)+main.offset,0.1,uniq,border=NULL,bty="n",fill=rainbow(length(uniq)))
					for(j in 1:length(uniq)) lines(range+main.offset,pred.d.p[which(pred.grid[,2]==uniq[j]),i],col=rainbow(length(uniq))[j],lwd=2)
				dev.off()
			}

			## derivative pvals (monte carlo)
			if(figs[6]){
				pdf(paste(path,"pred.d.pboot",paste(i,"pdf",sep="."),sep="/"))
					plot(0,0,pch=30,main=ynames[i],xlim=range(range)+main.offset,xlab="age(y)",ylim=c(0,0.1),ylab="pval")
					lines(range(range)+main.offset,c(0.05,0.05),lty=3)
					lines(range(range)+main.offset,c(0.01,0.01),lty=3)
					lines(range(range)+main.offset,c(0.001,0.001),lty=3)
					legend(min(range)+main.offset,0.1,uniq,border=NULL,bty="n",fill=rainbow(length(uniq)))
					for(j in 1:length(uniq)) lines(range+main.offset,pred.d.pboot[which(pred.grid[,2]==uniq[j]),i],col=rainbow(length(uniq))[j],lwd=2)
				dev.off()
			}
		}

	}

	if(type>1){
		## all y in 1 graph
		if(dim(pred)[2]>1){

			if(figs[1]){
				pdf(paste(path,"pred.dif.pdf",sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(min(pred.dif),max(pred.dif),length.out=20)
					image(matrix(range.,1,length(range.)),col=topo.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,3))
					title(ylab="difference from null model")
					image(range+main.offset,1:dim(pred)[2],pred.dif,col=topo.colors(length(range.)),xlim=range(range+main.offset),xlab="age(y)",yaxt="n",ylab="Y")
					axis(2,1:length(ynames),ynames,las=1)	
				dev.off()
			}

			if(figs[2]){
				pdf(paste(path,"pred.dif.p.pdf",sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(0,0.1,length.out=20)
					image(matrix(range.,1,length(range.)),col=heat.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,2))
					title(ylab="p")
					image(range+main.offset,1:dim(pred)[2],pred.dif.p,col=heat.colors(length(range.)),zlim=c(0,0.1),xlim=range(range+main.offset),xlab="age(y)",yaxt="n",ylab="Y")
					axis(2,1:length(ynames),ynames,las=1)	
				dev.off()
			}

			if(figs[3]){
				pdf(paste(path,"pred.dif.pboot.pdf",sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(0,0.1,length.out=20)
					image(matrix(range.,1,length(range.)),col=heat.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,2))
					title(ylab="p")
					image(range+main.offset,1:dim(pred)[2],pred.dif.pboot,col=heat.colors(length(range.)),zlim=c(0,0.1),xlim=range(range+main.offset),xlab="age(y)",yaxt="n",ylab="Y")
					axis(2,1:length(ynames),ynames,las=1)	
				dev.off()
			}

			if(figs[4]){
				pdf(paste(path,"pred.d.dif.pdf",sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(min(pred.d.dif),max(pred.d.dif),length.out=20)
					image(matrix(range.,1,length(range.)),col=topo.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,3))
					title(ylab="maturation index derivative")
					image(range+main.offset,1:dim(pred)[2],pred.d.dif,col=topo.colors(length(range.)),xlim=range(range+main.offset),xlab="age(y)",yaxt="n",ylab="Y")
					axis(2,1:length(ynames),ynames,las=1)	
				dev.off()
			}

			if(figs[5]){
				pdf(paste(path,"pred.d.dif.p.pdf",sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(0,0.1,length.out=20)
					image(matrix(range.,1,length(range.)),col=heat.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,2))
					title(ylab="p")
					image(range+main.offset,1:dim(pred)[2],pred.d.dif.p,col=heat.colors(length(range.)),zlim=c(0,0.1),xlim=range(range+main.offset),xlab="age(y)",yaxt="n",ylab="Y")
					axis(2,1:length(ynames),ynames,las=1)	
				dev.off()
			}

			if(figs[6]){
				pdf(paste(path,"pred.d.dif.pboot.pdf",sep="/"))
					layout(matrix(c(2,1),ncol=2L),widths=c(4,1))
					range.<-seq(0,0.1,length.out=20)
					image(matrix(range.,1,length(range.)),col=heat.colors(length(range.)),xaxt="n",yaxt="n")
					axis(2,seq(0,1,length.out=length(range.)),round(range.,2))
					title(ylab="p")
					image(range+main.offset,1:dim(pred)[2],pred.d.dif.pboot,col=heat.colors(length(range.)),zlim=c(0,0.1),xlim=range(range+main.offset),xlab="age(y)",yaxt="n",ylab="Y")
					axis(2,1:length(ynames),ynames,las=1)
				dev.off()
			}
		}
	}

}
