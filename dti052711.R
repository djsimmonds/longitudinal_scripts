load("atlas.ind.Rframe")
y.ind<-cbind(c(1,1,1,1,1,5,5,5,5,5,6,6,6,6,6,7),c(2,3,4,5,6,1,3,4,5,6,2,3,4,5,6,2))
Y.<-sapply(1:dim(y.ind)[1], function(i) rowMeans(Y[,atlas.ind[[y.ind[i,1]]][[y.ind[i,2]]]]))

range<-data.frame(x=seq(ceiling(min(demo.$age)*10)/10,floor(max(demo.$age)*10)/10,0.1))
pdf("/home/danisimmonds/Dropbox/fa_roi16.pdf")
plot(0,0,xlim=c(min(demo.$age),max(demo.$age)),xlab="age(y)",ylim=c(min(data..),max(data..)),ylab="FA MI")
	for(i in 1:dim(Y)[2]){ ## 16
		d<-data.frame(y=data..[,i],x=demo.$age,w=W[[2]][,i])
		m<-lm(y~ns(x,k=c(13,18,23)),d,weights=w)
		points(demo.$age,data..[,i],cex=0.5,col=rainbow(16)[i])
		lines(range[,1],predict(m,range),col=rainbow(16)[i])
	}
dev.off()

pred<-sapply(
	1:dim(data.)[2], 
	function(i){
		d<-data.frame(y=data.[,i],x=demo.$age,w=W[[2]][,i])
		m<-lm(y~ns(x,k=c(13,18,23)),d,weights=w)
		predict(m,range)
	}
)

range.<-rep(range[,1],each=2)
ind<-which((1:length(range.))%%2==0)
range.[ind]<-range.[ind]+0.0001
range.<-data.frame(x=range.)

pred.d<-sapply(
	1:dim(data.)[2], 
	function(i){
		d<-data.frame(y=data.[,i],x=demo.$age,w=W[[2]][,i])
		m<-lm(y~ns(x,k=c(13,18,23)),d,weights=w)
		pred.<-predict(m,range.)
		diff(as.numeric(pred.))[ind-1]*10000
	}
)

pdf("/home/danisimmonds/Dropbox/fa_roi16_pred_deriv.pdf")
heatmap(t(pred.d),Rowv=NA,Colv=NA,labRow=ynames,labCol=range[,1],scale="none")
dev.off()


ynames<-sapply(1:16, function(x) paste(names(atlas.ind)[y.ind[x,1]],names(atlas.ind[[y.ind[x,1]]])[y.ind[x,2]]))

## I SHOULD CONSIDER CALCULATING BOTH TYPES OF P-VALUES FOR BOOTSTRAPS
	#pnorm(x,mean(sim),sd(sim)) --> is this valid?
	#rowMeans(sim<x) --> significance can only go up to # simulations
