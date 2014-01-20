## libraries
#library(splines)

## load data and demographics
#load("/home/danisimmonds/Dani/dti_0511/tbss/analysis/cross/vox/data.Rdata")
#load("/home/danisimmonds/Dani/dti_0511/tbss/analysis/cross/vox/demo.Rdata")
#demo<-demo.
#load("/home/danisimmonds/Dani/dti_0511/tbss/atlas.ind.Rframe")

## data for all WM and specific tract types
#data.all<-rowMeans(data.[,-atlas.ind$L1.type[[1]]])
#data.P<-rowMeans(data.[,atlas.ind$L1.type$P])
#data.Ce<-rowMeans(data.[,atlas.ind$L1.type$Ce])
#data.AL<-rowMeans(data.[,atlas.ind$L1.type$AL])
#data.A<-rowMeans(data.[,atlas.ind$L1.type$A])
#data.Ca<-rowMeans(data.[,atlas.ind$L1.type$Ca])
#data<-data.frame(all=data.all,P=data.P,Ce=data.Ce,AL=data.AL,A=data.A,Ca=data.Ca)

## PLOTS

## anti.percErr
pdf("anti.percErr.pdf")
par(mfrow=c(3,2),oma=c(1,1,3,1))
## males blue, females red
## children circles, adolescents triangles, adults diamonds
demo.<-demo[!is.na(demo$anti.percErr),]
data.<-data[!is.na(demo$anti.percErr),]
ind.m<-which(demo.$sex=="m")
ind.f<-which(demo.$sex=="f")
ind.chi<-which(demo.$age<13)
ind.ado<-intersect(which(demo.$age>=13),which(demo.$age<18))
ind.adu<-which(demo.$age>=18)

## data.all
plot(0,0,main="all WM tracts",xlab="% anti errors",xlim=c(min(demo.$anti.percErr),max(demo.$anti.percErr)),ylab="MI",ylim=c(min(data.[,1]),max(data.[,1])))
points(demo.$anti.percErr[intersect(ind.m,ind.chi)],data.[intersect(ind.m,ind.chi),1],col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.m,ind.ado)],data.[intersect(ind.m,ind.ado),1],pch=2,col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.m,ind.adu)],data.[intersect(ind.m,ind.adu),1],pch=5,col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.chi)],data.[intersect(ind.f,ind.chi),1],col="red",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.ado)],data.[intersect(ind.f,ind.ado),1],pch=2,col="red",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.adu)],data.[intersect(ind.f,ind.adu),1],pch=5,col="red",cex=0.5)
abline(lm(data.[ind.m,1]~demo.$anti.percErr[ind.m]),lty=2,col="blue")
abline(lm(data.[ind.f,1]~demo.$anti.percErr[ind.f]),lty=2,col="red")

## data.P
plot(0,0,main="projection tracts",xlab="% anti errors",xlim=c(min(demo.$anti.percErr),max(demo.$anti.percErr)),ylab="MI",ylim=c(min(data.[,2]),max(data.[,2])))
points(demo.$anti.percErr[intersect(ind.m,ind.chi)],data.[intersect(ind.m,ind.chi),2],col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.m,ind.ado)],data.[intersect(ind.m,ind.ado),2],pch=2,col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.m,ind.adu)],data.[intersect(ind.m,ind.adu),2],pch=5,col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.chi)],data.[intersect(ind.f,ind.chi),2],col="red",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.ado)],data.[intersect(ind.f,ind.ado),2],pch=2,col="red",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.adu)],data.[intersect(ind.f,ind.adu),2],pch=5,col="red",cex=0.5)
abline(lm(data.[ind.m,2]~demo.$anti.percErr[ind.m]),lty=2,col="blue")
abline(lm(data.[ind.f,2]~demo.$anti.percErr[ind.f]),lty=2,col="red")

## data.Ce
plot(0,0,main="cerebellar tracts",xlab="% anti errors",xlim=c(min(demo.$anti.percErr),max(demo.$anti.percErr)),ylab="MI",ylim=c(min(data.[,3]),max(data.[,3])))
points(demo.$anti.percErr[intersect(ind.m,ind.chi)],data.[intersect(ind.m,ind.chi),3],col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.m,ind.ado)],data.[intersect(ind.m,ind.ado),3],pch=2,col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.m,ind.adu)],data.[intersect(ind.m,ind.adu),3],pch=5,col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.chi)],data.[intersect(ind.f,ind.chi),3],col="red",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.ado)],data.[intersect(ind.f,ind.ado),3],pch=2,col="red",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.adu)],data.[intersect(ind.f,ind.adu),3],pch=5,col="red",cex=0.5)
abline(lm(data.[ind.m,3]~demo.$anti.percErr[ind.m]),lty=2,col="blue")
abline(lm(data.[ind.f,3]~demo.$anti.percErr[ind.f]),lty=2,col="red")

## data.AL
plot(0,0,main="assoc/limbic tracts",xlab="% anti errors",xlim=c(min(demo.$anti.percErr),max(demo.$anti.percErr)),ylab="MI",ylim=c(min(data.[,4]),max(data.[,4])))
points(demo.$anti.percErr[intersect(ind.m,ind.chi)],data.[intersect(ind.m,ind.chi),4],col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.m,ind.ado)],data.[intersect(ind.m,ind.ado),4],pch=2,col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.m,ind.adu)],data.[intersect(ind.m,ind.adu),4],pch=5,col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.chi)],data.[intersect(ind.f,ind.chi),4],col="red",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.ado)],data.[intersect(ind.f,ind.ado),4],pch=2,col="red",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.adu)],data.[intersect(ind.f,ind.adu),4],pch=5,col="red",cex=0.5)
abline(lm(data.[ind.m,4]~demo.$anti.percErr[ind.m]),lty=2,col="blue")
abline(lm(data.[ind.f,4]~demo.$anti.percErr[ind.f]),lty=2,col="red")

## data.A
plot(0,0,main="association tracts",xlab="% anti errors",xlim=c(min(demo.$anti.percErr),max(demo.$anti.percErr)),ylab="MI",ylim=c(min(data.[,5]),max(data.[,5])))
points(demo.$anti.percErr[intersect(ind.m,ind.chi)],data.[intersect(ind.m,ind.chi),5],col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.m,ind.ado)],data.[intersect(ind.m,ind.ado),5],pch=2,col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.m,ind.adu)],data.[intersect(ind.m,ind.adu),5],pch=5,col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.chi)],data.[intersect(ind.f,ind.chi),5],col="red",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.ado)],data.[intersect(ind.f,ind.ado),5],pch=2,col="red",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.adu)],data.[intersect(ind.f,ind.adu),5],pch=5,col="red",cex=0.5)
abline(lm(data.[ind.m,5]~demo.$anti.percErr[ind.m]),lty=2,col="blue")
abline(lm(data.[ind.f,5]~demo.$anti.percErr[ind.f]),lty=2,col="red")

## data.Ca
plot(0,0,main="callosal tracts",xlab="% anti errors",xlim=c(min(demo.$anti.percErr),max(demo.$anti.percErr)),ylab="MI",ylim=c(min(data.[,6]),max(data.[,6])))
points(demo.$anti.percErr[intersect(ind.m,ind.chi)],data.[intersect(ind.m,ind.chi),6],col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.m,ind.ado)],data.[intersect(ind.m,ind.ado),6],pch=2,col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.m,ind.adu)],data.[intersect(ind.m,ind.adu),6],pch=5,col="blue",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.chi)],data.[intersect(ind.f,ind.chi),6],col="red",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.ado)],data.[intersect(ind.f,ind.ado),6],pch=2,col="red",cex=0.5)
points(demo.$anti.percErr[intersect(ind.f,ind.adu)],data.[intersect(ind.f,ind.adu),6],pch=5,col="red",cex=0.5)
abline(lm(data.[ind.m,6]~demo.$anti.percErr[ind.m]),lty=2,col="blue")
abline(lm(data.[ind.f,6]~demo.$anti.percErr[ind.f]),lty=2,col="red")
mtext("anti.percErr",outer=TRUE,line=1)
dev.off()


