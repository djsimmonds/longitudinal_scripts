## libraries
library(splines)
library(lme4)

## paths
paths<-list()
paths$study<-"/home/danisimmonds/Dani/dti_0511/tbss/analysis/behavior"
paths$scripts<-"/home/danisimmonds/Dropbox/scripts2"

## functions
llr.<-function(obj,obj.n,y){
	m<-refit(obj,y)
	m.n<-refit(obj.n,y)
	pchiqsq(-2*(logLik(m)-logLik(m.n)),obj$rank-obj.n$rank)
}

## load data and demographics
load(paste(paths$study,"demo.Rframe",sep="/"))
demo$age<-demo$age-8 ## to help with collinearity

## load outlier scripts
source(paste(paths$scripts,"outlier.R",sep="/"))
source(paste(paths$scripts,"out.sub.R",sep="/"))

## behavioral variables
vgs.mRT<-demo$vgs.mRT
vgs.cv<-demo$vgs.cv
anti.percErr<-demo$anti.percErr
anti.cv<-demo$anti.cv

## ranges
vgs.mRT.range<-quantile(vgs.mRT,1:10/10)
vgs.cv.range<-quantile(vgs.cv,1:10/10)
anti.percErr.range<-quantile(anti.percErr,1:10/10)
anti.cv.range<-quantile(anti.cv,1:10/10)

## indices for cross-sectional/longitudinal analyses
iC<-demo$nscans.ind==1
iL<-demo$nscans.tot>=4

## calculating outliers for each variable
vgs.mRT.cross.out<-outlier(vgs.mRT,sub=iC,within=FALSE,allY=FALSE,logfile="vgs.mRT.cross.out.log",outfile="vgs.mRT.cross.out.Rdata")
vgs.mRT.long.out<-outlier(vgs.mRT,sub=iL,allY=FALSE,logfile="vgs.mRT.long.out.log",outfile="vgs.mRT.long.out.Rdata")
vgs.cv.cross.out<-outlier(vgs.cv,sub=iC,within=FALSE,allY=FALSE,logfile="vgs.cv.cross.out.log",outfile="vgs.cv.cross.out.Rdata")
vgs.cv.long.out<-outlier(vgs.cv,sub=iL,allY=FALSE,logfile="vgs.cv.long.out.log",outfile="vgs.cv.long.out.Rdata")
anti.percErr.cross.out<-outlier(anti.percErr,sub=iC,within=FALSE,allY=FALSE,logfile="anti.percErr.cross.out.log",outfile="anti.percErr.cross.out.Rdata")
anti.percErr.long.out<-outlier(anti.percErr,sub=iL,allY=FALSE,logfile="anti.percErr.long.out.log",outfile="anti.percErr.long.out.Rdata")
anti.cv.cross.out<-outlier(anti.cv,sub=iC,within=FALSE,allY=FALSE,logfile="anti.cv.cross.out.log",outfile="anti.cv.cross.out.Rdata")
anti.cv.long.out<-outlier(anti.cv,sub=iL,allY=FALSE,logfile="anti.cv.long.out.log",outfile="anti.cv.long.out.Rdata")

## regression models


## PLOTS

## vgs.mRT
pdf("vgs.mRT.pdf")
par(mfrow=c(2,2),oma=c(1,1,3,1))
## cross
demo.<-demo[vgs.mRT.cross.out$ind$sub.na,]
w.<-vgs.mRT.cross.out$w==1
demo..<-demo.[w.,]
## males blue, females red
## children circles, adolescents triangles, adults diamonds
ind.m<-which(demo.$sex=="m")
ind.f<-which(demo.$sex=="f")
ind.chi<-which(demo.$age<13)
ind.ado<-intersect(which(demo.$age>=13),which(demo.$age<18))
ind.adu<-which(demo.$age>=18)
ind.m.<-which(demo..$sex=="m")
ind.f.<-which(demo..$sex=="f")
ind.chi.<-which(demo..$age<13)
ind.ado.<-intersect(which(demo..$age>=13),which(demo..$age<18))
ind.adu.<-which(demo..$age>=18)
## top left
plot(0,0,main="cross",xlab="age(y)",xlim=c(min(demo.$age),max(demo.$age)),ylab="vgs.mRT(ms)",ylim=c(min(demo.$vgs.mRT),max(demo.$vgs.mRT)))
points(demo.$age[intersect(ind.m,ind.chi)],demo.$vgs.mRT[intersect(ind.m,ind.chi)],col="blue",cex=0.5)
points(demo.$age[intersect(ind.m,ind.ado)],demo.$vgs.mRT[intersect(ind.m,ind.ado)],pch=2,col="blue",cex=0.5)
points(demo.$age[intersect(ind.m,ind.adu)],demo.$vgs.mRT[intersect(ind.m,ind.adu)],pch=5,col="blue",cex=0.5)
points(demo.$age[intersect(ind.f,ind.chi)],demo.$vgs.mRT[intersect(ind.f,ind.chi)],col="red",cex=0.5)
points(demo.$age[intersect(ind.f,ind.ado)],demo.$vgs.mRT[intersect(ind.f,ind.ado)],pch=2,col="red",cex=0.5)
points(demo.$age[intersect(ind.f,ind.adu)],demo.$vgs.mRT[intersect(ind.f,ind.adu)],pch=5,col="red",cex=0.5)
abline(lm(demo.$vgs.mRT[ind.m]~demo.$age),lty=2,col="blue")
abline(lm(demo.$vgs.mRT[ind.f]~demo.$age),lty=2,col="red")
## top right
plot(0,0,main="cross (out rem)",xlab="age(y)",xlim=c(min(demo.$age),max(demo.$age)),ylab="vgs.mRT(ms)",ylim=c(min(demo.$vgs.mRT),max(demo.$vgs.mRT)))
points(demo..$age[intersect(ind.m.,ind.chi.)],demo..$vgs.mRT[intersect(ind.m.,ind.chi.)],col="blue",cex=0.5)
points(demo..$age[intersect(ind.m.,ind.ado.)],demo..$vgs.mRT[intersect(ind.m.,ind.ado.)],pch=2,col="blue",cex=0.5)
points(demo..$age[intersect(ind.m.,ind.adu.)],demo..$vgs.mRT[intersect(ind.m.,ind.adu.)],pch=5,col="blue",cex=0.5)
points(demo..$age[intersect(ind.f.,ind.chi.)],demo..$vgs.mRT[intersect(ind.f.,ind.chi.)],col="red",cex=0.5)
points(demo..$age[intersect(ind.f.,ind.ado.)],demo..$vgs.mRT[intersect(ind.f.,ind.ado.)],pch=2,col="red",cex=0.5)
points(demo..$age[intersect(ind.f.,ind.adu.)],demo..$vgs.mRT[intersect(ind.f.,ind.adu.)],pch=5,col="red",cex=0.5)
abline(lm(demo..$vgs.mRT[ind.m.]~demo..$age),lty=2,col="blue")
abline(lm(demo..$vgs.mRT[ind.f.]~demo..$age),lty=2,col="red")
## cross
demo.<-demo[vgs.mRT.long.out$ind$sub.na,]
w.<-vgs.mRT.long.out$w==1
demo..<-demo.[w.,]
## males blue, females red
## children circles, adolescents triangles, adults diamonds
ind.m<-which(demo.$sex=="m")
ind.f<-which(demo.$sex=="f")
ind.chi<-which(demo.$age<13)
ind.ado<-intersect(which(demo.$age>=13),which(demo.$age<18))
ind.adu<-which(demo.$age>=18)
ind.m.<-which(demo..$sex=="m")
ind.f.<-which(demo..$sex=="f")
ind.chi.<-which(demo..$age<13)
ind.ado.<-intersect(which(demo..$age>=13),which(demo..$age<18))
ind.adu.<-which(demo..$age>=18)
## bottom left
plot(0,0,main="cross",xlab="age(y)",xlim=c(min(demo.$age),max(demo.$age)),ylab="vgs.mRT(ms)",ylim=c(min(demo.$vgs.mRT),max(demo.$vgs.mRT)))
points(demo.$age[intersect(ind.m,ind.chi)],demo.$vgs.mRT[intersect(ind.m,ind.chi)],col="blue",cex=0.5)
points(demo.$age[intersect(ind.m,ind.ado)],demo.$vgs.mRT[intersect(ind.m,ind.ado)],pch=2,col="blue",cex=0.5)
points(demo.$age[intersect(ind.m,ind.adu)],demo.$vgs.mRT[intersect(ind.m,ind.adu)],pch=5,col="blue",cex=0.5)
points(demo.$age[intersect(ind.f,ind.chi)],demo.$vgs.mRT[intersect(ind.f,ind.chi)],col="red",cex=0.5)
points(demo.$age[intersect(ind.f,ind.ado)],demo.$vgs.mRT[intersect(ind.f,ind.ado)],pch=2,col="red",cex=0.5)
points(demo.$age[intersect(ind.f,ind.adu)],demo.$vgs.mRT[intersect(ind.f,ind.adu)],pch=5,col="red",cex=0.5)
abline(lm(demo.$vgs.mRT[ind.m]~demo.$age),lty=2,col="blue")
abline(lm(demo.$vgs.mRT[ind.f]~demo.$age),lty=2,col="red")
## bottom right
plot(0,0,main="cross (out rem)",xlab="age(y)",xlim=c(min(demo.$age),max(demo.$age)),ylab="vgs.mRT(ms)",ylim=c(min(demo.$vgs.mRT),max(demo.$vgs.mRT)))
points(demo..$age[intersect(ind.m.,ind.chi.)],demo..$vgs.mRT[intersect(ind.m.,ind.chi.)],col="blue",cex=0.5)
points(demo..$age[intersect(ind.m.,ind.ado.)],demo..$vgs.mRT[intersect(ind.m.,ind.ado.)],pch=2,col="blue",cex=0.5)
points(demo..$age[intersect(ind.m.,ind.adu.)],demo..$vgs.mRT[intersect(ind.m.,ind.adu.)],pch=5,col="blue",cex=0.5)
points(demo..$age[intersect(ind.f.,ind.chi.)],demo..$vgs.mRT[intersect(ind.f.,ind.chi.)],col="red",cex=0.5)
points(demo..$age[intersect(ind.f.,ind.ado.)],demo..$vgs.mRT[intersect(ind.f.,ind.ado.)],pch=2,col="red",cex=0.5)
points(demo..$age[intersect(ind.f.,ind.adu.)],demo..$vgs.mRT[intersect(ind.f.,ind.adu.)],pch=5,col="red",cex=0.5)
abline(lm(demo..$vgs.mRT[ind.m.]~demo..$age),lty=2,col="blue")
abline(lm(demo..$vgs.mRT[ind.f.]~demo..$age),lty=2,col="red")


mtext("vgs.mRT",outer=TRUE,line=1)
dev.off()


## regression

## vgs.mRT

## cross-sectional
demo.<-demo[vgs.mRT.cross.out$ind$sub.na,]
demo.$w<-vgs.mRT.cross.out$w

m0<-lm(vgs.mRT~1,demo.,weights=w) ## null
m1<-lm(vgs.mRT~age,demo.,weights=w) ## lin
m2<-lm(vgs.mRT~ns(age,2),demo.,weights=w) ## nonlinear (spline)

## lin vs. null
pchisq(2*(logLik(m1)-logLik(m0)),m1$rank-m0$rank)
## spl vs. null
pchisq(2*(logLik(m2)-logLik(m0)),m2$rank-m0$rank)
## spl vs. lin
pchisq(2*(logLik(m2)-logLik(m1)),m2$rank-m1$rank)

