## paths
paths=list()
paths$study<-"~/Dani/dti_0511"
paths$analysis<-paste(paths$study,"tbss.dtitk","analysis","cross","vox",sep="/")

## load data
load(paste(paths$analysis,"data.Rdata",sep="/"))
load(paste(paths$analysis,"demo.Rdata",sep="/"))
load(paste(paths$analysis,"weights.Rdata",sep="/"))

## functions
aic.<-function(obj,y){
	m<-refit(obj,y)
	-2*logLik(m)+2*obj$rank
}

llr.<-function(obj,obj.n,y){
	m<-refit(obj,y)
	m.n<-refit(obj.n,y)
	pchiqsq(-2*(logLik(m)-logLik(m.n)),obj$rank-obj.n$rank)
}

aic.q<-function(y,w,X){
	if(!is.null(w)) m<-lm.wfit(X,y,w) else m<-lm.fit(X,y)
	-2*logLik.(m)+2*m$rank
}

llr.q<-function(y,w,X,X.n){
	if(!is.null(w)){
		m<-lm.wfit(X,y,w)
		m.n<-lm.wfit(X.n,y,w)
	}else{
		m<-lm.fit(X,y)
		m.n<-lm.fit(X.n,y)
	}
	pchisq(2*(logLik.(m)-logLik.(m.n)),m$rank-m.n$rank)
}

logLik.<-function(obj){
  res <- obj$residuals
  N <- length(res)
  if(is.null(w <- obj$weights)) {
    w <- rep.int(1, N)
  } else {
    excl <- w == 0
    if (any(excl)) {
		  res <- res[!excl]
		  N <- length(res)
		  w <- w[!excl]
    }
  }
	0.5*(sum(log(w))-N*(log(2*pi)+1-log(N)+log(sum(w*res^2))))
}

## age distribution for different analyses

## ALL
#> quantile(demo$age)
#   0%   25%   50%   75%  100% 
# 8.14 12.69 15.46 18.45 28.17 

## LONGITUDINAL
#> quantile(demo$age[which(demo$nscans.tot>=4)])
#     0%     25%     50%     75%    100% 
# 9.4600 14.1850 17.0850 19.4675 24.7200 

## CROSS-SECTIONAL
#> quantile(demo$age[which(demo$nscans.ind==1)])
#   0%   25%   50%   75%  100% 
# 8.14 11.59 13.93 16.53 28.17 

## k=(13,18)
	## long: (.152,.588)
	## cross: (.442,.849)

## finding best spline

## 1) knots based on percentiles

X.=list(as.matrix(rep(1,dim(data.)[1])),cbind(1,demo.$age))
for(i in 2:10) X.[[i+1]]<-cbind(1,ns(demo.$age,df=i))
pt<-proc.time()[[3]]
aic.mat=sapply(
	1:length(X.),
	function(i){
		X=X.[[i]]
		sapply(1:dim(data.)[2], function(j) aic.q(data.[,j],w.[,j],X)) 
	}
)
proc.time()[[3]]-pt
#[1] 103.175

f<-function(x) x-min(x)
aic.mat<-t(apply(aic.mat,1,f))
aic.mean<-colMeans(aic.mat)
aic0<-sapply(1:length(aic.mean),function(i) length(which(aic.mat[,i]==0)))
aic2<-sapply(1:length(aic.mean),function(i) length(which(aic.mat[,i]<2)))
#> cbind(aic.mean,aic0,aic2)
#				 aic.mean  aic0  aic2
# [1,] 14.6271153  1254  2663
# [2,]  0.7362049 16503 21880
# [3,]  1.6092503  2937 21476
# [4,]  2.4453610  1975  8933
# [5,]  3.5375374   754  4785
# [6,]  4.6904032   426  2607
# [7,]  5.8875242   259  1379
# [8,]  7.0476460   172   882
# [9,]  8.0968833   151   677
#[10,]  9.0619733   220   618
#[11,]  9.9913104   229   643

## 2) knots based on set range
	## 1- 16
	## 2- 13,18
	## 3- 12,15,18
	## 4- 12,14,17,19
	## 5- 11,13,15,17,20

X.=list(as.matrix(rep(1,dim(data.)[1])),cbind(1,demo.$age))
X.[[3]]<-cbind(1,ns(demo.$age,k=16))
X.[[4]]<-cbind(1,ns(demo.$age,k=c(13,18)))
X.[[5]]<-cbind(1,ns(demo.$age,k=c(12,15,18)))
X.[[6]]<-cbind(1,ns(demo.$age,k=c(12,14,17,19)))
X.[[7]]<-cbind(1,ns(demo.$age,k=c(11,13,15,17,20)))

aic.mat=sapply(
	1:length(X.),
	function(i){
		X=X.[[i]]
		sapply(1:dim(data.)[2], function(j) aic.q(data.[,j],w.[,j],X)) 
	}
)

f<-function(x) x-min(x)
aic.mat<-t(apply(aic.mat,1,f))
aic.mean<-colMeans(aic.mat)
aic0<-sapply(1:length(aic.mean),function(i) length(which(aic.mat[,i]==0)))
aic2<-sapply(1:length(aic.mean),function(i) length(which(aic.mat[,i]<2)))
#> cbind(aic.mean,aic0,aic2)
#       aic.mean  aic0  aic2
#[1,] 14.5194529  1280  2738
#[2,]  0.6285424 17155 22300
#[3,]  1.5842272  2575 21734
#[4,]  2.3572621  2344  9182
#[5,]  3.5228911   811  4847
#[6,]  4.7298246   343  2274
#[7,]  5.8731928   372  1453

## df = N*T-N-1 (# total scans - # subjects - 1)
## organization within analysis folder

## FOR ROI
> model
	> t/chisq file (1 entry/voxel or roi)
	> df file (only need 1 value)
	> pval file (1 entry/voxel or roi)\
	> deriv
		> file of predicted values
		> mean of null predicted values
		> sd of null predicted values
		> pvals of actual vs. null predicted values
		
## FOR VOXELWISE
> model
	> t/chisq file (1 entry/voxel or roi)
	> df file (only need 1 value)
	> pval file (1 entry/voxel or roi)
