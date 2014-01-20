ll.test<-function(Y=data.,W=Wts,W.=wts.,m1.=models$model[[m]],m0.=models$model[[models$con[[m]][c]]],mixed=setup$mixed,ynames=setup$ynames){
	## gets indices of rows which are excluded in m1
	na<-m1.$na
	## if any excluded, adjust data and weight matrices
	## note: W. (an additional weighting matrix, e.g. for behavioral variables) already has NAs removed by the outlier script
	if(length(na)>0){
		Y<-Y[-na,]
		W<-as.matrix(W[-na,])
	}
	if(!is.null(W.)) W<-apply(W,2,"*",W.)
	Y<-as.matrix(Y)
	W<-as.matrix(W)
	## data frame for m1 (NAs already removed)
	d<-m1.$data

	## updates models with appropriate data frame
	if(mixed==TRUE){
		m0<-update(m0.$fit,data=d,REML=FALSE)
		m1<-update(m1.$fit,data=d,REML=FALSE)
	}else{
		m0<-update(m0.$fit,data=d)
		m1<-update(m1.$fit,data=d)
	}

	## estimates

	## longitudinal
	if(mixed==TRUE){
		sapply(
			1:dim(Y)[2],
			function(i){
				w<-W[,i]
				m0@pWt<-w ## weights
				m0<-refit(m0,Y[,i])
				m1@pWt<-w ## weights
				m1<-refit(m1,Y[,i])
				## log likelihood test
				m0.ll<-logLik(m0)
				m1.ll<-logLik(m1)
				llr<-2*(as.numeric(m1.ll)-as.numeric(m0.ll))
				df<-attr(m1.ll,"df")-attr(m0.ll,"df")
				p<-1-pchisq(llr,df)
				c(llr,df,p)
			}
		)
		
	## cross-sectional
	}else{
		X0<-model.matrix(m0)
		X<-model.matrix(m1)
		sapply(
			1:dim(Y)[2],
			function(i){
				m0<-lm.wfit(X0,Y[,i],W[,i])
				m1<-lm.wfit(X,Y[,i],W[,i])
				llr<-2*(logLik.(m1)-logLik.(m0))
				df<-m1$rank-m0$rank
				p<-1-pchisq(llr,df)
				c(llr,df,p)
			}
		)
	}
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
