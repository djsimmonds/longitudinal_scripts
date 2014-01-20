ll.test<-function(
	y=Y,
	m1.=models$model[[m]],
	m0.=models$model[[models$con[[m]][c]]],
	mixed=setup$mixed,
	ynames=setup$ynames
){

	## gets indices of rows which are excluded in m1; if any excluded, adjust data matrix
	if(length(m1.$all.exc)>0) Y<-as.matrix(y[-m1.$all.exc,]) else Y<-y
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
		temp<-sapply(
			1:dim(Y)[2],
			function(i){
				m0<-refit(m0,Y[,i])
				m1<-refit(m1,Y[,i])
				## log likelihood test
				m0.ll<-logLik(m0)
				m1.ll<-logLik(m1)
				c(2*(as.numeric(m1.ll)-as.numeric(m0.ll)),attr(m1.ll,"df")-attr(m0.ll,"df"))
			}		
		)
		
	## cross-sectional
	}else{
		X0<-model.matrix(m0)
		X<-model.matrix(m1)
		temp<-sapply(
			1:dim(Y)[2],
			function(i){
				m0<-lm.fit(X0,Y[,i])
				m1<-lm.fit(X,Y[,i])
				c(2*(logLik.(m1)-logLik.(m0)),m1$rank-m0$rank)
			}
		)
	}
	llr<-temp[1,]
	df<-temp[2,]
	p<-1-pchisq(llr,df)
	cbind(llr,df,p)
}

logLik.<-function(obj){
  res <- obj$residuals
  N <- length(res)
	0.5*(-N*(log(2*pi)+1-log(N)+log(sum(res^2))))
}
