out.sub<-function(y,w,X,id,type,within){

	w.ind<-which(w==1)
	
	## whole, t
	if(type=="mean") res<-(y-weighted.mean(y,w))^2 else res<-lm.wfit(X,y,w)$resid^2
	## normalizes squared residuals to sd=1, for easier calculation and comparison across region
	res.whole.t<-res/sd(res,na.rm=TRUE)
	
	## longitudinal calculations
	if(within==TRUE){
		res.whole.id<-array(NA,length(unique(id)))
		res.whole.id.t<-array(NA,length(y))
		res.within.t<-array(NA,length(y))
		res.within.id<-array(NA,length(unique(id)))
		res.within.id.t<-array(NA,length(y))
		
		for(j in 1:length(unique(id))){
			ind<-which(id==unique(id)[j])
			## whole, id
			res.whole.id[j]<-mean(res[ind])
			res.whole.id.t[ind]<-res.whole.id[j]
			## within, t (if there is 1 scan, puts in NA, and if 2, just calculates the residuals based on the mean)
			if(sum(w[ind])==1) res.<-NA else if(sum(w[ind])==2|type=="mean") res.<-(y[ind]-weighted.mean(y[ind],w[ind]))^2 else res.<-lm.wfit(as.matrix(X[ind,]),y[ind],w[ind])$resid^2
			## within, id
			res.within.id[j]<-mean(res.)
			res.within.id.t[ind]<-res.within.id[j]
		}
		
		## normalizes squared residuals to sd=1, for easier calculation and comparison across regions
		res.whole.id<-res.whole.id/sd(res.whole.id[w.ind],na.rm=TRUE)
		res.whole.id.t<-res.whole.id.t/sd(res.whole.id[w.ind],na.rm=TRUE)
		res.within.t<-res.within.t/sd(res.within.t[w.ind],na.rm=TRUE)
		res.within.id<-res.within.id/sd(res.within.id[w.ind],na.rm=TRUE)
		res.within.id.t<-res.within.id.t/sd(res.within.id.t,na.rm=TRUE)
	
		list(res.whole.t,res.whole.id,res.whole.id.t,res.within.t,res.within.id,res.within.id.t)
	}else{
		res.whole.t
	}
}
