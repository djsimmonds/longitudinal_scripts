out.sub<-function(y,w,X,id,x.apply,within){

	w.ind<-which(w==1)
	w.ind.id<-which(is.element(unique(id),unique(id[which(w==1)])))
	
	## whole, t
	if(x.apply[1]==FALSE) res<-(y-weighted.mean(y,w))^2 else res<-lm.wfit(X,y,w)$resid^2
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
			if(sum(w[ind])<=1) res.<-NA else if(sum(w[ind])==2|x.apply[2]==FALSE) res.<-(y[ind]-weighted.mean(y[ind],w[ind]))^2 else res.<-lm.wfit(as.matrix(X[ind,]),y[ind],w[ind])$resid^2
			## within, id
			res.within.t[ind]<-res.
			res.within.id[j]<-mean(res.,na.rm=TRUE)
			res.within.id.t[ind]<-res.within.id[j]
		}
		
		## normalizes squared residuals to sd=1, for easier calculation and comparison across regions
		temp.sd<-sd(res.whole.id[w.ind.id],na.rm=TRUE)
		res.whole.id<-res.whole.id/temp.sd
		res.whole.id.t<-res.whole.id.t/temp.sd
		res.within.t<-res.within.t/sd(res.within.t[w.ind],na.rm=TRUE)
		temp.sd<-sd(res.within.id[w.ind.id],na.rm=TRUE)		
		res.within.id<-res.within.id/temp.sd
		res.within.id.t<-res.within.id.t/temp.sd
	
		list(res.whole.t,res.whole.id,res.whole.id.t,res.within.t,res.within.id,res.within.id.t)
	}else{
		res.whole.t
	}
}
