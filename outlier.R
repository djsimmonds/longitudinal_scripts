outlier<-function(y,x=demo$age,x.apply=c(FALSE,TRUE),id=demo$id,sub=subset,within=TRUE,allY=TRUE,sd.thr=2,basename){
## function for dealing with outliers
## details:
	## Y - (required) - matrix of data, dimensions (# of scans, # of y measurements (such as voxels, ROIs, etc...)) 
	## X - variables to be used in regression for calculation of residuals (i.e. age) (each an array of length # of scans)
		## only implemented using linear regression
	## x.apply - use X for calculating residuals across sample (default=FALSE) and/or within subject (default=TRUE)
	## id - id for within-subject (or group) residual calculation (array of length # of scans) if applicable
	## sub - subset of scans to calculate outliers
	## within - should within subject residuals be calculated (TRUE/FALSE) (default=TRUE)
	## allY - should outliers across all y variables be calculated and removed before assessing within each y (TRUE/FALSE) (default=TRUE)
	## sd.thr - standard deviation threshold at which weighting is applied (default=2)
	## logfile, outfile

	## write to logfile
	sink(paste(basename,"log",sep="_"))

	Y<-as.matrix(y) ## in case of array

	## indexes NAs
	na.ind<-ifelse(is.na(rowMeans(Y)),0,1) ## Y
	if(!is.null(x)){
		X<-as.matrix(x) ## in case of array
		na.ind<-na.ind*ifelse(is.na(rowMeans(X)),0,1) ## X and Y
	}
	if(sum(na.ind)<3) stop("too many NAs - check your variables")

	## taking subset of subjects for calculations
	sub.ind<-which(sub*na.ind>0)
	if(!is.null(x)) X<-as.matrix(X[sub.ind,]) else X<-NULL
	Y<-as.matrix(Y[sub.ind,])
	if(within==TRUE) id<-id[sub.ind]

	## list to be returned with all outlier calculations
	out<-list()
	## variables
	out$X<-X
	out$Y<-Y
	out$id<-id
	## indices of included scans
	out$ind<-list()
	out$ind$sub<-which(sub>0)
	out$ind$na<-which(na.ind>0)
	out$ind$sub.na<-sub.ind
	
	## weighting based on magnitude of normalized residuals
	w<-matrix(1,dim(Y)[1],dim(Y)[2])
	w.all<-array(1,dim(Y)[1])
	w.all.bin<-array(1,dim(Y)[1])
	
	## defining list containing residuals and enclosed elements (different residual calculations)
	## whole=residuals calculated across whole sample, within=residuals calculated within subjects
	## t=residuals calculated for each time point, id=residuals calculated for each subject/group, id.t=residuals calculated within subject/group, but then assigned to each time point
	res.whole.t<-matrix(NA,dim(Y)[1],dim(Y)[2])
	if(within==TRUE){
		res.whole.id<-matrix(NA,length(unique(id)),dim(Y)[2])
		res.whole.id.t<-matrix(NA,dim(Y)[1],dim(Y)[2])
		res.within.t<-matrix(NA,dim(Y)[1],dim(Y)[2])
		res.within.id<-matrix(NA,length(unique(id)),dim(Y)[2])
		res.within.id.t<-matrix(NA,dim(Y)[1],dim(Y)[2])
	}

	## calculate for each y column
	for(i in 1:dim(Y)[2]){

		if(i>1) dt<-proc.time()[[3]]-t1
		if(i%%100==0) cat(i,"\t",dt,"\n")
		t1<-proc.time()[[3]]
	
		temp<-out.sub(Y[,i],rep(1,dim(Y)[1]),X,id,x.apply,within)
		if(within==TRUE){
			res.whole.t[,i]<-temp[[1]]
			res.whole.id[,i]<-temp[[2]]
			res.whole.id.t[,i]<-temp[[3]]
			res.within.t[,i]<-temp[[4]]
			res.within.id[,i]<-temp[[5]]
			res.within.id.t[,i]<-temp[[6]]
		}else{
			res.whole.t[,i]<-temp
		}
	}

	## if allY is TRUE, calculate across y residuals and recalculate within y excluding outliers
	if(allY==TRUE){
	
		## list for calculations across y's
		mean.res.whole.t<-rowMeans(res.whole.t)
		mean.res.whole.t[union(which(is.na(mean.res.whole.t)),which(is.nan(mean.res.whole.t)))]<-0 ## NAs converted to 0 so the time points aren't excluded (for example, when some subjects have only 1 scan but within==TRUE)
		w.all<-ifelse(mean.res.whole.t>sd.thr,1/mean.res.whole.t,1)
		w.all.bin<-ifelse(w.all==1,1,0)
		w.all.ind<-which(w.all.bin==1)
		
		## longitudinal calculations
		if(within==TRUE){
			mean.res.whole.id<-rowMeans(res.whole.id)
			mean.res.whole.id[union(which(is.na(mean.res.whole.id)),which(is.nan(mean.res.whole.id)))]<-0 ## NAs converted to 0
			mean.res.whole.id.t<-rowMeans(res.whole.id.t)
			mean.res.whole.id.t[union(which(is.na(mean.res.whole.id.t)),which(is.nan(mean.res.whole.id.t)))]<-0 ## NAs converted to 0
			mean.res.within.t<-rowMeans(res.within.t)
			mean.res.within.t[union(which(is.na(mean.res.within.t)),which(is.nan(mean.res.within.t)))]<-0 ## NAs converted to 0
			mean.res.within.id<-rowMeans(res.within.id)
			mean.res.within.id[union(which(is.na(mean.res.within.id)),which(is.nan(mean.res.within.id)))]<-0 ## NAs converted to 0
			mean.res.within.id.t<-rowMeans(res.within.id.t)
			mean.res.within.id.t[union(which(is.na(mean.res.within.id.t)),which(is.nan(mean.res.within.id.t)))]<-0 ## NAs converted to 0
			## multiplies all weights together 
			w.all<-w.all*ifelse(mean.res.whole.id.t>sd.thr,1/mean.res.whole.id.t,1)*ifelse(mean.res.within.t>sd.thr,1/mean.res.within.t,1)*ifelse(mean.res.within.id.t>sd.thr,1/mean.res.within.id.t,1)
			w.all.bin<-ifelse(w.all==1,1,0)
			w.all.ind<-which(w.all.bin==1)
			w.all.id.bin<-ifelse(mean.res.whole.id>sd.thr,0,1)*ifelse(mean.res.within.id>sd.thr,0,1)
			w.all.id.ind<-which(w.all.id.bin==1)
		
			## saves 1st iteration calculations into nested list
			out$all.res.whole.t<-res.whole.t
			out$all.res.whole.id<-res.whole.id
			out$all.res.whole.id.t<-res.whole.id.t
			out$all.res.within.t<-res.within.t
			out$all.res.within.id<-res.within.id
			out$all.res.within.id.t<-res.within.id.t
			out$mean.res.whole.t<-mean.res.whole.t
			out$mean.res.whole.id<-mean.res.whole.id
			out$mean.res.whole.id.t<-mean.res.whole.id.t
			out$mean.res.within.t<-mean.res.within.t
			out$mean.res.within.id<-mean.res.within.id
			out$mean.res.within.id.t<-mean.res.within.id.t
			out$w.all<-w.all
			out$w.all.bin<-w.all.bin
			out$w.all.ind<-w.all.ind
			out$w.all.id.bin<-w.all.id.bin
			out$w.all.id.ind<-w.all.id.ind
			
		}else{
			out$all.res.whole.t<-res.whole.t
			out$mean.res.whole.t<-mean.res.whole.t
			out$w.all<-w.all
			out$w.all.bin<-w.all.bin
			out$w.all.ind<-w.all.ind
		}
		
		## recalculate for each y column, excluding 2sd outliers across all regions in regression estimation
		for(i in 1:dim(Y)[2]){

			if(i>1) dt<-proc.time()[[3]]-t1
			if(i%%100==0) cat(i,"\t",dt,"\n")
			t1<-proc.time()[[3]]
	
			temp<-out.sub(Y[,i],w.all.bin,X,id,x.apply,within)
			if(within==TRUE){
				res.whole.t[,i]<-temp[[1]]
				res.whole.id[,i]<-temp[[2]]
				res.whole.id.t[,i]<-temp[[3]]
				res.within.t[,i]<-temp[[4]]
				res.within.id[,i]<-temp[[5]]
				res.within.id.t[,i]<-temp[[6]]
			}else{
				res.whole.t[,i]<-temp
			}
		}
	}
	
	## removing NAs if present
	res.whole.t[union(which(is.na(rowMeans(res.whole.t))),which(is.nan(rowMeans(res.whole.t)))),]<-0 ## NAs converted to 0
	if(within==TRUE){
		res.whole.id[union(which(is.na(rowMeans(res.whole.id))),which(is.nan(rowMeans(res.whole.id)))),]<-0 ## NAs converted to 0
		res.whole.id.t[union(which(is.na(rowMeans(res.whole.id.t))),which(is.nan(rowMeans(res.whole.id.t)))),]<-0 ## NAs converted to 0
		res.within.t[union(which(is.na(rowMeans(res.within.t))),which(is.nan(rowMeans(res.within.t)))),]<-0 ## NAs converted to 0
		res.within.id[union(which(is.na(rowMeans(res.within.id))),which(is.nan(rowMeans(res.within.id)))),]<-0 ## NAs converted to 0
		res.within.id.t[union(which(is.na(rowMeans(res.within.id.t))),which(is.nan(rowMeans(res.within.id.t)))),]<-0 ## NAs converted to 0
	}
	
	## calculation of weights
	w<-as.numeric(w.all)*ifelse(res.whole.t>sd.thr,1/res.whole.t,1)
	if(within==TRUE) w<-w*ifelse(res.whole.id.t>sd.thr,1/res.whole.id.t,1)*ifelse(res.within.t>sd.thr,1/res.within.t,1)*ifelse(res.within.id.t>sd.thr,1/res.within.id.t,1)

	if(within==TRUE){
		out$res.whole.t<-res.whole.t
		out$res.whole.id<-res.whole.id
		out$res.whole.id.t<-res.whole.id.t
		out$res.within.t<-res.within.t
		out$res.within.id<-res.within.id
		out$res.within.id.t<-res.within.id.t
		out$w<-w
	}else{
		out$res.whole.t<-res.whole.t
		out$w<-w
	}
	print(out)
	sink()
	save(out,file=paste(basename,"data",sep="_"))
	weights=list(ind=out$ind$sub.na,w=out$w)
	save(weights,file=paste(basename,"weights",sep="_"))
	weights
}
