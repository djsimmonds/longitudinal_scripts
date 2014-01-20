out.time<-function(y,x=demo$age,id=demo$id,sub=subset,within=TRUE,allY=TRUE,sd.thr=2){
## function for dealing with outliers
## details:
	## Y - (required) - matrix of data, dimensions (# of scans, # of y measurements (such as voxels, ROIs, etc...)) 
	## X - variables to be used in regression for calculation of residuals (i.e. age) (each an array of length # of scans)
		## only implemented using linear regression
	## id - id for within-subject (or group) residual calculation (array of length # of scans) if applicable
	## sub - subset of scans to calculate outliers
	## within - should within subject residuals be calculated (TRUE/FALSE) (default=TRUE)
	## allY - should outliers across all y variables be calculated and removed before assessing within each y (TRUE/FALSE) (default=TRUE)
	## sd.thr - standard deviation threshold at which weighting is applied (default=2)

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
	X<-as.matrix(X[sub.ind,])
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

	## types for residual calculation
	if(is.null(X)) type="mean" else type="lm"

	## calculate for each y column
	for(i in 1:110){
		if(i==10) t1<-proc.time()[[3]] ## short burn-in before estimate
		temp<-out.sub(Y[,i],w[,i],X,id,type,within)
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

	dt<-proc.time()[[3]]-t1
	
	tot.<-(dt/100)*dim(Y)[2]*(allY+1)
	if(tot.<60){
		cat("Estimated time:",tot.,"seconds\n")
	}else if(tot.<60*60){
		cat("Estimated time:",tot./60,"minutes\n")
	}else if(tot.<60*60*24){
		cat("Estimated time:",tot./60/60,"hours\n")
	}else{
		cat("Estimated time:",tot./60/60/24,"days\n")
	}
}
