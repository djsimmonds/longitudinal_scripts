outlier<-function(y,x=demo$age,id=demo$id,sub=subset,within=TRUE,allY=TRUE,sd.thr=2){
## function for dealing with outliers
## details:
	## Y - (required) - matrix of data, dimensions (# of scans, # of y measurements (such as voxels, ROIs, etc...)) 
	## X - variables to be used in regression for calculation of residuals (i.e. age) (each an array of length # of scans)
		## only implemented using linear regression
	## id - id for within-subject (or group) residual calculation (array of length # of scans) if applicable
	## within - should within subject residuals be calculated (TRUE/FALSE) (default=TRUE)
	## allY - should outliers across all y variables be calculated and removed before assessing within each y (TRUE/FALSE) (default=TRUE)
	## sd.thr - standard deviation threshold at which weighting is applied (default=2)

	Y<-as.matrix(y) ## in case of array

	## indexes NAs
	na.ind<-ifelse(is.na(rowMeans(Y)),0,1) ## Y
	if(!is.null(x)){
		X<-as.matrix(x)
		na.ind<-na.ind*ifelse(is.na(rowMeans(X)),0,1) ## X and Y
	}
	sub.<-which(sub*na.ind>0)
	if(sum(na.ind)<3) stop("too many NAs - check your variables")
	X<-as.matrix(X[sub.,])
	Y<-as.matrix(Y[sub.,])
	id<-id[na.ind]

	## taking subset of subjects for calculations
	if(!is.null(X)){
		X<-as.matrix(x) ## in case it is an array
		X<-as.matrix(X[sub,]) ## subset, again with array precautions
	}
	Y<-as.matrix(y)
	Y<-as.matrix(y[sub,])
	id<-id[sub]


	## list to be returned with all outlier calculations
	out<-list()
	## weighting based on magnitude of normalized residuals
	out$w<-matrix(1,dim(Y)[1],dim(Y)[2])
	## defining list containing residuals and enclosed elements (different residual calculations)
	## whole=residuals calculated across whole sample, within=residuals calculated within subjects
	## t=residuals calculated for each time point, id=residuals calculated for each subject/group, id.t=residuals calculated within subject/group, but then assigned to each time point
	out$res<-list()
	out$res$whole<-list()
	out$res$whole$t<-matrix(NA,dim(Y)[1],dim(Y)[2])
	if(within==TRUE){
		out$res$whole$id<-matrix(NA,length(unique(id)),dim(Y)[2])
		out$res$whole$id.t<-matrix(NA,dim(Y)[1],dim(Y)[2])
		out$res$within<-list()
		out$res$within$t<-matrix(NA,dim(Y)[1],dim(Y)[2])
		out$res$within$id<-matrix(NA,length(unique(id)),dim(Y)[2])
		out$res$within$id.t<-matrix(NA,dim(Y)[1],dim(Y)[2])
	}
	
	## formulas for residual calculation
	if(!is.null(X)) f<-as.formula("Y[,i]~X") else f<-as.formula("Y[,i]~1")
	if(within==TRUE) if(!is.null(X)) f.w<-as.formula("Y[ind,i]~X[ind]") else f.w<-as.formula("Y[ind,i]~1")

	t1=proc.time()[[3]]

	## calculate for each y column
	for(i in 1:dim(Y)[2]){

		dt<-proc.time()[[3]]-t1
		print(dt)
		t1<-proc.time()[[3]]
	
		## whole, t
		out$res$whole$t[,i]<-resid(lm(f))^2
		## normalizes squared residuals to sd=1, for easier calculation and comparison across regions
		out$res$whole$t[,i]<-out$res$whole$t[,i]/sd(out$res$whole$t[,i])
		
		## longitudinal calculations
		if(within==TRUE){
			for(j in 1:length(unique(id))){
				ind<-which(id==unique(id)[j])
				## whole, id
				out$res$whole$id[j,i]<-mean(out$res$whole$t[ind,i],na.rm=TRUE)
				out$res$whole$id.t[ind,i]<-out$res$whole$id[j,i]
				## within, t (if there is 1 scan, puts in NA, and if 2, just calculates the residuals based on the mean)
				if(length(ind)==1) out$res$within$t[ind,i]<-NA else if(length(ind)==2) out$res$within$t[ind,i]<-resid(lm(as.formula("Y[ind,i]~1")))^2 else out$res$within$t[ind,i]<-resid(lm(f.w))^2
				## within, id
				out$res$within$id[j,i]<-mean(out$res$within$t[ind,i],na.rm=TRUE)
				out$res$within$id.t[ind,i]<-out$res$within$id[j,i]
			}
			## normalizes squared residuals to sd=1, for easier calculation and comparison across regions
			out$res$whole$id[,i]<-out$res$whole$id[,i]/sd(out$res$whole$id[,i],na.rm=TRUE)
			out$res$whole$id.t[,i]<-out$res$whole$id.t[,i]/sd(out$res$whole$id[,i],na.rm=TRUE)
			out$res$within$t[,i]<-out$res$within$t[,i]/sd(out$res$within$t[,i],na.rm=TRUE)
			out$res$within$id[,i]<-out$res$within$id[,i]/sd(out$res$within$id[,i],na.rm=TRUE)
			out$res$within$id.t[,i]<-out$res$within$id.t[,i]/sd(out$res$within$id[,i],na.rm=TRUE)
		}
	}

	## if allY is TRUE, calculate across y residuals and recalculate within y excluding outliers
	if(allY==TRUE){
	
		## list for calculations across y's
		out$resAll<-list()
		out$resAll$whole<-list()
		out$resAll$whole$t<-rowMeans(out$res$whole$t,na.rm=TRUE)
		out$resAll$whole$t[which(is.na(out$resAll$whole$t))]<-0 ## NAs converted to 0 so the time points aren't excluded (for example, when some subjects have only 1 scan but within==TRUE)
		out$w[,]<-ifelse(out$resAll$whole$t>sd.thr,1/out$resAll$whole$t,1)
		
		## longitudinal calculations
		if(within==TRUE){
			out$resAll$whole$id<-rowMeans(out$res$whole$id,na.rm=TRUE)
			out$resAll$whole$id[which(is.na(out$resAll$whole$id))]<-0 ## NAs converted to 0
			out$resAll$whole$id.t<-rowMeans(out$res$whole$id.t,na.rm=TRUE)
			out$resAll$whole$id.t[which(is.na(out$resAll$whole$id.t))]<-0 ## NAs converted to 0
			out$resAll$within<-list()
			out$resAll$within$t<-rowMeans(out$res$within$t,na.rm=TRUE)
			out$resAll$within$t[which(is.na(out$resAll$within$t))]<-0 ## NAs converted to 0
			out$resAll$within$id<-rowMeans(out$res$within$id,na.rm=TRUE)
			out$resAll$within$id[which(is.na(out$resAll$within$id))]<-0 ## NAs converted to 0
			out$resAll$within$id.t<-rowMeans(out$res$within$id.t,na.rm=TRUE)
			out$resAll$within$id.t[which(is.na(out$resAll$within$id.t))]<-0 ## NAs converted to 0
			## multiplies all weights together 
			out$w[,]<-ifelse(out$resAll$whole$t>sd.thr,1/out$resAll$whole$t,1)*ifelse(out$resAll$whole$id.t>sd.thr,1/out$resAll$whole$id.t,1)*ifelse(out$resAll$within$t>sd.thr,1/out$resAll$within$t,1)*ifelse(out$resAll$within$id.t>sd.thr,1/out$resAll$within$id.t,1)
			out$w.bin.id<-ifelse(out$resAll$whole$id>sd.thr,0,1)*ifelse(out$resAll$within$id>sd.thr,0,1)
			out$w.ind.id<-which(out$w.bin.id==1)
		}
		
		## outlier indices to be excluded for second within-region outlier calculation
		out$w.bin<-ifelse(out$w[,1]==1,1,0)
		out$w.ind<-which(out$w.bin==1)
		
		## saves 1st iteration calculations into nested list
		out$allY<-out

		## recalculate for each y column, excluding 2sd outliers across all regions in regression estimation
		for(i in 1:dim(Y)[2]){

			## whole, t
			out$res$whole$t[,i]<-resid(lm(f),weights=out$w.bin)^2
			## normalizes squared residuals to sd=1 (excluding 2sd outliers), for easier calculation and comparison across regions
			out$res$whole$t[,i]<-out$res$whole$t[,i]/sd(out$res$whole$t[out$w.ind,i],na.rm=TRUE)
			
			## longitudinal calculations
			if(within==TRUE){
			
				for(j in 1:length(unique(id))){
					ind<-which(id==unique(id)[j])
					## across, id
					out$res$whole$id[j,i]<-mean(out$res$whole$t[ind,i],na.rm=TRUE)
					out$res$whole$id.t[ind,i]<-out$res$whole$id[j,i]
					## within, t
					if(sum(out$w.bin[ind])<2) out$res$within$t[ind,i]<-NA else if(sum(out$w.bin[ind])==2) out$res$within$t[ind,i]<-resid(lm(as.formula("Y[ind,i]~1")),weights=out$w.bin[ind])^2 else out$res$within$t[ind,i]<-resid(lm(f.w),weights=out$w.bin[ind])^2
					## within, id
					out$res$within$id[j,i]<-mean(out$res$within$t[ind,i],na.rm=TRUE)
					out$res$within$id.t[ind,i]<-out$res$within$id[j,i]
				}
				
				## normalizes squared residuals to sd=1, for easier calculation and comparison across regions
				out$res$whole$id[,i]<-out$res$whole$id[,i]/sd(out$res$whole$id[,i],na.rm=TRUE)
				out$res$whole$id.t[,i]<-out$res$whole$id.t[,i]/sd(out$res$whole$id[,i],na.rm=TRUE)
				out$res$within$t[,i]<-out$res$within$t[,i]/sd(out$res$within$t[,i],na.rm=TRUE)
				out$res$within$id[,i]<-out$res$within$id[,i]/sd(out$res$within$id[,i],na.rm=TRUE)
				out$res$within$id.t[,i]<-out$res$within$id.t[,i]/sd(out$res$within$id[,i],na.rm=TRUE)
			}
		}
	}
	
	## loop through regions again, this time constructing weighting matrix
	for(i in 1:dim(Y)[2]){

		## removing NAs if present
		out$res$whole$t[which(is.na(out$res$whole$t))]<-0 ## NAs converted to 0
		if(within==TRUE){
			out$res$whole$id[which(is.na(out$res$whole$id))]<-0 ## NAs converted to 0
			out$res$whole$id.t[which(is.na(out$res$whole$id.t))]<-0 ## NAs converted to 0
			out$res$within$t[which(is.na(out$res$within$t))]<-0 ## NAs converted to 0
			out$res$within$id[which(is.na(out$res$within$id))]<-0 ## NAs converted to 0
			out$res$within$id.t[which(is.na(out$res$within$id.t))]<-0 ## NAs converted to 0
		
		## calculation of weights
			out$w[,i]<-out$w[,i]*ifelse(out$res$whole$t[,i]>sd.thr,1/out$res$whole$t[,i],1)*ifelse(out$res$whole$id.t[,i]>sd.thr,1/out$res$whole$id.t[,i],1)*ifelse(out$res$within$t[,i]>sd.thr,1/out$res$within$t[,i],1)*ifelse(out$res$within$id.t[,i]>sd.thr,1/out$res$within$id.t[,i],1)
		}else{
			out$w[,i]<-out$w[,i]*ifelse(out$res$whole$t[,i]>sd.thr,1/out$res$whole$t[,i],1)
		}	
	}
	out
}
