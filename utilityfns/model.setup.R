model.setup<-function(setup){
## set up model structure with:
	## total number of models to examine
	## to which models they need to be compared
	## which models to run derivatives analysis for
## details:
	## F (required) - fixed effects of model (indicated in setup file)
	## D (required) - which models to compute derivatives for

	F<-list()
	if(setup$main==TRUE) F[[1]]<-setup$fixef$main
	if(setup$numeric==TRUE) F[[length(F)+1]]<-setup$fixef$numeric
	if(length(setup$fixef$categorical)>0) for(i in 1:length(setup$fixef$categorical)) F[[length(F)+1]]<-setup$fixef$categorical[[i]]
	if(length(F)==0) stop("ERROR: NO VARIABLES IN MODEL")

	## initialize models structure with null model as first row
	models<-list(
		L=numeric(1),
		M=numeric(1),
		X=array(NA,c(1,length(F))),
		X.var=array("",c(1,length(F))),
		deriv.do=NA,
		coef.do=NA,
		con=list(numeric(0))
	)
	m<-2 ## model index (rows of models)
	M<-1 ## model group
	## first loops over levels (# predictor variable groups)
	for(L in 1:length(F)){
		grid<-list()
		for(l in 1:L) grid[[l]]<-c(1:length(F))
		grid<-expand.grid(grid)
		grid.ind<-numeric(0)
		if(L>1) for(l in 2:L) grid.ind<-union(grid.ind,which(grid[,l]<=grid[,l-1]))
		if(length(grid.ind)>0) grid<-grid[-grid.ind,]
		## next loop over individual predictor variables
		for(i in 1:dim(grid)[1]){
			grid2<-list()
			for(j in 1:dim(grid)[2]){
				grid2[[j]]<-1:dim(F[[grid[i,j]]])[2]
			}
			grid2<-expand.grid(grid2)
			## set up models
			for(j in 1:dim(grid2)[1]){
				models$L[m]<-L
				models$M[m]<-M
				models$X<-rbind(models$X,models$X[1,])
				models$X.var<-rbind(models$X.var,models$X.var[1,])
				for(k in 1:dim(grid2)[2]){
					models$X[m,grid[i,k]]<-grid2[j,k]
					models$X.var[m,grid[i,k]]<-pred.str(F[[grid[i,k]]][,grid2[j,k]])
				}
				models$deriv.do[m]<-deriv.check(models$X[m,],setup)
				models$coef.do[m]<-coef.check(models$X[m,],setup)
				## set up contrasts (which rows to compare model to)
				if(L==1){
					models$con[[m]]<-1
				}else{
					temp<-array(0,c(m-1,length(F)))
					for(k in 1:length(F)) if(!is.na(models$X[m,k])) temp[,k]<-(models$L[1:(m-1)]==L-1)*(models$X[1:(m-1),k]==models$X[m,k])
					models$con[[m]]<-which(rowSums(temp,na.rm=TRUE)==max(rowSums(temp,na.rm=TRUE)))
				}
				m<-m+1
			}
			M<-M+1
		}
	}
	save(models,file=paste(setup$path,"models",sep="/"))
	models
}
