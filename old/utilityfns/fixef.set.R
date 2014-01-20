fixef.set<-function(X=models$X[m,],fixef=setup$fixef,numeric=setup$numeric){
## formats fixed effects for lmer.est()
## details:
	## X (required) - list of numbers corresponding to level and variable of fixed effect
	## fixef (required) - list from setup file containing all fixed effects
	vars<-which(!is.na(X))
	f<-matrix("",3,length(vars))
	if(numeric==TRUE) cat.ind<-3 else cat.ind<-2
	for(i in 1:length(vars)){
		if(vars[i]<cat.ind){
			f[,i]<-fixef[[vars[i]]][,X[vars[i]]]
		}else{
			f[,i]<-fixef[[3]][[vars[i]-cat.ind+1]][,X[vars[i]]]
		}
	}
	f
}

