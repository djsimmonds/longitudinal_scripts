fixef.set<-function(X=models$X[m,],fixef=setup$fixef){
## formats fixed effects for lmer.est()
## details:
	## X (required) - list of numbers corresponding to level and variable of fixed effect
	## fixef (required) - list from setup file containing all fixed effects
	vars<-which(!is.na(X))
	f<-matrix("",3,length(vars))
	for(i in 1:length(vars)) f[,i]<-fixef[[vars[i]]][,X[vars[i]]]
	f
}

