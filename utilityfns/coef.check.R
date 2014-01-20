## check against setup file to see if coefficients analysis should be performed
coef.check<-function(vars,setup){
## details:
	## vars (required) - for each level of fixed effects, indicates which variables are included in model

	D<-list()
	if(setup$main==TRUE) D[[1]]<-setup$coef$vars$main else return(0) ## need main variable for derivatives analysis
	if(setup$numeric==TRUE) D[[length(D)+1]]<-setup$coef$vars$numeric
	if(length(setup$coef$vars$categorical)>0) for(i in 1:length(setup$coef$vars$categorical)) D[[length(D)+1]]<-setup$coef$vars$categorical[[i]]
	if(length(D)==0) return(0)

	## only perform analysis if main is in model
	if(is.na(vars[1])) return(0)
	## all variables should be included in derivatives list in setup file
	for(i in 1:length(vars)) if(!is.na(vars[i]) & length(which(D[[i]]==vars[i]))==0) return(0)
	## if all conditions are met, do analysis
	
	## different types of derivatives analysis and figures
		## 1) y~main
		## 2) y~main*numeric
		## 3) y~main*categorical(1)
		## 4) y~main*numeric*categorical
		## 5) y~main*categorical(>1)
	if(length(which(!is.na(vars)))==1){
		return(1)
	}else if(setup$numeric==TRUE & !is.na(vars[2])){
		if(length(which(!is.na(vars)))==2){
			return(2)
		}else{
			return(4)
		}
	}else{
		if(length(which(!is.na(vars)))==2){
			return(3)
		}else{
			return(5)
		}
	}
}
