## check against setup file to see if derivative analysis should be performed
deriv.check<-function(vars,setup,p=1){
## details:
	## vars (required) - for each level of fixed effects, indicates which variables are included in model
	## p - p-value of analysis (if you want to skip the derivatives analysis for non-significant log-likelihood tests) (default = 0, none skipped)
	## p.thr - threshold that p must exceed in order for derivatives analysis to proceed
	## inc - all the variables to be included in derivatives analyses (cross checks against this to ensure that analysis should be performed) 

	D<-list()
	if(setup$main==TRUE) D[[1]]<-setup$deriv$vars$main else return(0) ## need main variable for derivatives analysis
	if(setup$numeric==TRUE) D[[length(D)+1]]<-setup$deriv$vars$numeric
	if(length(setup$deriv$vars$categorical)>0) for(i in 1:length(setup$deriv$vars$categorical)) D[[length(D)+1]]<-setup$deriv$vars$categorical[[i]]
	if(length(D)==0) return(0)

	## only perform analysis if main is in model
	if(is.na(vars[1])) return(0)
	## only perform analysis if log-likelihood test is significant
	if(p<1-setup$deriv$p.thr) return(0)
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
