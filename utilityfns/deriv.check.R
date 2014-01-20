## check against setup file to see if derivative analysis should be performed
deriv.check<-function(vars=NA,p,p.thr=setup$deriv$p.thr,inc=setup$deriv$vars){
## details:
	## vars (required) - for each level of fixed effects, indicates which variables are included in model
	## p - p-value of analysis (if you want to skip the derivatives analysis for non-significant log-likelihood tests) (default = 0, none skipped)
	## p.thr - threshold that p must exceed in order for derivatives analysis to proceed
	## inc - all the variables to be included in derivatives analyses (cross checks against this to ensure that analysis should be performed) 

	## need to have main variable for analysis
	if(is.na(vars[1])) return(FALSE)
	## only perform analysis if log-likelihood test is significant
	if(p>p.thr) return(FALSE)
	## all variables should be included in derivatives list in setup file
	for(i in 1:length(vars)) if(!is.na(vars[i]) & length(which(inc[[i]]==vars[i]))==0) return(FALSE)
	## if all conditions are met, do analysis
	TRUE
}
