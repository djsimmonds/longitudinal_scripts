## utility functions for analysis written by Dani
utilityfns<-function(file,path=paths$study){
	source(paste(path,"scripts","utilityfns",file,sep="/"))
	cat(date(),"\t-",file,"loaded\n") ## an error check might be useful here
}


utilityfns("outlier.R")
## function for dealing with outliers
#outlier<-function(Y=get(setup$data)[subset,],X=get(setup$demo)$age[subset,],id=get(setup$demo)$id[subset,],within=setup$out$within,allY=setup$out$allY,sd.thr=setup$out$sd.thr)
## details:
	## Y - (required) - matrix of data, dimensions (# of scans, # of y measurements (such as voxels, ROIs, etc...)) 
	## id - id for within-subject (or group) residual calculation (array of length # of scans) if applicable
	## X - variables to be used in regression for calculation of residuals (i.e. age) (each an array of length # of scans)
		## only implemented using linear regression
	## within - should within subject residuals be calculated (TRUE/FALSE) (default=TRUE)
	## allY - should outliers across all y variables be calculated and removed before assessing within each y (TRUE/FALSE) (default=TRUE)
	## sd.thr - standard deviation threshold at which weighting is applied (default=2)
	## returns list with all outlier calculations


utilityfns("model.setup.R")
#model.setup<-function(F=setup$fixef,W=setup$within,D=setup$deriv$vars)
## set up model structure with:
	## total number of models to examine
	## to which models they need to be compared
	## which models to run derivatives analysis for
## details:
	## F (required) - fixed effects of model (indicated in setup file)
	## W (required) - compare main models to each other?
	## D (required) - which models to compute derivatives for


utilityfns("deriv.check.R")
#deriv.check<-function(vars=NA,p=0,p.thr=setup$deriv$p.thr,inc=setup$deriv$vars)
## check against setup file to see if derivative analysis should be performed
## details:
	## vars (required) - for each level of fixed effects, indicates which variables are included in model
	## p - p-value of analysis (if you want to skip the derivatives analysis for non-significant log-likelihood tests) (default = 0, none skipped)
	## p.thr - threshold that p must exceed in order for derivatives analysis to proceed
	## inc - all the variables to be included in derivatives analyses (cross checks against this to ensure that analysis should be performed) 


utilityfns("subset.index.R")
#subset.index<-function(subset=setup$subset)
## function to select subset of subjects based on variables in demographics frame
## intersection across all items in list
## details:
	## subset (required) - list(variable name, operator(e,ne,g,l,ge,le), value)


utilityfns("fixef.set.R")
#fixef.set<-function(X=models$X[m,]),fixef=setup$fixef)
## formats fixed effects for lmer.est()
## details:
	## X (required) - list of numbers corresponding to level and variable of fixed effect
	## fixef (required) - list from setup file containing all fixed effects


utilityfns("pred.str.R")
#pred.str<-function(x)
## formats fixed effect into string appropriate for formula
	## details:
	## x = fixed effect for formatting c(var,type,param)


utilityfns("lm.est.R")
#lm.est<-function(f=cbind(c("","null","")),y=setup$resp,data=data.)
## format formula strings and estimate models (regular)
## details:
	## f (required) - fixed effects
	## y (required) - name of response variable
	## data (required) - data frame with response variable, demographics variables and weights (w)


utilityfns("lmer.est.R")
#lmer.est<-function(f=cbind(c("","null","")),r=setup$ranef,y=setup$resp,data=data.)
## format formula strings and estimate models (mixed)
## details:
	## f (required) - fixed effects
	## r (required) - random effects
	## y (required) - name of response variable
	## data (required) - data frame with response variable, demographics variables and weights (w)


utilityfns("ll.test.R")
#ll.test<-function(m1=models$model[[m]],m0=models$model[[models$con[[m]][c]]],data=data.,mixed=setup$mixed,anova.<-anova.txt)
## calculate LL test (if there are NAs, rerun models and then calculate LL test)
## details:
	## m1 (required) - model of interest
	## m2 (required) - model being tested against (null model)
	## data (required) - data frame if NAs and models need to be rerun
	## mixed (required) - TRUE if mixed model, FALSE if regular regression model


utilityfns("cont.to.fact.R")
#cont.to.fact<-function(vars=setup$cont.to.fact$vars,demo=demo.,cut=setup$cont.to.fact$cut)
## convert continuous variable to factor (returns data frame with continuous variables replaced by ordered factors
	## vars (required) - variables to convert
	## demo (required) - demographics frame where variable is
	## cut (required) - percentiles at which to categorize


utilityfns("deriv.est.R")
#deriv.est<-function(m1=models$model[[m]],m0=models$model[[models$con[[m]][c]]],path=path..,deriv=setup$deriv,mixed=setup$mixed,data=data.)
## derivatives analysis
## details:
	## need to fill in
