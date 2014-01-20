## functions 
model.est2<-function(
	model=models$model[[m]],
	y=Y,
	mixed=setup$mixed,
	r=setup$ranef
){

## format formula strings and estimate models (regular)
## details:
	## f (required) - fixed effects
		## format (variable, type(null,lin,ns,poly,etc...), extra parameters)
	## x (required) - data frame containing demographics and other non-response variables
	## y (required) - response variable matrix
	## mixed - regular regression (lm) or mixed model regression (lmer)? (default=FALSE, for lm)
	## r - random effects if mixed model regression

	d<-model$data
	f<-model$fixef

	if(dim(f)[2]==1){
		f.<-"1"
	}else if(dim(f)[2]==2){
		f.<-paste(pred.str(f[,1]),pred.str(f[,2]),sep="+")
	}else if(dim(f)[2]==3){
		f.<-""
		for(i in 1:dim(f)[2]){
			for(j in 1:dim(f)[2]){
				if(f.=="") sep="" else sep="+"
				if(j>i) f.<-paste(f.,paste(pred.str(f[,i]),pred.str(f[,j]),sep="*"),sep=sep)
			}
		}
	}else if(dim(f)[2]==4){
		f.<-""
		for(i in 1:dim(f)[2]){
			for(j in 1:dim(f)[2]){
				for(k in 1:dim(f)[2]){
					if(f.=="") sep="" else sep="+"
					if(k>j&j>i) f.<-paste(f.,paste(pred.str(f[,i]),pred.str(f[,j]),pred.str(f[,k]),sep="*"),sep=sep)
				}
			}
		}
	}
	## this is as many levels as i need right now, but more can be added later
	
	## ranef
	if(mixed==TRUE){
		r.<-paste("(",pred.str(r[1:3,1]),"|",r[4,1],")",sep="") ## formats string for random effects portion of formula
		if(dim(r)[2]>1){ ## if more than one random effect, add to formula string and data frame
			for(i in 2:dim(r)[2]) r.<-paste(r.,"+","(",pred.str(r[1:3,i]),"|",r[4,i],")",sep="")
		}
			
	## formula string
		formula<-paste("y~",f.,"+",r.,sep="")
	}else{
		formula<-paste("y~",f.,sep="")
	}

	## calculate model fit
	if(mixed==TRUE) m.call<-paste("lmer(",formula,",data=d)",sep="") else m.call<-paste("lm(",formula,",data=d)",sep="")
	list(fit=eval(parse(text=m.call)),formula=formula)
}
