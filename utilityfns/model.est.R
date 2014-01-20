## functions 
model.est<-function(f=cbind(c("","null","")),x=demo.,y=data.[,1],w=Wts[,1],mixed=setup$mixed,r=setup$ranef){
## format formula strings and estimate models (regular)
## details:
	## f (required) - fixed effects
		## format (variable, type(null,lin,ns,poly,etc...), extra parameters)
	## x (required) - data frame containing demographics and other non-response variables
	## y (required) - response variable vector (default to first column of data matrix)
	## w - weights vector (default to first column of weights matrix)
	## mixed - regular regression (lm) or mixed model regression (lmer)? (default=FALSE, for lm)
	## r - random effects if mixed model regression

	## indices of NAs for removal in analysis
	na<-numeric(0)
	
	## fixef
	f.<-pred.str(f[,1]) ## for formula string
	col<-which(names(x)==f[1,1]) ## column index for variable
	if(nchar(f[1,1])>0){ ## if there is a variable, add to data frame
		na<-c(na,which(is.na(x[,col])))
		d<-data.frame(x[,col])
		colnames(d)[dim(d)[2]]<-f[1,1]
	}else{
		d<-NULL
	}
	if(dim(f)[2]>1){ ## if more than 1 fixed effect, add to formula string and data frame
		for(i in 2:dim(f)[2]){
			f.<-paste(f.,"*",pred.str(f[,i]),sep="")
			col<-which(names(x)==f[1,i])
			na<-c(na,which(is.na(x[,col])))
			d<-cbind(d,x[,col])
			colnames(d)[dim(d)[2]]<-f[1,i]
		}
	}
	
	## ranef
	if(mixed==TRUE){
		r.<-paste("(",pred.str(r[1:3,1]),"|",r[4,1],")",sep="") ## formats string for random effects portion of formula
		if(length(which(names(d)==r[1,1]))==0){ ## checks if variable already in data frame, adds if not
			col<-which(names(x)==r[1,1])
			na<-c(na,which(is.na(x[,col])))
			if(is.null(d)) d=data.frame(x[,col]) else d<-cbind(d,x[,col])
			colnames(d)[dim(d)[2]]<-r[1,1]		
		}
		if(length(which(names(d)==r[4,1]))==0){
			col<-which(names(x)==r[4,1])
			na<-c(na,which(is.na(x[,col])))
			d<-cbind(d,x[,col])
			colnames(d)[dim(d)[2]]<-r[4,1]		
		}		
		if(dim(r)[2]>1){ ## if more than one random effect, add to formula string and data frame
			for(i in 2:dim(r)[2]){
				r.<-paste(r.,"+","(",pred.str(r[1:3,i]),"|",r[4,i],")",sep="")
				if(length(which(names(d)==r[1,i]))==0){
					col<-which(names(x)==r[1,i])
					na<-c(na,which(is.na(x[,col])))
					d<-cbind(d,x[,col])
					colnames(d)[dim(d)[2]]<-r[1,i]		
				}
				if(length(which(names(d)==r[4,i]))==0){
					col<-which(names(x)==r[4,i])
					na<-c(na,which(is.na(x[,col])))
					d<-cbind(d,x[,col])
					colnames(d)[dim(d)[2]]<-r[4,i]		
				}		
			}
		}
			
	## formula string
		formula<-paste("y~",f.,"+",r.,sep="")
	}else{
		formula<-paste("y~",f.,sep="")
	}

	## adding response and weights to data frame, also indexing NAs
	if(is.null(d)) d=data.frame(y) else d<-cbind(d,y)
	na<-c(na,which(is.na(y)))
	if(!is.null(w)){
		d<-cbind(d,w)
		na<-c(na,which(is.na(w)))
	}
	
	## exclude NAs from data frame
	if(length(na)>0) d<-d[-unique(na),]

	## calculate model fit
	if(mixed==TRUE){
		m.call<-paste("lmer(",formula,",data=d",sep="")
		if(!is.null(w)) m.call<-paste(m.call,",weights=w)",sep="") else m.call<-paste(m.call,")",sep="")
		fit<-eval(parse(text=m.call))
	}else{
		m.call<-paste("lm(",formula,",data=d",sep="")
		if(!is.null(w)) m.call<-paste(m.call,",weights=w)",sep="") else m.call<-paste(m.call,")",sep="")
		fit<-eval(parse(text=m.call))
	}
	
	## returns list with fixed effects, formula, na indices, model matrix (X) and model fit
	list(fixef=f,formula=formula,na=unique(na),data=d,fit=fit)
}
