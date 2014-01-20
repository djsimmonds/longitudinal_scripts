## functions 
model.est<-function(
	f=cbind(c("","null","")),
	x=X,
	y=Y,
	mixed=setup$mixed,
	r=setup$ranef,
	analytics.first=TRUE ## if true, only runs analytics on first Y variable; otherwise, runs on all
){

## format formula strings and estimate models (regular)
## details:
	## f (required) - fixed effects
		## format (variable, type(null,lin,ns,poly,etc...), extra parameters)
	## x (required) - data frame containing demographics and other non-response variables
	## y (required) - response variable matrix
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

	## indexing NAs in response variables, combine all NAs into unique vector and exclude from data frame and y variable
	na<-unique(c(na,which(is.na(rowMeans(as.matrix(y))))))
	if(length(na)>0){
		d<-d[-na,]
		if(is.null(dim(y))) y<-as.matrix(y[-na]) else y<-as.matrix(y[-na,])
	}

	## add response to data frame
	if(is.null(d)) d=data.frame(y=y[,1]) else d<-cbind(d,y=y[,1])

	## calculate model fit
	if(mixed==TRUE) m.call<-paste("lmer(",formula,",data=d)",sep="") else m.call<-paste("lm(",formula,",data=d)",sep="")
	fit<-eval(parse(text=m.call))

	## loop through responses and calculate influence parameters, exclude subjects with cook's distance >4/n for any response
	## NOTE: REST OF SCRIPT IS FLEXIBLE, BUT THIS IS RIGHT NOW ONLY SET UP FOR A SINGLE RANDOM EFFECT, NAMED "id"
	if(analytics.first==TRUE){
		m.est<-influence(fit,"id")
		subj.exc.ind<-which(cooks.distance(m.est)>4/length(unique(d$id)))
	}else{
		subj.exc.ind<-numeric(0)
		for(i in 1:dim(y)[2]){
			fit<-refit(fit,y[,i])
			m.est<-influence(fit,"id")
			subj.exc.ind<-c(subj.exc.ind,which(cooks.distance(m.est)>4/length(unique(d$id))))
		}
	}
	subj.exc<-unique(d$id)[unique(subj.exc.ind)]
	scan.exc<-which(d$id %in% subj.exc)
	if(length(scan.exc)>0) d<-d[-scan.exc,]
	fit<-eval(parse(text=m.call))

	## index of all excluded scans (NA and influence outliers)
	if(length(na)>0) all.exc<-unique(c(na,(1:dim(x)[1])[-na][scan.exc])) else all.exc<-(1:dim(x)[1])[scan.exc]

	## returns list with fixed effects, formula, na indices, model matrix (X) and model fit
	list(
		fixef=f,
		formula=formula,
		na=na,
		subj.exc=subj.exc,
		scan.exc=scan.exc,
		all.exc=all.exc,
		data=d,
		fit=fit
	)
}
