lmer.est<-function(f=cbind(c("","null","")),r=setup$ranef,y=setup$resp,data=data.){
## format formula strings and estimate models (mixed)
## details:
	## f (required) - fixed effects
	## r (required) - random effects
	## y (required) - name of response variable
	## data (required) - data frame with response variable, demographics variables and weights (w)
	na<-numeric(0)
	## fixef
	f.<-pred.str(f[,1])
	if(nchar(f[1,1])>0) na<-c(na,which(is.na(data[,which(names(data)==f[1,1])])))
	if(dim(f)[2]>1){
		for(i in 2:dim(f)[2]){
			f.<-paste(f.,"*",pred.str(f[,i]),sep="")
			na<-c(na,which(is.na(data[,which(names(data)==f[1,i])])))
		}
	}
	## ranef
	r.<-paste("(",pred.str(r[1:3,1]),"|",r[4,1],")",sep="")
	if(dim(r)[2]>1) for(i in 2:dim(r)[2]) r.<-paste(r.,"+","(",pred.str(r[1:3,i]),"|",r[4,i],")",sep="")
	## formula
	formula<-paste(y,"~",f.,"+",r.,sep="")
	if(length(na)>0) data=data[-unique(na),]
	list(fixef=f,formula=formula,na=unique(na),fit=lmer(as.formula(formula),data,REML=FALSE,weights=w))
}

