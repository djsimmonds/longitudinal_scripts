lm.est<-function(f=cbind(c("","null","")),y=setup$resp,data=data.){
## format formula strings and estimate models (regular)
## details:
	## f (required) - fixed effects
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
	## formula
	formula<-paste(y,"~",f.,sep="")
	if(length(na)>0) data=data[-unique(na),]
	list(fixef=f,formula=formula,na=unique(na),fit=lm(as.formula(formula),data,weights=w))
}

