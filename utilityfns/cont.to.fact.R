cont.to.fact<-function(vars=setup$cont.to.fact$vars,demo=demo.,cut=setup$cont.to.fact$cut){
## convert continuous variable to factor (returns data frame with continuous variables replaced by ordered factors
	## vars (required) - variables to convert
	## demo (required) - demographics frame where variable is
	## cut (required) - percentiles at which to categorize
	for(v in 1:length(vars)){
		cat(date(),"\tvar =",vars[v],"\n")
		col<-which(names(demo)==vars[v])
		var.new<-numeric(length(demo[,col]))+1
		for(i in 1:length(cut[[v]])){
			cat(date(),"\t\tperc =",round(cut[[v]][i],2),", value =",round(quantile(demo[,col],cut[[v]][i],na.rm=TRUE),2),"\n")
			var.new<-var.new+ifelse(demo[,col]>=quantile(demo[,col],cut[[v]][i],na.rm=TRUE),1,0)
		}
		var.new<-ifelse(is.na(demo[,col]),NA,var.new)
		ord<-character(length(cut[[v]])+1)
		ord[1]<-paste(0,round(cut[[v]][1]*100,0),sep="-")
		ord[length(ord)]<-paste(round(cut[[v]][length(cut[[v]])]*100,0),100,sep="-")	
		if(length(ord)>2) for(i in 1:(length(cut[[v]])-1)) ord[i+1]<-paste(round(cut[[v]][i]*100,0),round(cut[[v]][i+1]*100,0),sep="-")
		ord<-as.factor(ord) ## replaced ordered with factor because was having trouble with prediction
		demo[,col]<-ord[var.new]
	}
	demo
}
