subset.index<-function(subset=setup$subset,demo=demo.){
## function to select subset of subjects based on variables in demographics frame
## details:
	## subset (required) - list(variable name, operator(e,ne,g,l,ge,le), value)
	## intersection across all items in list
	inc<-numeric()
	for(i in 1:length(subset)){
		if(subset=="all"){
			inc<-1:dim(demo)[1]
			break
		}
		ind<-which(names(demo)==subset[[i]][[1]])
		inc.<-switch(subset[[i]][[2]],
			e=which(demo[,ind]==subset[[i]][[3]]),
			ne=which(demo[,ind]!=subset[[i]][[3]]),
			g=which(demo[,ind]>subset[[i]][[3]]),
			l=which(demo[,ind]<subset[[i]][[3]]),
			ge=which(demo[,ind]>=subset[[i]][[3]]),
			le=which(demo[,ind]<=subset[[i]][[3]])
		)
		if(i==1) inc<-inc.
		inc<-intersect(inc,inc.)
	}
	inc
}
