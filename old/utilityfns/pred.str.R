pred.str<-function(x){
## formats fixed effect into string appropriate for formula
	## details:
	## x = fixed effect for formatting c(var,type,param)
	if(x[2]=="null") "1" else
	if(x[2]=="lin") x[1] else
	paste(x[2],"(",x[1],",",x[3],")",sep="")
}
