## llr
filename=paste(setup$path,"llr",sep="/")
sink(filename)
cat("",setup$ynames,"\n")
sink()
write.table(llr,col.names=FALSE,row.names=labels,file=filename,append=TRUE)

## df
filename=paste(setup$path,"llr.df",sep="/")
write.table(llr.df,col.names=FALSE,row.names=labels,file=filename,append=TRUE)

## pvals
filename=paste(setup$path,"llr.p",sep="/")
sink(filename)
cat("",setup$ynames,"\n")
sink()
write.table(llr.p,col.names=FALSE,row.names=labels,file=filename,append=TRUE)

## significance stars
filename=paste(setup$path,"llr.p.star",sep="/")
sink(filename)
cat("",setup$ynames,"\n")
sink()
write.table(
	ifelse(llr.p<.001,"***",ifelse(llr.p<.01,"**",ifelse(llr.p<.05,"*",""))),
	col.names=FALSE,row.names=labels,file=filename,append=TRUE
)
