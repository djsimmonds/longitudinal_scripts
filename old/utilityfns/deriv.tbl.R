## pred.grid
filename=paste(path.tables,"pred.grid",sep="/")
sink(filename)
cat("",names(pred.grid),"\n")
sink()
write.table(pred.grid,col.names=FALSE,file=filename,append=TRUE)

## predicted values
filename=paste(path.tables,"pred",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(pred,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## mean of simulated predictions
filename=paste(path.tables,"sim.pred.mean",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(pred.sim.m,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## sd of simulated predictions
filename=paste(path.tables,"sim.pred.sd",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(pred.sim.sd,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## predicted derivative values
filename=paste(path.tables,"pred.d",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(pred.d,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## mean of simulated derivative predictions
filename=paste(path.tables,"sim.pred.d.mean",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(pred.sim.d.m,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## sd of simulated derivative predictions
filename=paste(path.tables,"sim.pred.d.sd",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(pred.sim.d.sd,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## pred.p
filename=paste(path.tables,"pred.p",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(p,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## pred.p significance stars
filename=paste(path.tables,"pred.p.star",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*",""))),
	col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## pred.p sign (-1 for low/1 for high)
filename=paste(path.tables,"pred.p.sign",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(p.sign,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## pred.d.p
filename=paste(path.tables,"pred.d.p",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(p.d,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## pred.d.p significance stars
filename=paste(path.tables,"pred.d.p.star",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(ifelse(p.d<.001,"***",ifelse(p.d<.01,"**",ifelse(p.d<.05,"*",""))),
	col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## pred.d.p sign (-1 for low/1 for high)
filename=paste(path.tables,"pred.d.p.sign",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(p.d.sign,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## pred.pboot
filename=paste(path.tables,"pred.pboot",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(pboot,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## pred.pboot significance stars
filename=paste(path.tables,"pred.pboot.star",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(ifelse(pboot<.001,"***",ifelse(pboot<.01,"**",ifelse(pboot<.05,"*",""))),
	col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## pred.pboot sign (-1 for low/1 for high)
filename=paste(path.tables,"pred.pboot.sign",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(pboot.sign,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## pred.d.pboot
filename=paste(path.tables,"pred.d.pboot",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(pboot.d,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## pred.d.pboot significance stars
filename=paste(path.tables,"pred.d.pboot.star",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(ifelse(pboot.d<.001,"***",ifelse(pboot.d<.01,"**",ifelse(pboot.d<.05,"*",""))),
	col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

## pred.d.p sign (-1 for low/1 for high)
filename=paste(path.tables,"pred.d.pboot.sign",sep="/")
sink(filename)
cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
sink()
write.table(pboot.d.sign,col.names=FALSE,row.names=apply(pred.grid,1,paste,collapse=","),file=filename,append=TRUE)

if(deriv.do>1){

	## predicted.dif values
	filename=paste(path.tables,"pred.dif",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pred.dif,col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## mean of simulated predictions.dif
	filename=paste(path.tables,"sim.pred.dif.mean",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pred.sim.dif.m,col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## sd of simulated predictions.dif
	filename=paste(path.tables,"sim.pred.dif.sd",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pred.sim.dif.sd,col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## predicted.dif derivative values
	filename=paste(path.tables,"pred.dif.d",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pred.dif.d,col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## mean of simulated derivative predictions.dif
	filename=paste(path.tables,"sim.pred.dif.d.mean",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pred.sim.dif.d.m,col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## sd of simulated derivative predictions.dif
	filename=paste(path.tables,"sim.pred.dif.d.sd",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(pred.sim.dif.d.sd,col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## pred.dif.p
	filename=paste(path.tables,"pred.dif.p",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(dif.p,col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## pred.dif.p significance stars
	filename=paste(path.tables,"pred.dif.p.star",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(ifelse(dif.p<.001,"***",ifelse(dif.p<.01,"**",ifelse(dif.p<.05,"*",""))),
		col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## pred.dif.p sign (-1 for low/1 for high)
	filename=paste(path.tables,"pred.dif.p.sign",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(dif.p.sign,col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## pred.dif.d.p
	filename=paste(path.tables,"pred.dif.d.p",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(dif.p.d,col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## pred.dif.d.p significance stars
	filename=paste(path.tables,"pred.dif.d.p.star",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(ifelse(dif.p.d<.001,"***",ifelse(dif.p.d<.01,"**",ifelse(dif.p.d<.05,"*",""))),
		col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## pred.dif.d.p sign (-1 for low/1 for high)
	filename=paste(path.tables,"pred.dif.d.p.sign",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(dif.p.d.sign,col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## pred.dif.pboot
	filename=paste(path.tables,"pred.dif.pboot",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(dif.pboot,col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## pred.dif.pboot significance stars
	filename=paste(path.tables,"pred.dif.pboot.star",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(ifelse(dif.pboot<.001,"***",ifelse(dif.pboot<.01,"**",ifelse(dif.pboot<.05,"*",""))),
		col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## pred.dif.pboot sign (-1 for low/1 for high)
	filename=paste(path.tables,"pred.dif.pboot.sign",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(dif.pboot.sign,col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## pred.dif.d.pboot
	filename=paste(path.tables,"pred.dif.d.pboot",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(dif.pboot.d,col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## pred.dif.d.pboot significance stars
	filename=paste(path.tables,"pred.dif.d.pboot.star",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(ifelse(dif.pboot.d<.001,"***",ifelse(dif.pboot.d<.01,"**",ifelse(dif.pboot.d<.05,"*",""))),col.names=FALSE,row.names=range,file=filename,append=TRUE)

	## pred.dif.d.p sign (-1 for low/1 for high)
	filename=paste(path.tables,"pred.dif.d.pboot.sign",sep="/")
	sink(filename)
	cat(paste(names(pred.grid),collapse=","),setup.$ynames,"\n")
	sink()
	write.table(dif.pboot.d.sign,col.names=FALSE,row.names=range,file=filename,append=TRUE)

}
