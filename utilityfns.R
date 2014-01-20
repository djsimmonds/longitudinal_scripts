## utility functions for analysis written by Dani
utilityfns<-function(file,path=setup$path.scripts){
	source(paste(path,"utilityfns",file,sep="/"))
	cat(date(),"\t-",file,"loaded\n") ## an error check might be useful here
}

utilityfns("coef.check.R")
utilityfns("deriv.check.R")
utilityfns("model.setup.R")
utilityfns("pred.str.R")
utilityfns("fixef.set.R")
utilityfns("model.est.R")
utilityfns("ll.test.R")
utilityfns("coef.est.R")
utilityfns("deriv.est.R")
#utilityfns("ll.tbl.R")
#utilityfns("deriv.tbl.R")

