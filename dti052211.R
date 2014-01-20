## required libraries
#library(lme4) ## for mixed models
library(splines) ## for splines
#library(AICcmodavg) ## for predicting from mixed models

## write results to nifti
library(Rniftilib)
mask<-nifti.image.read(paste(paths$data,"mean_FA_skeleton_mask_2mm_thr02_bin.nii.gz",sep="/")
load(paste(paths$study,"tbss","vox.Rframe",sep="/"))

## functions
lm.est<-function(f=cbind(c("","null","")),x=demo.,y=data.[,1],w=w.[,1]){
## format formula strings and estimate models (regular)
## details:
	## f (required) - fixed effects
	## y (required) - name of response variable
	## data (required) - data frame with response variable, demographics variables and weights (w)
	## need to update comments
	d.<-data.frame(y,w)
	na<-numeric(0)
	## fixef
	f.<-pred.str(f[,1])
	col<-which(names(x)==f[1,1])
	if(nchar(f[1,1])>0){
		na<-c(na,which(is.na(x[,col])))
		d.<-cbind(d.,x[,col])
		colnames(d.)[dim(d.)[2]]<-f[1,1]
	}
	if(dim(f)[2]>1){
		for(i in 2:dim(f)[2]){
			f.<-paste(f.,"*",pred.str(f[,i]),sep="")
			col<-which(names(x)==f[1,i])
			na<-c(na,which(is.na(x[,col])))## null vs. linear effect of age	
			d.<-cbind(d.,x[,col])
			colnames(d.)[dim(d.)[2]]<-f[1,1]
		}
	}
	## formula
	formula<-paste("y~",f.,sep="")
	if(length(na)>0) d.<-d.[-unique(na),]
	fit<-lm(as.formula(formula),data=d.,weights=w)
	X.<-fit$model
	X<-matrix(NA,dim(X.)[1]+length(unique(na)),dim(X.)[2])
	if(length(na)>0) X[-unique(na),]<-X. else X<-X.
	if(dim(X)[2]>2) X=as.matrix(cbind(1,X[,2:dim(X)[2]-1])) else X=as.matrix(numeric(dim(X)[1])+1)
	list(fixef=f,formula=formula,na=unique(na),X=X)
}


#data<-get(setup$data)
#demo.<-get(setup$demo)
demo.$age<-demo.$age-8 ## converting to avoid singularity for correlation of random intercepts/slopes

## models
models<-model.setup()
sink(paste(setup$path,"models.txt",sep="/"))
cat(date(),"\tMODEL SETUP\n\n")
print(models)
sink()

## logfile
log.txt=paste(setup$path,"log.txt",sep="/")
sink(log.txt)
cat(date(),"\tLONGITUDINAL ANALYSIS\n\n")

models$model<-list()
models$model[[1]]<-lm.est()
#models$deriv<-list()

## loop through all models
for(m in 2:dim(models$X)[1]){

	## print progress through m
	cat(date(),"\t\tm =",m,"/",dim(models$X)[1],", vars =",models$X.var[m,],"\n")

	## creates directory for model
	dir.create(paste(setup$path,m,sep="/"))

	## estimate initial model to get design matrix and files to exclude
	models$model[[m]]<-lm.est(fixef.set(models$X[m,]))
	## design matrix
	X<-models$model[[m]]$X
	## scans to be excluded (NA)
	exc<-models$model[[m]]$na

	pt<-proc.time()[[3]]
	ll.test<-sapply(
		1:length(models$con[[m]]),
		function(c){
			## print progress through m
			cat(date(),"\t\tc =",c,"/",length(models$con[[m]]),", vars =",models$X.var[models$con[[m]][c],],"\n")
			cat(models$model[[m]]$formula,"vs",models$model[[models$con[[m]][c]]]$formula,"\n")
			X.n<-models$model[[models$con[[m]][c]]]$X
			exc.<-union(exc,models$model[[models$con[[m]][c]]]$na)
			sapply(
				1:dim(data.)[2],
				function(y) if(length(exc.>0)) llr.q(data.[-exc.,y],w.[-exc.,y],X[-exc.,],X.n[-exc,]) else llr.q(data.[,y],w.[,y],X,X.n)
			)
		}
	)
	save(ll.test,file=paste(setup$path,m,"ll.p.Rdata",sep="/"))
	write.table(ll.test,file=paste(setup$path,m,"ll.p.txt",sep="/"))
	for(c in 1:dim(ll.test)[2]){
		for(i in 1:dim(vox)[1]) mask[vox$i[i],vox$j[i],vox$k[i]]<-ll.test[i,c]
		nifti.set.filenames(mask,paste(setup$path,paste(m,c,"nii.gz",sep="."),sep="/"))
		nifti.image.write(mask)
	}
	cat(date(),"\t\tTime elapsed:",proc.time()[[3]]-pt,"\n")

	## derivatives analysis
	#if(c==1 & deriv.check(models$X[m,],models$ll.test[[m]][[c]]$p)){
	#	cat(date(),"\t\t\t\tderivatives analysis...\n")
	#	if(length(comb.na)==0){
	#		deriv.est()
	#	}else{
	#		if(setup$mixed==TRUE){
	#			m1<-lmer.est(models$model[[m]]$fixef,data=data.[-comb.na,])
	#			m0<-lmer.est(models$model[[models$con[[m]][c]]]$fixef,data=data.[-comb.na,])
	#		}else{
	#			m1<-lm.est(models$model[[m]]$fixef,data=data.[-comb.na,])
	#			m0<-lm.est(models$model[[models$con[[m]][c]]]$fixef,data=data.[-comb.na,])
	#		}
	#		deriv.est(m1,m0,data=data.[-comb.na,])
	#	}
	#	cat(date(),"\t\t\t\tderivatives analysis completed\n")

	## if within is true, calculate these LL tests also (in separate file)
	if(length(models$con.W[[m]])>0){
		pt<-proc.time()[[3]]
		ll.test.W<-sapply(
			1:length(models$con.W[[m]]),
			function(c){
				## print progress through m
				cat(date(),"\t\tc =",c,"/",length(models$con.W[[m]]),", vars =",models$X.var[models$con.W[[m]][c],],"\n")
				cat(models$model[[m]]$formula,"vs",models$model[[models$con.W[[m]][c]]]$formula,"\n")
				X.n<-models$model[[models$con.W[[m]][c]]]$X
				exc.<-union(exc,models$model[[models$con.W[[m]][c]]]$na)
				sapply(
					1:dim(data.)[2],
					function(y) if(length(exc.>0)) llr.q(data.[-exc.,y],w.[-exc.,y],X[-exc.,],X.n[-exc,]) else llr.q(data.[,y],w.[,y],X,X.n)
				)
			}
		)
		save(ll.test.W,file=paste(setup$path,m,"ll.W.p.Rdata",sep="/"))
		write.table(ll.test.W,file=paste(setup$path,m,"ll.W.p.txt",sep="/"))
		for(c in 1:dim(ll.test.W)[2]){
			for(i in 1:dim(vox)[1]) mask[vox$i[i],vox$j[i],vox$k[i]]<-ll.test.W[i,c]
			nifti.set.filenames(mask,paste(setup$path,paste(m,c,"W.nii.gz",sep="."),sep="/"))
			nifti.image.write(mask)
		}
		cat(date(),"\t\tTime elapsed:",proc.time()[[3]]-pt,"\n")
	}
}
sink()


## tract analyses (tbss.dtitk)

## 1
temp<-atlas.ind$L1.type
ind<-numeric(0)
for(t in 2:length(temp)) ind<-union(ind,temp[[t]])

y<-as.numeric(data.[,ind])
w<-as.numeric(w.[,ind])

## null vs. linear effect of age	
X<-matrix(1,length(y),1)
m0<-lm.wfit(X,y,w)
X<-as.matrix(cbind(1,rep(demo.$age,length(ind))))
m1<-lm.wfit(X,y,w)
#> 2*(logLik.(m1)-logLik.(m0))
#[1] 104730.1
#> m1$rank-m0$rank
#[1] 1
#> pchisq(2*(logLik.(m1)-logLik.(m0)),m1$rank-m0$rank)
#[1] 1

## effect of tract type
l1.type<-vox$L1.type....l2d.L1.type.[ind]
X<-as.matrix(cbind(1,rep(demo.$age,length(ind)),rep(l1.type,each=length(demo.$age))))
m2<-lm.wfit(X,y,w)
#> 2*(logLik.(m2)-logLik.(m1))
#[1] 2638.189
#> m2$rank-m1$rank
#[1] 1
#> pchisq(2*(logLik.(m2)-logLik.(m1)),m2$rank-m1$rank)
#[1] 1

## laterality
lat<-sign(vox$x[ind])
X<-as.matrix(cbind(1,rep(demo.$age,length(ind)),rep(l1.type,each=length(demo.$age)),rep(lat,each=length(demo.$age))))
m3<-lm.wfit(X,y,w)
#> 2*(logLik.(m3)-logLik.(m2))
#[1] 660.1963
#> m3$rank-m2$rank
#[1] 1
#> pchisq(2*(logLik.(m3)-logLik.(m2)),m3$rank-m2$rank)
#[1] 1

## within tract type
for(t in 2:length(temp)){
	print(names(temp)[t])
	ind<-temp[[t]]
	y<-as.numeric(data.[,ind])
	w<-as.numeric(w.[,ind])
	## null vs. linear effect of age	
	X<-matrix(1,length(y),1)
	m0<-lm.wfit(X,y,w)
	X<-as.matrix(cbind(1,rep(demo.$age,length(ind))))
	m1<-lm.wfit(X,y,w)
	print(2*(logLik.(m1)-logLik.(m0)))
	print(m1$rank-m0$rank)
	print(pchisq(2*(logLik.(m1)-logLik.(m0)),m1$rank-m0$rank))

	## effect of tract over tract type
	l1.tract<-vox$L1.tract....l2d.L1.tract.[ind]
	X<-as.matrix(cbind(1,rep(demo.$age,length(ind)),rep(l1.tract,each=length(demo.$age))))
	m2<-lm.wfit(X,y,w)
	print(2*(logLik.(m2)-logLik.(m1)))
	print(m2$rank-m1$rank)
	print(pchisq(2*(logLik.(m2)-logLik.(m1)),m2$rank-m1$rank))

	lat<-sign(vox$x[ind])
	X<-as.matrix(cbind(1,rep(demo.$age,length(ind)),rep(l1.tract,each=length(demo.$age)),rep(lat,each=length(demo.$age))))
	m3<-lm.wfit(X,y,w)
	print(2*(logLik.(m3)-logLik.(m2)))
	print(m3$rank-m2$rank)
	print(pchisq(2*(logLik.(m3)-logLik.(m2)),m3$rank-m2$rank))
}
#[1] "P"

#[1] 46293.31
#[1] 1
#[1] 1

#[1] 184.0063
#[1] 1
#[1] 1

#[1] "Ce"

#[1] 6824.062
#[1] 1
#[1] 1

#[1] 3.869281
#[1] 1
#[1] 0.9508224

#[1] "AL"

#[1] 5869.941
#[1] 1
#[1] 1

#[1] 1341.749
#[1] 1
#[1] 1

#[1] "A"

#[1] 25374.26
#[1] 1
#[1] 1

#[1] 52.79605
#[1] 1
#[1] 1

#[1] "Ca"

#[1] 21100.72
#[1] 1
#[1] 1

#[1] 4416.161
#[1] 1
#[1] 1

## tbss

## 1
temp<-atlas.ind$L1.type
ind<-numeric(0)
for(t in 2:length(temp)) ind<-union(ind,temp[[t]])

y<-as.numeric(data.[,ind])
w<-as.numeric(w.[,ind])

## null vs. linear effect of age	
X<-matrix(1,length(y),1)
m0<-lm.wfit(X,y,w)
X<-as.matrix(cbind(1,rep(demo.$age,length(ind))))
m1<-lm.wfit(X,y,w)
#> 2*(logLik.(m1)-logLik.(m0))
#[1] 84111.25
#> m1$rank-m0$rank
#[1] 1
#> pchisq(2*(logLik.(m1)-logLik.(m0)),m1$rank-m0$rank)
#[1] 1

## effect of tract type
l1.type<-vox$L1.type....l2d.L1.type.[ind]
X<-as.matrix(cbind(1,rep(demo.$age,length(ind)),rep(l1.type,each=length(demo.$age))))
m2<-lm.wfit(X,y,w)
#> 2*(logLik.(m2)-logLik.(m1))
#[1] 3954.695
#> m2$rank-m1$rank
#[1] 1
#> pchisq(2*(logLik.(m2)-logLik.(m1)),m2$rank-m1$rank)
#[1] 1

## laterality
lat<-sign(vox$x[ind])
X<-as.matrix(cbind(1,rep(demo.$age,length(ind)),rep(l1.type,each=length(demo.$age)),rep(lat,each=length(demo.$age))))
m3<-lm.wfit(X,y,w)
#> 2*(logLik.(m3)-logLik.(m2))
#[1] 128.3055
#> m3$rank-m2$rank
#[1] 1
#> pchisq(2*(logLik.(m3)-logLik.(m2)),m3$rank-m2$rank)
#[1] 1

## within tract type
for(t in 2:length(temp)){
	print(names(temp)[t])
	ind<-temp[[t]]
	y<-as.numeric(data.[,ind])
	w<-as.numeric(w.[,ind])
	## null vs. linear effect of age	
	X<-matrix(1,length(y),1)
	m0<-lm.wfit(X,y,w)
	X<-as.matrix(cbind(1,rep(demo.$age,length(ind))))
	m1<-lm.wfit(X,y,w)
	print(2*(logLik.(m1)-logLik.(m0)))
	print(m1$rank-m0$rank)
	print(pchisq(2*(logLik.(m1)-logLik.(m0)),m1$rank-m0$rank))

	## effect of tract over tract type
	l1.tract<-vox$L1.tract....l2d.L1.tract.[ind]
	X<-as.matrix(cbind(1,rep(demo.$age,length(ind)),rep(l1.tract,each=length(demo.$age))))
	m2<-lm.wfit(X,y,w)
	print(2*(logLik.(m2)-logLik.(m1)))
	print(m2$rank-m1$rank)
	print(pchisq(2*(logLik.(m2)-logLik.(m1)),m2$rank-m1$rank))

	lat<-sign(vox$x[ind])
	X<-as.matrix(cbind(1,rep(demo.$age,length(ind)),rep(l1.tract,each=length(demo.$age)),rep(lat,each=length(demo.$age))))
	m3<-lm.wfit(X,y,w)
	print(2*(logLik.(m3)-logLik.(m2)))
	print(m3$rank-m2$rank)
	print(pchisq(2*(logLik.(m3)-logLik.(m2)),m3$rank-m2$rank))
}
[1] "P"
[1] 32329.98
[1] 1
[1] 1
[1] 278.4300
[1] 1
[1] 1
[1] 326.6784
[1] 1
[1] 1

[1] "Ce"
[1] 4364.925
[1] 1
[1] 1
[1] 33.57017
[1] 1
[1] 1
[1] 130.1998
[1] 1
[1] 1

[1] "AL"
[1] 5346.018
[1] 1
[1] 1
[1] 170.7540
[1] 1
[1] 1
[1] 130.7468
[1] 1
[1] 1

[1] "A"
[1] 23108.64
[1] 1
[1] 1
[1] 34.95811
[1] 1
[1] 1
[1] 141.6889
[1] 1
[1] 1

[1] "Ca"
[1] 19710.64
[1] 1
[1] 1
[1] 10149.60
[1] 1
[1] 1
[1] 0.0002295685
[1] 1
[1] 0.01208870

X<-matrix(1,length(y),1)
m0<-lm.wfit(X,y,w)
X<-as.matrix(cbind(1,rep(demo.$age,length(ind))))
m1<-lm.wfit(X,y,w)
#> 2*(logLik.(m1)-logLik.(m0))
#[1] 104730.1
#> m1$rank-m0$rank
#[1] 1
#> pchisq(2*(logLik.(m1)-logLik.(m0)),m1$rank-m0$rank)
#[1] 1

## effect of tract type
l1.type<-vox$L1.type....l2d.L1.type.[ind]
X<-as.matrix(cbind(1,rep(demo.$age,length(ind)),rep(l1.type,each=length(demo.$age))))
m2<-lm.wfit(X,y,w)
#> 2*(logLik.(m2)-logLik.(m1))
#[1] 2638.189
#> m2$rank-m1$rank
#[1] 1
#> pchisq(2*(logLik.(m2)-logLik.(m1)),m2$rank-m1$rank)
#[1] 1

## laterality
lat<-sign(vox$x[ind])
X<-as.matrix(cbind(1,rep(demo.$age,length(ind)),rep(l1.type,each=length(demo.$age)),rep(lat,each=length(demo.$age))))
m3<-lm.wfit(X,y,w)
#> 2*(logLik.(m3)-logLik.(m2))
#[1] 660.1963
#> m3$rank-m2$rank
#[1] 1
#> pchisq(2*(logLik.(m3)-logLik.(m2)),m3$rank-m2$rank)
#[1] 1

## within tract type
for(t in 2:length(temp)){
	print(names(temp)[t])
	ind<-temp[[t]]
	y<-as.numeric(data.[,ind])
	w<-as.numeric(w.[,ind])
	## null vs. linear effect of age	
	X<-matrix(1,length(y),1)
	m0<-lm.wfit(X,y,w)
	X<-as.matrix(cbind(1,rep(demo.$age,length(ind))))
	m1<-lm.wfit(X,y,w)
	print(2*(logLik.(m1)-logLik.(m0)))
	print(m1$rank-m0$rank)
	print(pchisq(2*(logLik.(m1)-logLik.(m0)),m1$rank-m0$rank))

	## effect of tract over tract type
	l1.tract<-vox$L1.tract....l2d.L1.tract.[ind]
	X<-as.matrix(cbind(1,rep(demo.$age,length(ind)),rep(l1.tract,each=length(demo.$age))))
	m2<-lm.wfit(X,y,w)
	print(2*(logLik.(m2)-logLik.(m1)))
	print(m2$rank-m1$rank)
	print(pchisq(2*(logLik.(m2)-logLik.(m1)),m2$rank-m1$rank))

	lat<-sign(vox$x[ind])
	X<-as.matrix(cbind(1,rep(demo.$age,length(ind)),rep(l1.tract,each=length(demo.$age)),rep(lat,each=length(demo.$age))))
	m3<-lm.wfit(X,y,w)
	print(2*(logLik.(m3)-logLik.(m2)))
	print(m3$rank-m2$rank)
	print(pchisq(2*(logLik.(m3)-logLik.(m2)),m3$rank-m2$rank))
}
#[1] "P"

#[1] 46293.31
#[1] 1
#[1] 1

#[1] 184.0063
#[1] 1
#[1] 1

#[1] "Ce"

#[1] 6824.062
#[1] 1
#[1] 1

#[1] 3.869281
#[1] 1
#[1] 0.9508224

#[1] "AL"

#[1] 5869.941
#[1] 1
#[1] 1

#[1] 1341.749
#[1] 1
#[1] 1

#[1] "A"

#[1] 25374.26
#[1] 1
#[1] 1

#[1] 52.79605
#[1] 1
#[1] 1

#[1] "Ca"

#[1] 21100.72
#[1] 1
#[1] 1

#[1] 4416.161
#[1] 1
#[1] 1
