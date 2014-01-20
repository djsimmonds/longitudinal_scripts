#####################################
# longitudinal DTI study            #
#   ~	~	~	~	~	~	~	~                 #
# written by:       when:           #
#   Dani Simmonds     12/10 - 5/11  #
#####################################

########################################
# this file is important for setting up demographics data as a precursor to analysis
#
# 1) extract subject info (ID, age, sex, puberty, IQ)
# 2) match up subject info to scan info:
#	a) array of scans (n = 425)
#	b) array of unique, included subjects (n = ?)
# 3) process variables to create factors and other precursors to analysis
########################################

## paths
paths<-list()
paths$study<-"/home/danisimmonds/Dani/dti_0511"
paths$data<-"/mnt/Schwarzenagger/Governator/DTI_STUDY/dtitk/icbm_dti81/tbss/origdata"
paths$demo<-paste(paths$study,"demographics",sep="/")

## import demographic/behavioral data from LNCD database (thanks, david)
demo<-list()

demo$unique<-read.table(paste(paths$demo,"dtiInfo_unique.csv",sep="/"),sep=",",header=TRUE,fill=TRUE)
#names(demo$unique)
#[1] "LunaID"     "DateOfBirth"  "SexID" (1=male)       "Exclude"      "notesExclude"

demo$puberty<-read.table(paste(paths$demo,"dtiInfo_puberty.csv",sep="/"),sep=",",header=TRUE,fill=TRUE)
#names(demo$puberty)
# [1] "VisitID"             "LunaID"              "SexID"              
# [4] "DateOfBirth"         "VisitDate"           "Age"                
# [7] "pps1growth"          "pps2bodyhair"        "pps3skin"           
#[10] "Fpps4breast"         "FppsMenstr"          "FppsMenstrAgeYR"    
#[13] "FppsMenstrAgeMO"     "Mpps4voice"          "Mpps5hairface"      
#[16] "pps6devt"            "heightft"            "heightin"           
#[19] "weightlbs"           "pubcomplete"         "agepubcompYRS"      
#[22] "agepubcompMO"        "pps1to5"             "Tanner1hair"        
#[25] "Tanner2TSP_orBreast" "Tanner3m"            "tsr"                
#[28] "Mpubic"              "Mgenital"            "Fage_Menstration"   
#[31] "Fage_Breast"         "Fpubic"              "Notes"   

demo$wasi<-read.table(paste(paths$demo,"dtiInfo_wasi.csv",sep="/"),sep=",",header=TRUE,fill=TRUE)
#names(demo$wasi)
# [1] "VisitID"     "LunaID"      "SexID"       "DateOfBirth" "VisitDate"  
# [6] "Age"         "wvocrw"      "wvocT"       "wbdrw"       "wbdT"       
#[11] "wsimrw"      "wsimT"       "wmrrw"       "wmrT"        "wverb_t"    
#[16] "wverb_iq"    "wperf_t"     "wperf_iq"    "wfull4_t"    "wfull4"     
#[21] "wfull2_t"    "wfull2"      "Notes"      

## import my behavioral analyses of VGS/AS
demo$vgsDemo<-read.table(paste(paths$demo,"vgsDemo.txt",sep="/"),sep="\t",header=TRUE,fill=TRUE)
#names(demo$vgsDemo)
#[1] "ID"            "session"       "age"           "sex"          
#[5] "ID_ind"        "Good.1"        "Notes_exclude"
demo$vgsData<-read.table(paste(paths$demo,"vgsData.txt",sep="/"),sep="\t",header=TRUE,fill=TRUE)
#names(demo$vgsData)
# [1] "numCorr"  "mRT1"     "sdRT1"    "CV1"      "mu1"      "sigma1"  
# [7] "tau1"     "freq05.1" "mRT2"     "sdRT2"    "CV2"      "mu2"     
#[13] "sigma2"   "tau2"     "freq05.2" "freq05.3" "freq05.4"

demo$antiDemo<-read.table(paste(paths$demo,"antiDemo.txt",sep="/"),sep="\t",header=TRUE,fill=TRUE)
#names(demo$antiDemo)
#[1] "ID"            "session"       "age"           "sex"          
#[5] "ID_ind"        "Good.1"        "Notes_exclude"
demo$antiData<-read.table(paste(paths$demo,"antiData.txt",sep="/"),sep="\t",header=TRUE,fill=TRUE)
#names(demo$antiData)
# [1] "numCorr"     "numErr"      "percErr"     "mRTcorr1"    "sdRTcorr1"  
# [6] "CVcorr1"     "mRTinc1"     "sdRTinc1"    "CVinc1"      "mu1"        
#[11] "sigma1"      "tau1"        "freq05corr1" "mRTcorr2"    "sdRTcorr2"  
#[16] "CVcorr2"     "mRTinc2"     "sdRTinc2"    "CVinc2"      "mu2"        
#[21] "sigma2"      "tau2"        "freq05corr2" "freq05corr3" "freq05all3" 
#[26] "freq05corr4" "freq05all4" 

## matching demographics data to scan data
## excluding ineligible subjects/scans
data<-list()

## from imaging data (note - demo$unique contains all unique lunaID's for all subjects included in study)
data$scans<-dir(paths$data)
data$id<-substr(data$scans,1,5)
data$date<-substr(data$scans,7,14)

## which scans to include
include<-logical(length(data$scans))==FALSE ## initializing to TRUE, which will make it easier to use
## NOTE: THIS SECTION IS NOW UNNECESSARY, SCANS ALREADY REMOVED
## subjects to exclude
#exclude<-list()
## NOTE: due to diagnosis of self or family member having neurological or psychiatric disorder; otherwise "not part of study"
#exclude$subj<-as.character(c(10134,10164,10165,10166,10171,10232,10233,10267,10268,10269,10275,10278,10279,10282,10297,10306,10309,10323,10330,10331,10342,10346,10353,10354,10362,10369,10380,10435,10437))
#for(i in 1:length(exclude$subj)) include[data$id==exclude$subj[i]]<-FALSE
## scans to exclude
## NOTE: excluding a single scan which was not properly oriented (10220_20060211)
#exclude$scan<-list(id=as.character(c(10220)),date=as.character(c(20060211)))
#for(i in 1:length(exclude$scan$id)) include[intersect(which(data$id==exclude$scan$id[i]),which(data$date==exclude$scan$date[i]))]<-FALSE
## 10148 is a subject with strange timing between scans, leading to 6 legit scans in 5 years, chose the most spaced out 5 scans
#include[which(data$id=="10148")]<-c(TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,TRUE,TRUE)

## lines up demographics data with scan data to get age, exclude scans <3 months after previous valid scan within-subject, also exclude subjects with only 1 scan (if exclude.single.scans==TRUE)
data$age<-as.numeric(rep(NA,length(data$scans)))
for(i in 1:length(data$scans)){
	## only pull data if subject not already excluded
	if(include[i]==TRUE){
		## age
		temp.age<-difftime(strptime(data$date[i],"%Y%m%d"),strptime(demo$unique$DateOfBirth[which(demo$unique$LunaID==data$id[i])],"%m/%d/%Y"))
		data$age[i]<-round(unclass(temp.age)/365.25,2)
		## exclusion based on time between scans
		## NOTE: FOR NOW, THIS IS ALSO UNNECESSARY, TAKEN CARE OF PREVIOUSLY	
		## NOTE: 6 months (1/2 year) here
		## NOTE 2: looks back in scans until either:
			## 1) a different subject is encountered (include[i]==TRUE)
			## 2) a valid scan <3 months previously is encountered (include[i]==FALSE)
			## 3) the scans are checked all the way back to the first (include[i]==TRUE)
		## then breaks	 
		#if(length(which(data$id==data$id[i]))>=1){
		#	j<-1
		#	while(i>j){
		#		if(data$id[i]==data$id[i-j]){
		#			if((data$age[i]-data$age[i-j]<0.5)&(include[i-j]==TRUE)){
		#				include[i]<-FALSE
		#				data$age[i]<-NA
		#				break
		#			}
		#		}else{
		#			break
		#		}
		#		j<-j+1
		#	}
		#}
	}
}

## 10428 is missing a birthday, was 25 at time of scan, temporarily making age 25 until birthday is found
data$age[324]<-25

## organizing/printing info about data/demographics
id<-data$id[include]
age<-data$age[include]
uniq.id<-unique(id)
age.start<-numeric(length(uniq.id))
nscans<-numeric(length(uniq.id))
nscans.ind<-numeric(length(id))
nscans.tot<-numeric(length(id))
for(i in 1:length(uniq.id)){
	ind.id<-which(id==uniq.id[i])
	age.start[i]<-age[min(ind.id)]
	nscans[i]<-length(ind.id)
	nscans.ind[ind.id]<-1:length(ind.id)
	nscans.tot[ind.id]<-nscans[i]
}

## write to file
sink(paste(paths$study,"age.info",sep="/"))
cat(sprintf(gettext("Participants: %i\nAge (start of study): %.1f +/- %.1f\n# scans per subject:\n" ),length(uniq.id),mean(age.start),sd(age.start)))
for(i in 1:max(nscans)) cat(sprintf(gettext("\t%i scans: %i\n"),i,length(which(nscans==i))))
#cat(sprintf(gettext("# included scans: %i\nAge (all included scans): %.1f +/- %.1f\n# excluded scans: %i\n" ),length(id),mean(age),sd(age),length(which(!include))))
cat("\nIf excluding subjects with n scans, here will be the new sample age ranges at start of study:\n")
for(i in 2:max(nscans)){
	temp<-which(nscans>=i)
	cat(sprintf(gettext("%i... Participants: %i\nAge (start of study): %.1f +/- %.1f\n# scans per subject:\n" ),i,length(uniq.id[nscans>=i]),mean(age.start[nscans>=i]),sd(age.start[nscans>=i])))
}
sink()

## setting up other variables
## variables from unique ID demographics file (ie, sex)
ind.uniq<-as.numeric(rep(NA,length(id)))
## puberty
ind.pub<-as.numeric(rep(NA,length(id)))
## IQ
ind.wasi<-as.numeric(rep(NA,length(id)))
## behavioral data
ind.vgs<-as.numeric(rep(NA,length(id)))
ind.anti<-as.numeric(rep(NA,length(id)))

## lines up demographic/behavioral data with scan data
for(i in 1:length(id)){

	## uniq
	ind.uniq[i]<-which(demo$unique$LunaID==id[i])

	## puberty
	temp.ind.puberty<-which(demo$puberty$LunaID==id[i])
	if(length(temp.ind.puberty)>0){
		ind.pub[i]<-temp.ind.puberty[which.min(abs(demo$puberty$Age[temp.ind.puberty]-age[i]))]
		## if >3 months, exclude
		if(abs(demo$puberty$Age[ind.pub[i]]-age[i])>.25) ind.pub[i]<-NA
	}else{
		ind.pub[i]<-NA
	}

	## iq
	temp.ind.wasi<-which(demo$wasi$LunaID==id[i])
	if(length(temp.ind.wasi)>0){
		ind.wasi[i]<-temp.ind.wasi[which.min(abs(demo$wasi$Age[temp.ind.wasi]-age[i]))]
		## if >3 months, exclude
		if(abs(demo$wasi$Age[ind.wasi[i]]-age[i])>.25) ind.wasi[i]<-NA
	}else{
		ind.wasi[i]<-NA
	}

	## behavioral
	
	## vgs
	temp.ind.vgs<-which(demo$vgsDemo$ID==id[i])
	if(length(temp.ind.vgs)>0){
		temp.ind.vgs.min<-which.min(abs(demo$vgsDemo$age[temp.ind.vgs]-age[i]))
		if(length(temp.ind.vgs.min)>0){ ## this might be nonzero if there is an NA in age; should fix this for rest as well (NOTE: need to use na.rm=TRUE)
			ind.vgs[i]<-temp.ind.vgs[temp.ind.vgs.min]
			## if >3 months, exclude
			if(abs(demo$vgsDemo$age[ind.vgs[i]]-age[i])>.25) ind.vgs[i]<-NA
		}else{
			ind.vgs[i]<-NA
		}
	}else{
		ind.vgs[i]<-NA
	}
	
	## anti
	temp.ind.anti<-which(demo$antiDemo$ID==id[i])
	if(length(temp.ind.anti)>0){
		temp.ind.anti.min<-which.min(abs(demo$antiDemo$age[temp.ind.anti]-age[i]))
		if(length(temp.ind.anti.min)>0){ ## this might be nonzero if there is an NA in age
			ind.anti[i]<-temp.ind.anti[temp.ind.anti.min]
			## if >3 months, exclude
			if(abs(demo$antiDemo$age[ind.anti[i]]-age[i])>.25) ind.anti[i]<-NA
		}else{
			ind.anti[i]<-NA
		}
	}else{
		ind.anti[i]<-NA
	}
}

## NOTE TO SELF: "custom code" section where minor changes can be made?
## need option to make continuous variables into variable size groups with cutoffs

## sex (converting to m/f from 1/2)
sex<-ifelse(demo$unique$SexID[ind.uniq]==1,"m","f")
## tanner scale (puberty)
tsr<-ifelse(ind.pub>0,demo$puberty$tsr[ind.pub],NA)
tsr.2<-as.factor(ifelse(tsr<4,0,1))
temp<-ifelse(tsr<3,-1,0)
tsr.3<-as.factor(ifelse(tsr<5,temp,1))
## verbal and performance iq (wasi)
viq<-ifelse(ind.wasi>0,demo$wasi$wverb_iq[ind.wasi],NA)
temp<-ifelse(viq<100,-1,0)
viq.3<-as.factor(ifelse(viq<120,temp,1))
piq<-ifelse(ind.wasi>0,demo$wasi$wperf_iq[ind.wasi],NA)
temp<-ifelse(piq<100,-1,0)
piq.3<-as.factor(ifelse(piq<120,temp,1))
## behavioral data
## for vgs, exclude if number of correct trials <32 (2/3)
vgs.numCorr<-ifelse(ind.vgs>0,demo$vgsData$numCorr[ind.vgs],NA)
ind.vgs[which(vgs.numCorr<32)]<-NA
vgs.mRT<-ifelse(ind.vgs>0,demo$vgsData$mRT1[ind.vgs],NA)
vgs.sdRT<-ifelse(ind.vgs>0,demo$vgsData$sdRT1[ind.vgs],NA)
vgs.cv<-ifelse(ind.vgs>0,demo$vgsData$CV1[ind.vgs],NA)
vgs.mu<-ifelse(ind.vgs>0,demo$vgsData$mu1[ind.vgs],NA)
vgs.sigma<-ifelse(ind.vgs>0,demo$vgsData$sigma1[ind.vgs],NA)
vgs.tau<-ifelse(ind.vgs>0,demo$vgsData$tau1[ind.vgs],NA)
vgs.slow4<-ifelse(ind.vgs>0,demo$vgsData$freq05.4[ind.vgs],NA)
## for anti, exclude if number of correct+error trials <32 (2/3)
anti.numCorr<-ifelse(ind.anti>0,demo$antiData$numCorr[ind.anti],NA)
anti.numErr<-ifelse(ind.anti>0,demo$antiData$numErr[ind.anti],NA)
ind.anti[which(anti.numCorr + anti.numErr<32)]<-NA
anti.percErr<-ifelse(ind.anti>0,demo$antiData$percErr[ind.anti],NA)
## also excluding 2sd outliers for percent errors (~2/3, 67.6%)
ind.anti[which(anti.percErr>mean(anti.percErr, na.rm=TRUE)+2*sd(anti.percErr, na.rm=TRUE))]<-NA
anti.numErr<-ifelse(ind.anti>0,demo$antiData$numErr[ind.anti],NA)
anti.percErr<-ifelse(ind.anti>0,demo$antiData$percErr[ind.anti],NA)
anti.mRT<-ifelse(ind.anti>0,demo$antiData$mRTcorr1[ind.anti],NA)
anti.sdRT<-ifelse(ind.anti>0,demo$antiData$sdRTcorr1[ind.anti],NA)
anti.cv<-ifelse(ind.anti>0,demo$antiData$CVcorr1[ind.anti],NA)
anti.mu<-ifelse(ind.anti>0,demo$antiData$mu1[ind.anti],NA)
anti.sigma<-ifelse(ind.anti>0,demo$antiData$sigma1[ind.anti],NA)
anti.tau<-ifelse(ind.anti>0,demo$antiData$tau1[ind.anti],NA)
anti.slow4corr<-ifelse(ind.anti>0,demo$antiData$freq05corr4[ind.anti],NA)
anti.slow4all<-ifelse(ind.anti>0,demo$antiData$freq05all4[ind.anti],NA)

## organizing/printing info about data/demographics
var1<-data.frame(tsr.3, viq.3, piq.3, vgs.numCorr, vgs.mRT, vgs.sdRT, vgs.cv, vgs.mu, vgs.sigma, vgs.tau, vgs.slow4, anti.numCorr, anti.numErr, anti.percErr, anti.mRT, anti.sdRT, anti.cv, anti.mu, anti.sigma, anti.tau, anti.slow4corr, anti.slow4all)
var2<-data.frame(sex)
text.file<-paste(paths$study,"var.info",sep="/")

## prints info about variables
sample.info<-function(var1,var2,text.file,long=age,id=id){
	## only info where there is longitudinal info
	ind.long<-which(!is.na(long))
	id<-id[ind.long]
	long<-long[ind.long]
	var1<-data.frame(var1[ind.long,])
	var2<-data.frame(var2[ind.long,])
	## write to text file
	sink(text.file)
	if(length(var1)>= 1){		
		for(i in 1:dim(var1)[2]){
			ind<-which(!is.na(var1[,i]))
			cat(sprintf(gettext("%s: n = %i, n.scans = %i (%s = %.1f +/- %.1f %s)\n"),names(var1)[i],length(unique(id[ind])),length(ind),"age",mean(long[ind]),sd(long[ind]),"years"))
			if(length(var2)>=1){
				for(j in 1:dim(var2)[2]){
					var.lev<-levels(var2[ind,j])
					if(length(var.lev)>=1){
						for(k in 1:length(var.lev)){
							ind.lev<-ind[which(var2[ind,j]==var.lev[k])]
							cat(sprintf(gettext("\t%i %s (%s = %.1f +/- %.1f %s)\n"),length(var2[ind.lev, j]),var.lev[k],"age",mean(long[ind.lev]),sd(long[ind.lev]),"years"))
						}
					} #else {## error message here}
				}
			} #else {## error message here}
			cat("\n")
		}
	} #else {## error message here}
	sink()
}
