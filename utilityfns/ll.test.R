ll.test<-function(m1=models$model[[m]],m0=models$model[[models$con[[m]][c]]],na=comb.na,data=data.,mixed=setup$mixed,anova.=anova.txt){
## calculate LL test (if there are NAs, rerun models and then calculate LL test)
## details:
	## m1 (required) - model of interest
	## m0 (required) - model being tested against (null model)
	## data (required) - data frame if NAs and models need to be rerun
	## mixed (required) - TRUE if mixed model, FALSE if regular regression model

	## re-estimate models if NAs present
	if(length(na)>0){
		if(mixed==TRUE){
			m1<-lmer.est(m1$fixef,data=data[-na,])
			m0<-lmer.est(m0$fixef,data=data[-na,])
		}else{
			m1<-lm.est(m1$fixef,data=data[-na,])
			m0<-lm.est(m0$fixef,data=data[-na,])
		}
	}

	## log-likelihood test
	ll.test<-list()
	ll.test$m1<-logLik(m1$fit)
	ll.test$m0<-logLik(m0$fit)
	ll.test$LLR<--2*(ll.test$m0-ll.test$m1)
	ll.test$df<-attr(ll.test$LLR,"df")
	ll.test$p<-pchisq(ll.test$LLR,ll.test$df)

	## copy p-value into anova.txt
	sink(anova.,append=TRUE)
	cat(ll.test$p,"")
	sink()
	ll.test
}
