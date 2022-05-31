
gee_func<-function(ds,ses.key ){
  ds<-merge(ds,ses.key, by='id')
  #reg.data$col<-as.factor(reg.data$col)
  ds$date<-as.Date('2019/9/30') + ds$sample.index
  ds$quarter<-as.factor(quarter(ds$date))
  mod1<-geeglm(  col~ quarter+ ses.grp.vector  , data=ds, family='binomial', id=id,corstr='independence')
  est <- esticon(mod1, diag(length(mod1$coefficients)))
  covar<-names(mod1$coefficients)
  OR.CI <- as.data.frame(round(exp(cbind(est$Estimate, est$Lower, est$Upper)),2))
  names(OR.CI)<-c('or', 'or.lcl', 'or.ucl')
  OR.CI$covar<-covar
  return(OR.CI)
}