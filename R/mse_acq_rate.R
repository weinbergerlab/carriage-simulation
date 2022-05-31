mse_acq_rate<-function(ds,set.acq.rate,set.ses.effect){
  ds1a<-  lapply(ds,function(x) return_transition_prob(x) ) 
  ds2.lo<- lapply(ds1a, function(x) x[x$ses=='Low',1:3])   
  ds2.lo<- do.call(rbind, ds2.lo)
  ds2.hi<- lapply(ds1a, function(x) x[x$ses=='high',1:3])   
  ds2.hi<- do.call(rbind, ds2.hi)
  
  mse.lo<- mean(  (ds2.lo$median - set.acq.rate)^2)
  mse.hi<- mean(  (ds2.hi$median - set.acq.rate*1/set.ses.effect)^2)
  
  mse.combo<-c(mse.lo,mse.hi)
  return(mse.combo)
}
