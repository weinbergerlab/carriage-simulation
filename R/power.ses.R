power.ses<-function(ds){
  ses.effect<- as.data.frame(t(sapply(ds, function(x) as.matrix(x[x$covar=='ses.grp.vector',1:3]))))
  names(ses.effect)<-names(ds[[1]])[1:3]
  sig.effect<- ses.effect$or.lcl>1
  return(sig.effect)
}