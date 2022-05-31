#Calculate the colonization state in the pair at t-1
prev.lag.func<-function(ds){
  prev<-ds
  prev.pair<-prev[c(2,1),] #flip the order
  prev.pair.lag<- cbind( matrix(NA, nrow=2, ncol=1),prev.pair)[,-(ncol(prev)+1)] #lag it by 1--does prevalence at previous visit predict current acquisition?
  prev.pair.lag.m<-melt(prev.pair.lag)
  prev.pair.lag.c<- acast(prev.pair.lag.m,Var1+Var2~.)
}
