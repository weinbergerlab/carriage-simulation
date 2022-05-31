#This function generates the markov chain for whole population
generate_replicates<-function( acq.rate,clear.rate , ses.rr,init.prev){
  unobs.truth<- Matrix(sapply(ses.grp,function(x) generate_data_func(ses.grp=x, acq.rate= acq.rate,clear.rate= clear.rate ,ses.rr=ses.rr, init.prev=init.prev),simplify=TRUE),sparse=TRUE)#Generate nsim versions of the dataset
  return(unobs.truth)
}