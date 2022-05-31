#Extract transition.prob
return_transition_prob<-function(ds){
  high.ses<-cbind.data.frame(median=ds$prob.transition.high.ses$estimates[1,2],L=ds$prob.transition.high.ses$L[1,2],U=ds$prob.transition.high.ses$U[1,2], ses="high")
  low.ses<-cbind.data.frame(median=ds$prob.transition.low.ses$estimates[1,2],L=ds$prob.transition.low.ses$L[1,2],U=ds$prob.transition.low.ses$U[1,2], ses="Low")
  duration<-cbind.data.frame(median=1/ds$prob.transition.low.ses$estimates[2,1],L=1/ds$prob.transition.low.ses$L[2,1],U=1/ds$prob.transition.low.ses$U[2,1], ses="Duration")
  acquisition.rates<-rbind.data.frame(high.ses,low.ses,duration)
  return(acquisition.rates)
}