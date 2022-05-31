#Alternative specification of the markov model
markov.func.hh.cluster<-function(ds){
  
  #When setting up Q matrix, should calculate duration of time spent in each state (D_s), then inverse of this divided by number of transitions from the state gives the off-diagonal values
  #Also can set up the diagonal values to be 0 --they are defined in model as -(sum other values in row)
  #See PAGE 18 in the msm() manual
  # time in State 1 is 1-acq.rate
  
  state2.transition<-(1/duration+  (1-set.acq.rate))/2
  Q <- rbind ( c(0, 1, 0),   #go from uncolonized to singly colonized
               c(1, 0, 1),  #singly colonized to uncolonized or doubly colonized
               c(0, 1, 0)  #Double colonized to singly colonized
  ) 
  ds$date<-as.Date('2019/9/30') + ds$sample.index
  n.times<-nrow(ds[ds$id==1,])
  ds$hh.id<-rep(1:(n.people/2), each=2*n.times)
  
  ds.agg<-aggregate(ds$col, by=list('date'=ds$date, 'hh.id'=ds$hh.id), FUN=sum)
  names(ds.agg)[3]<-'states'
  ds.agg$qtr<-quarter(ds.agg$date)
  ds.agg$qtr2<-0
  ds.agg$qtr3<-0
  ds.agg$qtr4<-0
  ds.agg$qtr2[ds.agg$qtr==2]<-1
  ds.agg$qtr3[ds.agg$qtr==3]<-1
  ds.agg$qtr4[ds.agg$qtr==4]<-1
  ds.split<-split(ds.agg, ds.agg$hh.id)
  n.obs<-sapply(ds.split, nrow)
  ds.split<-ds.split[n.obs>1]
  ds2<-lapply(ds.split, day.index.func.date )
  ds2<-do.call(rbind.data.frame, ds2)
  ds2$day.index<-as.numeric(ds2$day.index)
  ds2$states<-ds2$states+1
  ds2$ses<-ses.hh.vector
  test.msm <- msm( states ~ day.index, subject=hh.id, covariates=list("1-2"= ~ ses, '2-3'=~ ses), data = ds2, 
                   qmatrix=Q, gen.inits=TRUE, opt.method='bobyqa') #control=list(fnscale=1000)  
  
  #sojourn.msm(test.msm, covariates=list(inc.grp.vec=1))
  #est.transition.comp<-qmatrix.msm(test.msm)
  acq.rate.ratio<-hazard.msm(test.msm)[['ses']]
  
  acq.rate.low<- 
    if(length(acq.rate.ratio)<3){
      acq.rate.ratio<-c(acq.rate.ratio,NA,NA)
    }
  max.t<-max(ds2$day.index)
  prev.state.low.ses<-totlos.msm(test.msm, tot=max.t, covariates=list(ses=0), ci='normal')/max.t 
  prev.state.hi.ses<-totlos.msm(test.msm, tot=max.t, covariates=list(ses=1), ci='normal') /max.t 
  prob.transition.low.ses<-pmatrix.msm(test.msm, t=1, covariates=list(ses=0), ci='normal')
  prob.transition.high.ses<-pmatrix.msm(test.msm, t=1, covariates=list(ses=1), ci='normal')
  out.list<-list(q.mat=test.msm$Qmatrices, acq.rate.ratio=acq.rate.ratio,prob.transition.low.ses=prob.transition.low.ses,prob.transition.high.ses=prob.transition.high.ses ,prev.state.low.ses=prev.state.low.ses,prev.state.hi.ses=prev.state.hi.ses )
  return(out.list)
}