##Use markov transition model
markov.func<-function(ds, ses.key, hh.txn=F){
  Q <- rbind ( c(-1, 1),
               + c(1, -1))
  ds<-merge(ds,ses.key, by='id')
  n.people<-length(unique(ds$id))
  ds$date<-as.Date('2019/9/30') + ds$sample.index
  ds$qtr<-quarter(ds$date)
  ds$qtr2<-0
  ds$qtr3<-0
  ds$qtr4<-0
  ds$qtr2[ds$qtr==2]<-1
  ds$qtr3[ds$qtr==3]<-1
  ds$qtr4[ds$qtr==4]<-1
  ds.split<-split(ds, ds$id)
  n.obs<-sapply(ds.split, nrow)
  ds.split<-ds.split[n.obs>1]
  ds2<-lapply(ds.split, day.index.func )
  ds2<-do.call(rbind.data.frame, ds2)
  ds2$state<-ds2$col+1
  if(hh.txn==F){
    #test.msm <- msm( state ~ day.index, subject=id, covariates=list("1-2"= ~ ses.grp.vector, "1-2"=~qtr2,"1-2"=~qtr4), data = ds2, qmatrix=Q, gen.inits=TRUE)  
    test.msm <- msm( state ~ day.index, subject=id, covariates=list("1-2"= ~ ses.grp.vector), data = ds2, qmatrix=Q, gen.inits=TRUE)  
    acq.rate.ratio<-hazard.msm(test.msm)[['ses.grp.vector']]['State 1 - State 2',]       
    acq.rate.low<-
      if(length(acq.rate.ratio)<3){
        acq.rate.ratio<-c(acq.rate.ratio,NA,NA)
      }
    max.t<-max(ds2$day.index)
    prev.state.low.ses<-totlos.msm(test.msm, tot=max.t, covariates=list(ses.grp.vector=0), ci='normal')/max.t
    prev.state.hi.ses<-totlos.msm(test.msm, tot=max.t, covariates=list(ses.grp.vector=1), ci='normal') /max.t
    prob.transition.low.ses<-pmatrix.msm(test.msm, t=1, covariates=list(ses.grp.vector=0), ci='normal')
    prob.transition.high.ses<-pmatrix.msm(test.msm, t=1, covariates=list(ses.grp.vector=1), ci='normal')
    acq.rate.ratio.hh<- hazard.msm(test.msm)[['pair.col']]['State 1 - State 2',]
    
  }else if(hh.txn==T){
    
    test1<-ds
    time.points<- sum(ds$id==1) # how many samples per person
    test1$t<-rep(1:time.points, times=n.people)
    test2<-test1[,c('t','id','col')]
    test1.m<-melt(test2, id.vars=c('id','t'))
    test1.c<-acast(test1.m,id~t)
    pair<-rep(c(1:(n.people/2)), each=2)
    test1.c.spl<-split(test1.c,pair)
    test1.c.spl<-lapply(test1.c.spl, function(x) matrix(x, nrow=2))
    
    pair.col.spl<-lapply(test1.c.spl, prev.lag.func) #colonization state of the paired member at previous time
    pair.col<-do.call(rbind,pair.col.spl)
    ds3<-cbind.data.frame(ds2,'pair.col'=pair.col[,1] )
    # ds4<-ds3[ds3$pair.col==1 & !is.na(ds3$pair.col),]
    #ds4<-ds3[!is.na(ds3$pair.col),] #drop first observation for each person
    ds4<-ds3
    ds4$state[ds3$pair.col==1 & ds3$col==1] <-3
    #ds5<-ds4[ds4$state %in% c(2,3),]
    #ds5$state<-ds5$state-1
    ds5<-ds4
    ds5$state<-ds5$col+1
    # Q <- rbind ( c(-1, 1, 0),
    #               c(1, -1,1),
    #               c(0,1, -1))
    Q <- rbind ( c(-1, 1),
                 + c(1, -1))
    # msm( state ~ day.index, subject=id, covariates=list("1-2"= ~ ses.grp.vector, "1-2"=~qtr2,"1-2"=~qtr4, "1-2"=~pair.col), data = ds4, qmatrix=Q, gen.inits=TRUE)  
    #test.msm <- msm( state ~ day.index, subject=id, covariates=list("1-2"= ~ ses.grp.vector, "1-2"=~pair.col), data = ds3, qmatrix=Q, gen.inits=TRUE)  
    test.msm <- msm( state ~ day.index, subject=id, covariates=list("1-2"=~pair.col), data = ds5, qmatrix=Q, gen.inits=TRUE)  
    acq.rate.ratio.hh<- hazard.msm(test.msm)[['pair.col']]['State 1 - State 2',]
    acq.rate.ratio<-NA
    prob.transition.low.ses<-NA
    prob.transition.high.ses<-NA
    prev.state.low.ses<-NA
    prev.state.hi.ses<-NA
    
  }
  
  #sojourn.msm(test.msm, covariates=list(inc.grp.vec=1))
  #est.transition.comp<-qmatrix.msm(test.msm)
  
  
  out.list<-list(acq.rate.ratio.hh=acq.rate.ratio.hh  ,acq.rate.ratio=acq.rate.ratio,prob.transition.low.ses=prob.transition.low.ses,prob.transition.high.ses=prob.transition.high.ses ,
                 prev.state.low.ses=prev.state.low.ses,prev.state.hi.ses=prev.state.hi.ses )
  return(out.list)
}