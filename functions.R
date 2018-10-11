#V8: add in Markov

####
#######
#YCCI LONG: 50*8=400
#YCCI CROSS: 300
#POA CROSS: 200
#SENIOR CENTERS: 50people*3 centers*8 samples=1200---2 low income, 1 high income
#Nursing home: 50*8=400
#Senior housing: 100*8=800

#This function generates a markov chain for 1 person, with acquisition rates varying by SES and by quarter
generate_data_func<-function(ses.grp, acq.rate,clear.rate , ses.rr, init.prev){
  col.vec<- rep(NA, times=ntimes)
  if(ses.grp==1){ses.effect=ses.rr
  }else{ses.effect=1}
  col.vec[1]<- rbinom(n=1, size=1, prob=init.prev*ses.effect )
  for(i in 2:ntimes){
    acquire.change<- rbinom(n=1, size=1, prob=acq.rate*ses.effect*qtr.rr.vec[i])
    clear.change <-   rbinom(n=1, size=1, prob=clear.rate)
    col.vec[i]<- col.vec[i-1] # Same state as previously
    col.vec[i][col.vec[i-1]==0 & acquire.change==1] <- 1 #If uncolonized previously, colonize?
    col.vec[i][col.vec[i-1]==1 & clear.change==1] <- 0 # if colonized previously, clear?
  }  
  return(col.vec)
}

#This function generates the markov chain for whole population
generate_replicates<-function( acq.rate,clear.rate , ses.rr,init.prev){
  unobs.truth<- Matrix(sapply(ses.grp,function(x) generate_data_func(ses.grp=x, acq.rate= acq.rate,clear.rate= clear.rate ,ses.rr=ses.rr, init.prev=init.prev),simplify=TRUE),sparse=TRUE)#Generate nsim versions of the dataset
  return(unobs.truth)
}

#Take samples
sample.col<-function(ds,sample.intervals,proportion.longitudinal,n.cross.people,n.long.people ){
  possible.start.dates<-c(seq(as.Date('2019/10/01'), as.Date('2020/04/01'), by="day"), seq(as.Date('2020/10/01'), as.Date('2021/04/01'), by="day"))
  all.dates<- seq(as.Date('2019/10/01'), length.out = ntimes , by="day")
  n.people.sampled=n.cross.people+n.long.people
  first.sample.dates<-sample(possible.start.dates, size=n.people.sampled, replace=TRUE)
  time.interval <- '2019/9/30' %--% first.sample.dates
  first.sample.index <- as.period(time.interval, unit='day')@day #index at which first sample is taken for each person
  repeated.sample<-rep(0,times=n.people) #assignment of person to be repeatedly or cross sectionally smapled
  repeated.sample[1:n.long.people]<-1 
  
  obs.col<-  vector("list", n.people.sampled) 
  for(i in 1:n.people.sampled){
    if(repeated.sample[i]==1){
      sample.indices<-  first.sample.index[i]+sample.intervals
    }else{
      sample.indices<-first.sample.index[i]
    }
    obs.col.select<- as.matrix(ds[sample.indices,i, drop=FALSE])
    obs.col[[i]]<-cbind.data.frame(id=rep(i, times=length(sample.indices)) ,sample.index=sample.indices,col=obs.col.select)
    obs.col[[i]]$repeated.samples<-repeated.sample[i]
  } 
  obs.col.all<-  do.call("rbind", obs.col)
  return(obs.col.all)
}

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

plot_gee<-function(ds){
  ses.effect<- as.data.frame(t(sapply(ds, function(x) as.matrix(x[x$covar=='ses.grp.vector',1:3]))))
  names(ses.effect)<-names(ds[[1]])[1:3]
  plot(1:nrow(ses.effect),ses.effect$or, ylim=c(0.2,2), bty='l')
  abline(h=1, lty=2, col='gray')
  abline(h=1/1.5, lty=3, col='red')
  arrows(1:nrow(ses.effect), ses.effect$or.lcl, 1:nrow(ses.effect), ses.effect$or.ucl, length=0.0, angle=90, code=3)
}
power.ses<-function(ds){
  ses.effect<- as.data.frame(t(sapply(ds, function(x) as.matrix(x[x$covar=='ses.grp.vector',1:3]))))
  names(ses.effect)<-names(ds[[1]])[1:3]
  sig.effect<- ses.effect$or.ucl<1
  return(sig.effect)
}

return_gee<-function(samps){
  mod<-lapply(samps, function(x) gee_func(ds=x, ses.key=ses.grp.key) ) 
  plot_gee(mod)
  #Power
  title(paste("Power ",100*mean(power.ses(mod)), "%"  ))
}


    ################################################################################
    ################################################################################
    ###############################################################################
    ###########NOW TAKE SOME SAMPLES FROM TE FULL POPULATION TO EVALUATE DIFFERNT SAMPLING SCHEMES
    #Which people are sampled on which dates?

    ##Use markov transition model
    markov.func<-function(ds, ses.key){
      Q <- rbind ( c(-1, 1),
                   + c(1, -1))
      ds<-merge(ds,ses.key, by='id')
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
      test.msm <- msm( state ~ day.index, subject=id, covariates=list("1-2"= ~ ses.grp.vector, "1-2"=~qtr2,"1-2"=~qtr4), data = ds2, qmatrix=Q, gen.inits=TRUE)  
      #sojourn.msm(test.msm, covariates=list(inc.grp.vec=1))
      #est.transition.comp<-qmatrix.msm(test.msm)
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
      out.list<-list(acq.rate.ratio=acq.rate.ratio,prob.transition.low.ses=prob.transition.low.ses,prob.transition.high.ses=prob.transition.high.ses ,prev.state.low.ses=prev.state.low.ses,prev.state.hi.ses=prev.state.hi.ses )
      return(out.list)
    }
   
    day.index.func<-function(ds1){
      ds1$day.index<- ds1$sample.index - min(ds1$sample.index)
      return(ds1)
    }
    #Extract hazard ratios
    return_acq_rate_ratio<-function(ds){
      ds$acq.rate.ratio
    }
    #Plot hazard ratios
    plot_acq_rate_ratio<-function(ds){
      ds2<- as.data.frame(t( sapply(ds,return_acq_rate_ratio) ))
      power<- mean(ds2$U<1)*100
      plot(1:nrow(ds2),ds2$HR, ylim=c(0.2,2), bty='l', ylab="Hazard Ratio", xlab="Simulation")
      abline(h=1, lty=2, col='gray')
      abline(h=1/1.5, lty=3, col='red')
      arrows(1:nrow(ds2), ds2$L, 1:nrow(ds2), ds2$U, length=0.0, angle=90, code=3)
      title(paste0("Power: ", power,"%"))
    }
    #Extract transition.prob
    return_transition_prob<-function(ds){
     high.ses<-cbind.data.frame(median=ds$prob.transition.high.ses$estimates[1,2],L=ds$prob.transition.high.ses$L[1,2],U=ds$prob.transition.high.ses$U[1,2], ses="high")
     low.ses<-cbind.data.frame(median=ds$prob.transition.low.ses$estimates[1,2],L=ds$prob.transition.low.ses$L[1,2],U=ds$prob.transition.low.ses$U[1,2], ses="Low")
     duration<-cbind.data.frame(median=1/ds$prob.transition.low.ses$estimates[2,1],L=1/ds$prob.transition.low.ses$L[2,1],U=1/ds$prob.transition.low.ses$U[2,1], ses="Duration")
     acquisition.rates<-rbind.data.frame(high.ses,low.ses,duration)
     return(acquisition.rates)
     }
    #Plot durations
    plot_duration<-function(ds,duration){
      ds2<-  lapply(ds,function(x) return_transition_prob(x) ) 
      ds2<- lapply(ds2, function(x) x[x$ses=='Duration',1:3])   
      ds2<- do.call(rbind, ds2)
      
      power<- mean(  ds2$L>=duration & ds2$U<=duration)*100
      plot(1:nrow(ds2),ds2$median, ylim=c(5,40), bty='l', ylab="Duration", xlab="Simulation")
      abline(h=duration, lty=3, col='red')
      arrows(1:nrow(ds2), ds2$L, 1:nrow(ds2), ds2$U, length=0.0, angle=90, code=3)
      title(paste0("Power: ", power,"%"))
    }
    #Plot acquistion rates
    plot_acq_rate<-function(ds,set.acq.rate,set.ses.effect){
      ds1a<-  lapply(ds,function(x) return_transition_prob(x) ) 
      ds2.lo<- lapply(ds1a, function(x) x[x$ses=='Low',1:3])   
      ds2.lo<- do.call(rbind, ds2.lo)
      ds2.hi<- lapply(ds1a, function(x) x[x$ses=='high',1:3])   
      ds2.hi<- do.call(rbind, ds2.hi)
      
      power.lo<- mean(  ds2.lo$U>=set.acq.rate & ds2.lo$L<=set.acq.rate)*100
      power.hi<- mean(  ds2.hi$U>=set.acq.rate*1/set.ses.effect & ds2.hi$L<=set.acq.rate*1/set.ses.effect)*100
      par(mfrow=c(1,2))
     #plot1
      plot.range<-range(cbind(ds2.hi,ds2.lo))
      plot(1:nrow(ds2.lo),ds2.lo$median, bty='l',ylim=plot.range, ylab="Duration", xlab="Simulation")
      abline(h=set.acq.rate, lty=3, col='red')
      arrows(1:nrow(ds2.lo), ds2.lo$L, 1:nrow(ds2.lo), ds2.lo$U, length=0.0, angle=90, code=3)
      title(paste0("Power Lo/High: ", power.lo,"%/ ", power.hi,"%"))
      #Plot2
      plot(1:nrow(ds2.hi),ds2.hi$median,ylim=plot.range, bty='l', ylab="Duration", xlab="Simulation")
      points(1:nrow(ds2.hi),ds2.hi$median, ylim=c(5,40), bty='l', ylab="Duration", xlab="Simulation")
      abline(h=set.acq.rate*1/set.ses.effect, lty=3, col='pink')
      arrows(1:nrow(ds2.hi), ds2.hi$L,1:nrow(ds2.hi), ds2.hi$U, length=0.0, angle=90, code=3)
      
    }
    
