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
sample.col<-function(ds,sample.intervals,proportion.longitudinal,n.cross.people,n.long.people , hh=T){
  possible.start.dates<-c(seq(as.Date('2019/10/01'), as.Date('2020/04/01'), by="day"), seq(as.Date('2020/10/01'), as.Date('2021/04/01'), by="day"))
  all.dates<- seq(as.Date('2019/10/01'), length.out = ntimes , by="day")
  n.people.sampled=n.cross.people+n.long.people
  #If performing HH sampling, ensure members of the HH have same sample dates
  if(hh==F){
    first.sample.dates<-sample(possible.start.dates, size=n.people.sampled, replace=TRUE)
  }else{
    first.sample.dates<-sample(possible.start.dates, size=n.people.sampled/2, replace=TRUE)
    first.sample.dates<-rep(first.sample.dates, each=2)
  }
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
  sig.effect<- ses.effect$or.lcl>1
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
    markov.func<-function(ds, ses.key, hh.txn=F){
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
      if(hh.txn==F){
    #test.msm <- msm( state ~ day.index, subject=id, covariates=list("1-2"= ~ ses.grp.vector, "1-2"=~qtr2,"1-2"=~qtr4), data = ds2, qmatrix=Q, gen.inits=TRUE)  
     test.msm <- msm( state ~ day.index, subject=id, covariates=list("1-2"= ~ ses.grp.vector), data = ds2, qmatrix=Q, gen.inits=TRUE)  
        
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
       
        pair.col.spl<-lapply(test1.c.spl, prev.lag.func) #colonization state of the paired member
        pair.col<-do.call(rbind,pair.col.spl)
        ds3<-cbind.data.frame(ds2,'pair.col'=pair.col[,1] )
        ds4<-ds3[!is.na(ds3$pair.col),] #drop first observation for each person
        #ds4[,c('id','col','pair.col')]
       # test.msm <- msm( state ~ day.index, subject=id, covariates=list("1-2"= ~ ses.grp.vector, "1-2"=~qtr2,"1-2"=~qtr4, "1-2"=~pair.col), data = ds3, qmatrix=Q, gen.inits=TRUE)  
         test.msm <- msm( state ~ day.index, subject=id, covariates=list("1-2"= ~ ses.grp.vector, "1-2"=~pair.col), data = ds3, qmatrix=Q, gen.inits=TRUE)  
        
      }
      
      #sojourn.msm(test.msm, covariates=list(inc.grp.vec=1))
      #est.transition.comp<-qmatrix.msm(test.msm)
      acq.rate.ratio<-hazard.msm(test.msm)[['ses.grp.vector']]['State 1 - State 2',]
      acq.rate.ratio.hh<- hazard.msm(test.msm)[['pair.col']]['State 1 - State 2',]
      acq.rate.low<- 
      if(length(acq.rate.ratio)<3){
        acq.rate.ratio<-c(acq.rate.ratio,NA,NA)
      }
      max.t<-max(ds2$day.index)
      prev.state.low.ses<-totlos.msm(test.msm, tot=max.t, covariates=list(ses.grp.vector=0), ci='normal')/max.t 
      prev.state.hi.ses<-totlos.msm(test.msm, tot=max.t, covariates=list(ses.grp.vector=1), ci='normal') /max.t 
      prob.transition.low.ses<-pmatrix.msm(test.msm, t=1, covariates=list(ses.grp.vector=0), ci='normal')
      prob.transition.high.ses<-pmatrix.msm(test.msm, t=1, covariates=list(ses.grp.vector=1), ci='normal')
      out.list<-list(acq.rate.ratio.hh=acq.rate.ratio.hh,acq.rate.ratio=acq.rate.ratio,prob.transition.low.ses=prob.transition.low.ses,prob.transition.high.ses=prob.transition.high.ses ,prev.state.low.ses=prev.state.low.ses,prev.state.hi.ses=prev.state.hi.ses )
      return(out.list)
    }

#Calculate the colonization state in the pair at t-1
prev.lag.func<-function(ds){
  prev<-ds
  prev.pair<-prev[c(2,1),] #flip the order
  prev.pair.lag<- cbind( matrix(NA, nrow=2, ncol=1),prev.pair)[,-(ncol(prev)+1)] #lag it by 1--does prevalence at previous visit predict current acquisition?
  prev.pair.lag.m<-melt(prev.pair.lag)
  prev.pair.lag.c<- acast(prev.pair.lag.m,Var1+Var2~.)
}

    day.index.func<-function(ds1){
      ds1$day.index<- ds1$sample.index - min(ds1$sample.index)+1
      return(ds1)
    }
    day.index.func.date<-function(ds1){
      ds1$day.index<- ds1$date - min(ds1$date)+1
      return(ds1)
    }
    #Extract hazard ratios
    return_acq_rate_ratio<-function(ds){
      ds$acq.rate.ratio
    }
    #Plot hazard ratios
    plot_acq_rate_ratio<-function(ds){
      ds2<- as.data.frame(t( sapply(ds,return_acq_rate_ratio) ))
      power<- mean(ds2$L>1)*100
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
      power.hi<- mean(  ds2.hi$U>=set.acq.rate*set.ses.effect & ds2.hi$L<=set.acq.rate*set.ses.effect)*100
      par(mfrow=c(2,2))
     #plot1
      plot.range<-range(cbind(ds2.hi,ds2.lo))
      plot(1:nrow(ds2.lo),ds2.lo$median, bty='l',ylim=plot.range, ylab="Acq. Rate", xlab="Simulation")
      abline(h=set.acq.rate, lty=3, col='red')
      arrows(1:nrow(ds2.lo), ds2.lo$L, 1:nrow(ds2.lo), ds2.lo$U, length=0.0, angle=90, code=3)
      title(paste0("Power Lo/High: ", power.lo,"%/ ", power.hi,"%"))
      #Plot2
      plot(1:nrow(ds2.hi),ds2.hi$median,ylim=plot.range, bty='l', ylab="Acq. Rate", xlab="Simulation")
      points(1:nrow(ds2.hi),ds2.hi$median, ylim=c(5,40), bty='l')
      abline(h=set.acq.rate*set.ses.effect, lty=3, col='pink')
      arrows(1:nrow(ds2.hi), ds2.hi$L,1:nrow(ds2.hi), ds2.hi$U, length=0.0, angle=90, code=3)
      
      #Plot3
      hist(ds2.lo$median)
      abline(v=set.acq.rate, lty=3, col='red', main="Observed vs estimates Low SES")
      
      #Plot4
      hist(ds2.lo$median)
      abline(v=set.acq.rate, lty=3, col='red', main="Observed vs estimates High SES")
      
      
    }
    
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
    
    #Generate data from a simple transition model
    gen.pair.data<-function(acq.rate= set.acq.rate,clear.rate= 1/duration,prob.transmit=prob.transmit,duration=duration){
      prob.transmit.day<-prob.transmit/duration
      state<-matrix(NA, ncol=ntimes, nrow=n.people)
      pair<-rep(c(1:(n.people/2)), each=2)
      state[,1]<-rbinom(n=n.people, size=1, prob=init.prev)
      state.spl<-split(state,  pair) #split by pairing
      state.spl<-lapply(state.spl, function(x) matrix(x, nrow=2))
      
      for(i in 1:length(state.spl)){ #loop through pairs
        for(t in 2:ncol(state.spl[[i]])){ #loop through the dates
          
          for(k in 1:2){ #each member of the pair
            if(state.spl[[i]][k,(t-1)]==0 & sum(state.spl[[i]][ ,(t-1)])==0 ){ #uncolonized and both members of pair uncolonized at t-1
              state.spl[[i]][k,t]<-rbinom(n=1,size=1, prob=(set.acq.rate*set.ses.effect^ses.hh.vector[i] ))
            }
            
            if(state.spl[[i]][k,(t-1)]==0 & sum(state.spl[[i]][ ,(t-1)])==1 ){ #uncolonized #and pair member was colonized at t-1
              state.spl[[i]][k,t]<-rbinom(n=1,size=1, prob=(set.acq.rate*set.ses.effect^ses.hh.vector[i]+prob.transmit.day ))
              #state.spl[[i]][k,t] <- 99
            }  
            
            if(state.spl[[i]][k,(t-1)]==1  ){ #colonized at t-1
              state.spl[[i]][k,t]<- 1 - rbinom(n=1,size=1, prob=clear.rate )
            }   
            
          }
        }
      }
      state2<-do.call(rbind, state.spl) #True states
      state2<-t(state2)
      return(state2)
    }
    
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
    
    
    format.jags.3state<-function(ds){
      
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
      return(ds2)
    }
    
    
  jags.func<-function(ds){
      n.times<-nrow(ds[ds$id==1,])
      ds$hh.id<-rep(1:(n.people/2), each=2*n.times)
      ds.split<-split(ds, ds$id)
      ds2<-lapply(ds.split, day.index.func )
      ds2<-do.call(rbind.data.frame, ds2)
      ds2$day.index<-as.numeric(ds2$day.index)
      ds2$index<- rep( rep(c(1,2), each=n.times) , times=max(ds$hh.id))   
      #3 dimensions: Household, individual, time
      test1.c<- acast( ds2[,c('index','hh.id', 'day.index','col')], hh.id~index~day.index,value.var='col' )
      days<-as.numeric(dimnames(test1.c)[[3]])
      time.increment<- days-lag(days)
      
      prev.pop<-mean(ds$col)
      
      model_string <-"
      model{
      
      for(i in 1:N.hh){
      for(j in 1:2){
      state[i,j,1] ~ dbern(p[i,j,1] )
      p[i,j,1]<- prev.pop

      for(k in 2:6){ 
      state[i,j,k] ~ dbern(p[i,j,k] )
      p[i,j,k]<- ( prob_1_2[k]*(1-state[i,j,k-1])*(1-state[i,switch[j],k-1]) #individual and partner uncolonized
      + prob_2_3[k]*(1-state[i,j,k-1])*(state[i,switch[j],k-1])  #individual uncolonized, partner colonized
      + (1-prob_2_1[k])*state[i,j,k-1]   # individual colonized, regardless of states
      )
      
      }
      }
      } 
      
      for(k in 2:6){
      prob_2_1[k]<- ilogit(beta_2_1)*time.increment[k]/14 #give probaility for 14 days
      prob_2_3[k]<- prob_1_2[k] + ilogit(beta_within_hh)*time.increment[k]/14
      prob_1_2[k]<- ilogit(beta_1_2)*time.increment[k]/14
      }
      #If truncate these priors at 0, it ensures that exp(beta_i_j) falls within the range [0,1]
      beta_1_2~dnorm(0,1e-2)T(,0) #log-acq rate
      beta_2_1~dnorm(0,1e-2)T(,0)
      #beta_2_3<- beta_1_2+ exp(beta_within_hh) #log_prob for 2_3 equal prob for 1_2 plus within HH effect
      beta_within_hh ~ dnorm(0 ,1e-2)T(,0)
      switch<-c(2,1)
      #prob_within_hh_day<-exp(beta_2_3-beta_1_2 )

    #exp(beta) gives prob change per 2 weeks
      prob_2_1_day<-exp(beta_2_1)/14
      prob_1_2_day<-exp(beta_1_2)/14
      #prob_2_3_day<-exp(beta_2_3)/14
      #logit_prob_within_hh_day<- logit(prob_2_3_day - prob_1_2_day)

      }
      "
      model_jags<-jags.model(textConnection(model_string),
                             data=list('N.hh' = max(ds2$hh.id),  
                                       'state'=test1.c,'prev.pop'=prev.pop,
                                       'time.increment'=time.increment
                             )) 
      
      update(model_jags, 
             n.iter=5000) 
      
      posterior_samples<-coda.samples(model_jags, 
                                      variable.names=c('prob_2_1_day','prob_1_2_day','beta_2_1','beta_1_2','beta_within_hh'),
                                      
                                      thin=1,
                                      n.iter=5000)
      
      post.q.prob.day<-t(apply(posterior_samples[[1]],2,quantile, probs=c(0.025,0.5,0.975)))
      #duration.in.state<-1/post.q.prob.day
      
      jags.results<-list('post.q.prob.day'=post.q.prob.day, 'posterior_samples'=posterior_samples)
   return(jags.results)
       }