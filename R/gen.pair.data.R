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