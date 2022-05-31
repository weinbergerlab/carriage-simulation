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