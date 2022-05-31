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