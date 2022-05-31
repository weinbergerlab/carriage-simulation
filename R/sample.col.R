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