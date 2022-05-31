day.index.func.date<-function(ds1){
  ds1$day.index<- ds1$date - min(ds1$date)+1
  return(ds1)
}