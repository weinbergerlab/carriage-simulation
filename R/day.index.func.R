day.index.func<-function(ds1){
  ds1$day.index<- ds1$sample.index - min(ds1$sample.index)+1
  return(ds1)
}