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