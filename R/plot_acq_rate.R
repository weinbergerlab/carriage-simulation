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