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