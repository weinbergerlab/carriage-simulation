
plot_gee<-function(ds){
  ses.effect<- as.data.frame(t(sapply(ds, function(x) as.matrix(x[x$covar=='ses.grp.vector',1:3]))))
  names(ses.effect)<-names(ds[[1]])[1:3]
  plot(1:nrow(ses.effect),ses.effect$or, ylim=c(0.2,2), bty='l')
  abline(h=1, lty=2, col='gray')
  abline(h=1/1.5, lty=3, col='red')
  arrows(1:nrow(ses.effect), ses.effect$or.lcl, 1:nrow(ses.effect), ses.effect$or.ucl, length=0.0, angle=90, code=3)
}