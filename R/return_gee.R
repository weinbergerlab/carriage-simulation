return_gee<-function(samps){
  mod<-lapply(samps, function(x) gee_func(ds=x, ses.key=ses.grp.key) ) 
  plot_gee(mod)
  #Power
  title(paste("Power ",100*mean(power.ses(mod)), "%"  ))
}