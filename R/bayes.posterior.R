bayes.posterior <-
function(x,...){
  media<-apply(x,1,mean)
  mediana<-apply(x,1,median)
  intervalos<-coda::HPDinterval(as.mcmc(t(x)),...)
  resumen<-cbind(BE_mean=media,BE_median=mediana,intervalos)
  return(resumen)
}
