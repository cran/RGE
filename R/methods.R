print.RGE<-function(x,...)
{
  sv<-x[grep("sv_",rownames(x)),]
  mu<-x[grep("Predicted_",rownames(x)),]
  bysi<-mu-1.96*sqrt(sv)
  namesbysi<-rownames(bysi)
  namesbysi<-gsub("Predicted","BYSI",namesbysi)
  rownames(bysi)<-namesbysi
  cat("--------------------------------------------\n")
  cat("Summary of predicted value of each genotypes\n")
  cat("--------------------------------------------\n")
  mu<-bayes.posterior(mu,...)
  print(mu)
  cat("\n")
  cat("\n")
  cat("-----------------------------------------------------\n")
  cat("Summary of stability variance value of each genotypes\n")
  cat("-----------------------------------------------------\n")
  sv<-bayes.posterior(sv,...)
  print(sv)
  cat("\n")
  cat("\n")
  cat("-----------------------------------------------------------\n")
  cat("Summary of bayesian yield stability index of each genotypes\n")
  cat("-----------------------------------------------------------\n")
  bysi<-bayes.posterior(bysi,...)
  print(bysi)

}

plot.RGE<-function(x,labelg="Predicted value",
                          labelsv="Stability variance",
                          labelby="Bayesian yield stability index",margin=c(1, 0.8, 0, 0.8),...)
  {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))  
  
  resumen<-summary.RGE(x,...)
  
  Y<-resumen[[1]]
  
  par(mai=margin)
  tmp<-Y[order(Y[,"BE_mean"],Y[,"lower"]),]
  delta<-(max(Y[,"upper"])-min(Y[,"lower"]))/50
  nm1<-nm2<-nm<-gsub("Predicted_","",rownames(tmp))
  tmp1<-tmp2<-tmp
  for (i in seq(1,length(nm),by=2)){nm1[i]<-NA;tmp1[i,]<-NA}
  for (i in seq(2,length(nm),by=2)){nm2[i]<-NA;tmp2[i,]<-NA}
  limites<-TRUE
  if (limites==TRUE){limites<-c(min(tmp[,"lower"],na.rm=TRUE)-delta,max(tmp[,"upper"],na.rm=TRUE)+delta)}
  plot(y=1:nrow(tmp1),x=tmp1[,"BE_mean"] , axes = FALSE,xlim=limites,xlab=labelg,ylab="",pch=18)
  axis(2, 1:nrow(tmp),nm1,las=2,cex.axis=0.9,tick=FALSE)
  axis(1)
  arrows ( y0=1:nrow(tmp1) ,x0=tmp1[,"upper"] , y1=1:nrow(tmp1) , x1=tmp1[,"lower"] , angle = 90 , code = 3 , length =0.01)
  points(y=1:nrow(tmp2),x=tmp2[,"BE_mean"] , pch=18, col="darkgreen")
  axis(4, 1:nrow(tmp2),nm2,las=2,cex.axis=0.9,tick=FALSE,col.axis="darkgreen")
  arrows ( y0=1:nrow(tmp2) ,x0=tmp2[,"upper"] , y1=1:nrow(tmp2) ,x1= tmp2[,"lower"] , angle = 90 , code = 3 , length =0.01,col="darkgreen")
  
  Y<-resumen[[2]]
  
  par(mai=margin)
  tmp<-Y[order(Y[,"BE_mean"],Y[,"lower"]),]
  delta<-(max(Y[,"upper"])-min(Y[,"lower"]))/50
  nm1<-nm2<-nm<-gsub("sv_","",rownames(tmp))
  tmp1<-tmp2<-tmp
  for (i in seq(1,length(nm),by=2)){nm1[i]<-NA;tmp1[i,]<-NA}
  for (i in seq(2,length(nm),by=2)){nm2[i]<-NA;tmp2[i,]<-NA}
  limites<-TRUE
  if (limites==TRUE){limites<-c(min(tmp[,"lower"],na.rm=TRUE)-delta,max(tmp[,"upper"],na.rm=TRUE)+delta)}
  plot(y=1:nrow(tmp1),x=tmp1[,"BE_mean"] , axes = FALSE,xlim=limites,xlab=labelsv,ylab="",pch=18)
  axis(2, 1:nrow(tmp),nm1,las=2,cex.axis=0.9,tick=FALSE)
  axis(1)
  arrows ( y0=1:nrow(tmp1) ,x0=tmp1[,"upper"] , y1=1:nrow(tmp1) , x1=tmp1[,"lower"] , angle = 90 , code = 3 , length =0.01)
  points(y=1:nrow(tmp2),x=tmp2[,"BE_mean"] , pch=18, col="darkgreen")
  axis(4, 1:nrow(tmp2),nm2,las=2,cex.axis=0.9,tick=FALSE,col.axis="darkgreen")
  arrows ( y0=1:nrow(tmp2) ,x0=tmp2[,"upper"] , y1=1:nrow(tmp2) ,x1= tmp2[,"lower"] , angle = 90 , code = 3 , length =0.01,col="darkgreen")

  Y<-resumen[[3]]
  
  par(mai=margin)
  tmp<-Y[order(Y[,"BE_mean"],Y[,"lower"]),]
  delta<-(max(Y[,"upper"])-min(Y[,"lower"]))/50
  nm1<-nm2<-nm<-gsub("BYSI_","",rownames(tmp))
  tmp1<-tmp2<-tmp
  for (i in seq(1,length(nm),by=2)){nm1[i]<-NA;tmp1[i,]<-NA}
  for (i in seq(2,length(nm),by=2)){nm2[i]<-NA;tmp2[i,]<-NA}
  limites<-TRUE
  if (limites==TRUE){limites<-c(min(tmp[,"lower"],na.rm=TRUE)-delta,max(tmp[,"upper"],na.rm=TRUE)+delta)}
  plot(y=1:nrow(tmp1),x=tmp1[,"BE_mean"] , axes = FALSE,xlim=limites,xlab=labelby,ylab="",pch=18)
  axis(2, 1:nrow(tmp),nm1,las=2,cex.axis=0.9,tick=FALSE)
  axis(1)
  arrows ( y0=1:nrow(tmp1) ,x0=tmp1[,"upper"] , y1=1:nrow(tmp1) , x1=tmp1[,"lower"] , angle = 90 , code = 3 , length =0.01)
  points(y=1:nrow(tmp2),x=tmp2[,"BE_mean"] , pch=18, col="darkgreen")
  axis(4, 1:nrow(tmp2),nm2,las=2,cex.axis=0.9,tick=FALSE,col.axis="darkgreen")
  arrows ( y0=1:nrow(tmp2) ,x0=tmp2[,"upper"] , y1=1:nrow(tmp2) ,x1= tmp2[,"lower"] , angle = 90 , code = 3 , length =0.01,col="darkgreen")

}

summary.RGE<-function(object,...)
{
  sv<-object[grep("sv_",rownames(object)),]
  mu<-object[grep("Predicted_",rownames(object)),]
  bysi<-mu-1.96*sqrt(sv)
  namesbysi<-rownames(bysi)
  namesbysi<-gsub("Predicted","BYSI",namesbysi)
  rownames(bysi)<-namesbysi
  
  mu<-bayes.posterior(mu)
  
  sv<-bayes.posterior(sv)
  
  bysi<-bayes.posterior(bysi)
  
  list(mu=mu,sv=sv,bysi=bysi)
  
}
