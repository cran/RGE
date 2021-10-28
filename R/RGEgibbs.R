RGEgibbs <-
function(data,gen_c,env_c,blk_c,y_c,prior.g=NULL,prior.vg=NULL,
              prior.b=NULL,prior.dfb=NULL,prior.sv=NULL,prior.dfsv=NULL,
              prior.se=NULL, prior.dfse=NULL, burnin=10, thin=5, niter=50,
              saveAt=10){
  
  #Columas donde estan las variables
  #gen_c<-5
  #env_c<-1
  #blk_c<-4
  #y_c<-9
  
  #Ordenando el archivo
  data0<-data[order(data[,gen_c],data[,env_c]),]
  
  #prior
  #prior.g<-NULL   #media
  #prior.vg<-NULL  #media
  #prior.b<-NULL   #rep dentro env
  #prior.dfb<-NULL #rep dentro env credibility of degree
  #prior.sv<-NULL  #stability variance 
  #prior.dfsv<-NULL #stability variance credibility of degree 
  #prior.se<-NULL  #error variance for each env
  #prior.dfse<-NULL #error variance credibility of degree
  
  #burnin<-100
  #thin<-5
  #niter<-100
  #saveAt<-10
  #dividiendo dataset
  data_i<-split(data0,data0[,gen_c])
  data_j<-split(data0,data0[,env_c])
  
  #Numero de ambientes y genotipos
  nenv<-nlevels(as.factor(data0[,env_c]))
  ngen<-nlevels(as.factor(data0[,gen_c]))
  namesgen<-levels(as.factor(data0[,gen_c]))
  
  #Prior calulos fuera de iteracion
  if (is.null(prior.g)) {prior.g<-rep(0,ngen)}
  if (is.null(prior.vg)) {prior.vg<-rep(1e8,ngen)}
  #if (is.null(prior.b)) {prior.b<-0.1}
  #if (is.null(prior.dfb)) {prior.dfb<-1}
  if (is.null(prior.b)) {prior.b<-rep(0.1,ngen)}
  if (is.null(prior.dfb)) {prior.dfb<-rep(1,ngen)}
  if (is.null(prior.sv)) {prior.sv<-rep(0.1,ngen)}
  if (is.null(prior.dfsv)) {prior.dfsv<-rep(1,ngen)}
  if (is.null(prior.se)) {prior.se<-rep(0.1,nenv)}
  if (is.null(prior.dfse)) {prior.dfse<-rep(1,nenv)}
  
  ###Calculos no iterativos varios
  nenvi<-lapply(data_i, function (X) ncol(model.Matrix(~as.factor(X[,env_c])-1)))
  nri<-lapply(data_i, function (X) ncol(model.Matrix(~as.factor(X[,blk_c])-1)))
  nj<-lapply(data_j, nrow) 
  
  dfsestar<-list()
  sedfse<-list()
  for (j in 1:nenv){
    dfsestar[[j]]<-prior.dfse[j]+nj[[j]]
    sedfse[[j]]<-prior.dfse[j]*prior.se[j]
  }
  
  dfbstar<-list()
  bdfb<-list()
  for (i in 1:ngen){
    dfbstar[[i]]<-prior.dfb[i]+nri[[i]]*nenvi[[i]]
    bdfb[[i]]<-prior.dfb[i]*prior.b[i]
  }
  
  dfsvstar<-list()
  svdfsv<-list()
  for (i in 1:ngen){
    dfsvstar[[i]]<-prior.dfsv[i]+nenvi[[i]]
    svdfsv[[i]]<-prior.dfsv[i]*prior.sv[i]
  }
  
  Wi<-lapply(data_i,function (X) {
    tmp<-model.Matrix(~0+paste(X[,blk_c],X[,env_c]))
    tmp1<-model.Matrix(~0+as.factor(X[,env_c]))
    tmp2<-rep(1,nrow(X))
    names(tmp2)<-names(X)
    cbind(tmp2,tmp,tmp1)
  })
  
  #Valores iniciales para GIBBS
  sigmaej<-rep(15,nenv)+runif(nenv)
  #assign("sigmaej",rep(15,nenv)+runif(nenv),pos=1)
  #sigmar<-10+runif(1)
  sigmari<-rep(10,ngen)+runif(ngen)
  #assign("sigmari",rep(10,ngen)+runif(ngen),pos=1)
  sigmasv<-rep(5,ngen)+runif(ngen)
  #assign("sigmasv",rep(5,ngen)+runif(ngen),pos=1)
  
  ### gibbs desde aqui
  gibbs<-function (iter){
    Ri<-lapply(data_i, function (X) {
      tmp<-mapply(rep,times= table(X[,env_c]), x=1/sigmaej, SIMPLIFY = FALSE)
      tmp<-lapply(tmp,diag)
      bdiag(tmp)
    }) 
    
    Di<-list()
    for(i in 1:ngen){
      #Di[[i]]<-Diagonal(x=c(1/prior.vg[i],rep(1/sigmar,nri[[i]]*nenvi[[i]]),rep(1/sigmasv[i],nenvi[[i]])))
      Di[[i]]<-Diagonal(x=c(1/prior.vg[i],rep(1/sigmari[i],nri[[i]]*nenvi[[i]]),rep(1/sigmasv[i],nenvi[[i]])))
    }
    
    Hi<-list()
    for (i in 1:ngen){
      Hi[[i]]<-crossprod(Wi[[i]],Ri[[i]])%*%Wi[[i]]+Di[[i]]
    }
    tethai<-list()
    for (i in 1:ngen){
      media<-crossprod(Wi[[i]],Ri[[i]])%*%data_i[[i]][,y_c]
      media[1]<-media[1]+prior.g[i]*(1/prior.vg[i])
      tethastar<-crossprod(chol(Hi[[i]]),rnorm(ncol(Hi[[i]])))+media
      tethai[[i]]<-solve(Hi[[1]],tethastar)
    }
    
    #Para el error
    ei<-list()
    for (i in 1:ngen){
      ei[[i]]<-as.vector(data_i[[i]][,y_c]-Wi[[i]]%*%tethai[[i]])
    }
    e<-unlist(ei)
    
    ej<-split(e,data0[,env_c])
    ssej<-lapply(ej,crossprod)
    sigmaej<-list()
    for (j in 1: nenv){
      sigmaej[[j]]<-(1/rchisq(1,dfsestar[[j]]))*(ssej[[j]]+sedfse[[j]])
    }
    sigmaej<<-unlist(sigmaej)
    
    #Para b/env by gen
    sigmari<-list()
    for (i in 1: ngen){
      bi<-tethai[[i]][grep("paste",rownames(tethai[[i]]))]
      ssbi<-crossprod(bi)
      sigmari[[i]]<-(1/rchisq(1,dfbstar[[i]]))*(ssbi+bdfb[[i]])
    }
    sigmari<<-unlist(sigmari)
    
    #Para stability variance
    sigmasv<-list()
    for (i in 1: ngen){
      svi<-tethai[[i]][grep("factor",rownames(tethai[[i]]))]
      sssvi<-crossprod(svi)
      sigmasv[[i]]<-(1/rchisq(1,dfsvstar[[i]]))*(sssvi+svdfsv[[i]])
    }
    sigmasv<<-unlist(sigmasv)
    names(sigmasv)<-paste("sv",namesgen,sep="_")
    
    mugen<-unlist(lapply(tethai,function (X) X[1]))
    names(mugen)<-paste("Predicted",namesgen,sep="_")
    
    salida<-c(as.numeric(mugen),as.numeric(sigmasv))
    names(salida)<-c(paste("Predicted",namesgen,sep="_"),names(sigmasv)<-paste("sv",namesgen,sep="_"))
    return(salida)
  }
  
  tmp<-sapply(1:burnin,gibbs)
  
  outt<-list()
  for (iter in 1:niter) {
    tmp<-sapply(1:(thin-1),gibbs)
    outt[[iter]]<-sapply(1,gibbs)
    if ((iter%%saveAt==0)){
    save(outt,file="outtS4.RData")
    print(iter)
    }
  }
  outt<-do.call(cbind,outt)
  class(outt) <- "RGE"
  return(outt)
}
