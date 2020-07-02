fitFreqdist=function (cell, distr){
  count=table(cell$Period)
  if(distr=="pois"){
    param=mean(as.numeric(count))
  }else if(distr=="nbinom"){
    nb_mle=as.numeric(MASS::fitdistr(as.numeric(count),"negative binomial")$estimate)
    param=c(nb_mle[1],nb_mle[1]/sum(nb_mle))
  }else{
    stop("The distribution has either to be 'pois' or 'nbinom'!")
  }
  return(buildFreqdist(distr, param))
}
