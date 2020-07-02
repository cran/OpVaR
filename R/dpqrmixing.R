dmixing=function (x, sevdist, log = FALSE){
  cmu=sevdist$par[[3]][[5]][1]
  ctau=sevdist$par[[3]][[5]][2]
  
  rx <- function(x, cmu, ctau) {
    (sevdist$par[[1]][[1]](x) - sevdist$par[[2]][[1]](x)) * 
      atan((x - cmu)/ctau)
  }
  d = x
  
  r = try(pracma::quadinf(function(x){rx(x, cmu=cmu, ctau = ctau)}, xa = 0, xb=Inf)$Q)
  
  pweights = pcauchy(x, cmu, ctau)
  d = ((1 - pweights) * sevdist$par[[2]][[1]](x) + pweights * 
    sevdist$par[[1]][[1]](x))/(1+r/pi)
  if (log) 
    d = log(d)
  return(d)
}

pmixing=function (q, sevdist, lower.tail = TRUE){
  sapply(q, function(q){
  cmu=sevdist$par[[3]][[5]][1]
  ctau=sevdist$par[[3]][[5]][2]
  
  rx <- function(x, cmu, ctau) {
    (sevdist$par[[1]][[1]](x) - sevdist$par[[2]][[1]](x)) * 
      atan((x - cmu)/ctau)
  }
  
  r = try(pracma::quadinf(function(x){rx(x, cmu = cmu, ctau = ctau)}, xa = 0,xb = Inf)$Q)
  
  p=NULL

    rx_body <- function(x, cmu, ctau) {
      (1 - pcauchy(x, cmu, ctau)) * sevdist$par[[2]][[1]](x)
    }
    rx_tail <- function(x, cmu, ctau) {
      pcauchy(x, cmu, ctau) * sevdist$par[[1]][[1]](x)
    }
    
    r1 = try(pracma::integral(function(x){rx_body(x,cmu = cmu, ctau = ctau)}, xmin = 0, xmax = q))
    
    r2 = try(pracma::integral(function(x){rx_tail(x, cmu = cmu, ctau = ctau)}, xmin = 0, xmax = q))
    
    p=c(p,(r1+r2)/(1+r/pi))
    
  if (!lower.tail)  p = 1 - p
    
  return(p)
  })
}

qmixing<-function (p, sevdist, lower.tail = TRUE) 
{
  sapply(p, function(p){
  if (!lower.tail) p = 1 - p

    q = c(q, uniroot(function(x) {pmixing(x,sevdist)-p}, c(sevdist$par[[2]][[3]](p),sevdist$par[[1]][[3]](p)), extendInt = "yes")$root)

  return(as.numeric(q))
  })
}

rmixing = function(n,sevdist){
  rslt=NULL
  while(length(rslt)<n){
    u=rbinom(1,1,.5)
    x=u*sevdist$par[[2]][[4]](1)+(1-u)*sevdist$par[[1]][[4]](1)
    p=sevdist$par[[3]][[2]](x)
    if(runif(1)<u*p+(1-u)*p) rslt=c(rslt,x)
  }
  return(rslt)
}
