pgpd <- function(q, loc, scale, shape){
  sapply(q, function(q){
    #define support 
    if((q < loc) || (shape < 0 && q > loc - scale/shape)){
      p = 0
    }else if(shape != 0){
      p = 1 - (1 + shape * (q - loc)/scale)^(-1/shape)
    }else if(shape == 0){
      p = 1- exp(-(q - loc)/scale)
    }
    return(p)
  })
}


dgpd=function (x, loc, scale, shape, log=FALSE){
  sapply(x, function(x) {
    if(!log){
      if ((x < loc) || (shape < 0 && x > loc - scale/shape)) {
        p = 0
      }
      else {
        if(shape!=0){
          p = 1/scale * (1 + shape * (x - loc)/scale)^(-1/shape - 1)
        }else{
          p=1/scale * exp(- (x - loc)/ scale)
        }
      }
    }else{
      if ((x < loc) || (shape < 0 && x > loc - scale/shape)) {
        p = -Inf
      }
      else {
        if(shape!=0){
          p = -log(scale) + (-1/shape - 1)*log(1 + shape * (x - loc)/scale)
        }else{
          p = -log(scale) - (x - loc)/ scale
        }
      }
    }
    return(p)
  })
}



rgpd=function (n, loc, scale, shape)
{
  return(replicate(n,ifelse(shape == 0, loc - scale * log(runif(n, 0, 1)),
                            loc + scale * (runif(n, 0, 1)^(-shape) - 1)/shape)))
}


qgpd <- function(p, loc, scale, shape){
  sapply(p, function(p){
    ifelse(shape == 0, loc - scale * log(1-p), loc + scale*((1-p)^(-shape)-1)/shape)
  })
}
