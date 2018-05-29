evaluate <- function(d, distr, param, z){
  return(do.call(paste(d, distr, sep=""),c(list(z), param)))
}

buildPlainSevdist <- function(distr, param){
  sevdist=list()
  sevdist$type = "plain"
  parlist=list()
  parlist[[1]]=list()
  parlist[[1]]=function(x) evaluate("d", distr, param, x) 
  parlist[[2]]=function(q) evaluate("p", distr, param, q) 
  parlist[[3]]=function(p) evaluate("q", distr, param, p) 
  parlist[[4]]=function(n) evaluate("r", distr, param, n) 
  parlist[[5]]=param
  parlist[[6]]=distr
  
  sevdist$par=parlist
  class(sevdist)="sevdist"
  return(sevdist) 
}


buildMixingSevdist=
  function(body.distr, body.param, tail.distr, tail.param, mixing.param){
    sevdist=list()
    sevdist$type = "mixing"
    parlist=list()
    parlist[[1]]=list()
    parlist[[2]]=list()
    parlist[[3]]=list()
    
    parlist[[1]][[1]]=function(x){evaluate("d", tail.distr, tail.param, x)}
    parlist[[1]][[2]]=function(q){evaluate("p", tail.distr, tail.param, q)}
    parlist[[1]][[3]]=function(p){evaluate("q", tail.distr, tail.param, p)}
    parlist[[1]][[4]]=function(n){evaluate("r", tail.distr, tail.param, n)}
    parlist[[1]][[5]]=tail.param
    parlist[[1]][[6]]=tail.distr
    
    parlist[[2]][[1]]=function(x){evaluate("d", body.distr, body.param, x)}
    parlist[[2]][[2]]=function(q){evaluate("p", body.distr, body.param, q)}
    parlist[[2]][[3]]=function(p){evaluate("q", body.distr, body.param, p)}
    parlist[[2]][[4]]=function(n){evaluate("r", body.distr, body.param, n)}
    parlist[[2]][[5]]=body.param
    parlist[[2]][[6]]=body.distr
    
    parlist[[3]][[1]]=function(x){evaluate("d", "cauchy", mixing.param, x)}
    parlist[[3]][[2]]=function(q){evaluate("p", "cauchy", mixing.param, q)}
    parlist[[3]][[3]]=function(p){evaluate("q", "cauchy", mixing.param, p)}
    parlist[[3]][[4]]=function(n){evaluate("r", "cauchy", mixing.param, n)}
    parlist[[3]][[5]]=mixing.param
    parlist[[3]][[6]]="cauchy"
    
    sevdist$par=parlist
    class(sevdist)="sevdist"
    return(sevdist)
  }
buildSplicedSevdist <- function(body.distr, body.param, tail.distr, tail.param, thresh, weight){
  sevdist=list()
  sevdist$type = "spliced"
  parlist=list()
  parlist[[1]]=list()
  parlist[[2]]=list()
  
  parlist[[1]][[1]]=function(x){ifelse(x>thresh, evaluate("d", tail.distr, tail.param, x)/(1-evaluate("p", tail.distr, tail.param, thresh)), 0)}
  parlist[[1]][[2]]=function(q){ifelse(q>thresh, (evaluate("p", tail.distr, tail.param, q) -  evaluate("p", tail.distr, tail.param, thresh)) /(1-evaluate("p", tail.distr, tail.param, thresh)),0)}
  parlist[[1]][[3]]=function(p){ifelse(p>weight, evaluate("q", tail.distr, tail.param, evaluate("p", tail.distr, tail.param, thresh) + p*(1-evaluate("p", tail.distr, tail.param, thresh))),0)}
  parlist[[1]][[5]]=tail.param
  parlist[[1]][[6]]=tail.distr
  
  parlist[[2]][[1]]=function(x){ifelse(x<=thresh, evaluate("d", body.distr, body.param, x)/evaluate("p", body.distr, body.param, thresh), 0)}
  parlist[[2]][[2]]=function(q){ifelse(q<=thresh, evaluate("p", body.distr, body.param, q)/evaluate("p", body.distr, body.param, thresh),0)}
  parlist[[2]][[3]]=function(p){ifelse(p<=weight, evaluate("q", body.distr, body.param, p*evaluate("p", body.distr, body.param, thresh)), 0)}
  parlist[[2]][[5]]=body.param
  parlist[[2]][[6]]=body.distr
  
  sevdist$par=parlist
  sevdist$thresh = thresh
  sevdist$weights = c(weight, 1-weight)
  class(sevdist)="sevdist"
  return(sevdist) 
}