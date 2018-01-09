
fitSpliced=function(cell, body, tail, method,thresh=NULL){
  if(method == "BestFit"){
    return(fitSplicedBestFit(cell, body, tail, thresh0 = 0.7, thresh.max = 0.98))
  } else {
    if(method != "Fixed"){
      thresh<-fitThreshold(cell, body, tail, method) 
    }
    pars=fitSplicedPar(cell, thresh, body, tail) 
    weight=fitWeights(cell,thresh)
    sevdist=buildSplicedSevdist(body, pars[[1]], tail, pars[[2]], thresh, weight)
    return(sevdist)
  }
}