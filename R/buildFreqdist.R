buildFreqdist<-
function(freq.distr, freq.param){
  freqdist = list()
  freqdist[[1]] = freq.distr
  freqdist[[2]] = freq.param
  class(freqdist)="freqdist"
  return(freqdist)
}