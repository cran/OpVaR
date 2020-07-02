fitWeights <- function(cell, thresh){
  n.body = length(which(cell$Loss <= thresh))
  n.total = length(cell$Loss)
  return(n.body/n.total)
}