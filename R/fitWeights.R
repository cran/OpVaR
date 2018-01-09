fitWeights <- function(cell, thresh){
  n.body = length(which(cell$Loss <= thresh))
  n.tail = length(which(cell$Loss > thresh))
  return(1/(n.tail/n.body +1))
}