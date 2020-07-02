fitSplicedBestFit <- function(cell, body, tail, thresh0 = 0.7, thresh.max = 0.98){
  N = length(cell$Loss)
  loss = sort(cell$Loss)
  thresh0 = as.numeric(quantile(cell$Loss, thresh0))
  thresh.max = as.numeric(quantile(cell$Loss, thresh.max))
  lower = max(which(loss<=thresh0))
  upper = max(which(loss<=thresh.max))
  results = list()
  test = c(rep(0, (upper-lower+1)))
  count = 1
  for(i in seq(lower, upper)){
    par <- fitSplicedPar(cell, loss[i], body, tail)
    weight <- fitWeights(cell, loss[i])
    sevdist.tmp <- buildSplicedSevdist(body, par[[1]], tail, par[[2]], loss[i], weight)
    pval <- suppressWarnings(stats::ks.test(loss, rspliced(N, sevdist.tmp))$p.value)
    results[[count]] = c(par, loss[i], weight)
    test[count] = pval
    count = count + 1 
  }
  loc = min(which(test == max(test[which(is.na(test) == FALSE)])))
  sevdist = buildSplicedSevdist(body, results[[loc]][[1]], tail, results[[loc]][[2]], results[[loc]][[3]], results[[loc]][[4]])
  return(sevdist)
}
