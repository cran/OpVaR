fitFreqdist=function (cell, distr) 
{
  distr2 <- switch(distr, pois = "poisson", nbinom = "nbinomial")
  if (is.null(distr2)) {
    stop("The distribution has either to be 'pois' or 'nbinom'!")
  }
  periods = unique(cell$Period)
  count = stats::ftable(c(periods, cell$Period), col.vars = 1, 
                        row.vars = NULL)
  param = as.numeric(vcd::goodfit(as.vector(count), type = distr2, 
                                  method = "ML")$par)
  return(buildFreqdist(distr,param))
}