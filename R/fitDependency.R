fitDependency=function(cell,copula_family,method="mle"){
  freq_x=as.numeric(table(cell$Period))
  freq_u=pobs(rep(freq_x,freq_x))
  sev_u=pobs(cell$Loss)
  copfit=BiCopEst(freq_u,sev_u,copula_family,method=method)
  copfit
}

