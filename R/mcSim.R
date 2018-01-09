
mcSim=function (opriskmodel, n_sim, verbose=TRUE) 
{
  total_loss = list()
  for (i in 1:length(opriskmodel)) {
    if(verbose){
      print(paste("Cell", i))
      pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
    }
    total_loss_sim_cell = NULL
    for (iter in 1:n_sim) {
      sevmodel = opriskmodel[[i]]$sevdist
      sevsim = 0
      if(opriskmodel[[i]]$freqdist[[1]]=="pois"){
        n=rpois(1,opriskmodel[[i]]$freqdist[[2]])
      }else if(opriskmodel[[i]]$freqdist[[1]]=="nbinom"){
        n=rnbinom(1,opriskmodel[[i]]$freqdist[[2]][1],opriskmodel[[i]]$freqdist[[2]][2])
      }
      
      if (is.null(opriskmodel[[i]]$dependency)||opriskmodel[[i]]$dependency$family==0){
        sevsim = rsevdist(n, sevmodel)
      }
      else {
        u1 = eval(parse(text = paste0("p", opriskmodel[[i]]$freqdist[[1]], 
                                      "(n,", paste0(opriskmodel[[i]]$freqdist[[2]], 
                                                    collapse = ","), ")")))
        u2 = BiCopCondSim(n,cond.val=u1,cond.var=1,opriskmodel[[i]]$dependency)
        sevsim = qsevdist(u2, sevmodel)
      }
      total_loss_sim_cell = c(total_loss_sim_cell, sum(sevsim))
      if(verbose){
        setTxtProgressBar(pb, iter)
      }
    }
    total_loss[[i]] = total_loss_sim_cell
  }
  total_loss
}


VaR=function(mc_out,alpha){
  varvec=NULL
  for(i in 1:length(mc_out)){
    varvec=c(varvec,quantile(mc_out[[i]],alpha))
  }
  varvec
}
