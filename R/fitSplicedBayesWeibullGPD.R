fitSplicedBayesWeibullGPD<-function(cell,prior,burnin=10,niter=100,proposal_scale=evmix::fweibullgpd(cell,method="Nelder-Mead")$se,
                            start=evmix::fweibullgpd(cell,method="Nelder-Mead")$optim$par){
  
  # initialization of vectors of estimators where values will be assigned later
  # burnin length of the burn-in-phase, niter+burnin iterations -> assign values to function
  Sample_xi=numeric(niter+burnin)
  Sample_tau=numeric(niter+burnin)
  Sample_beta=numeric(niter+burnin)
  Sample_wscale=numeric(niter+burnin)
  Sample_wshape=numeric(niter+burnin)
  
  
  # prior distributions (usually quite flat)
  
  prior_xi<-function(x){ 
    fun=truncnorm::dtruncnorm(x,0,Inf,prior$xi[1],prior$xi[2])
    return(fun)
  }
  
  prior_tau<-function(x){ # preferably informative, otherwise too large distortion can occur
    fun=truncnorm::dtruncnorm(x,0,Inf,prior$tau[1],prior$tau[2])
    return(fun)
  }
  
  prior_beta<-function(x){
    fun=truncnorm::dtruncnorm(x,0,Inf,prior$beta[1],prior$beta[2])
    return(fun)
  } 
  
  prior_wshape<-function(x){
    fun=truncnorm::dtruncnorm(x,0,Inf,prior$wshape[1],prior$wshape[2])
    return(fun)
  }
  
  prior_wscale<-function(x){
    fun=truncnorm::dtruncnorm(x,0,Inf,prior$wscale[1],prior$wscale[2])
    return(fun)
  }

  
  
  ## log-likelihood function for the distribution model
  loglikelihood_weibullgpd<-function(x,vwshape,vwscale,vtau,vbeta,vxi){
    likeli=sum(evmix::dweibullgpd(x,vwshape,vwscale,vtau,vbeta,vxi,log=TRUE))
    return(likeli)}
  
  ## logarithmic acceptance rate for the symmetric proposal densities
  acceptance_ratesymm<-function(vtheta_prop,vtheta_old,vlogposterior){
    a=min(log(1),(vlogposterior(vtheta_prop)-vlogposterior(vtheta_old)))
    return(a)
  }
  
  
  #### MH-step 1: sampling xi
  
  mhstep1<-function(vinitialval,cell){ 
    
    wshape=vinitialval[1]  # current parameter values
    wscale=vinitialval[2]
    tau=vinitialval[3]
    xi=vinitialval[4]
    beta=vinitialval[5]
    
    xi_old=xi # current value for xi
    
    ## log-likelihood function depending on xi
    loglikelihoodfunction_xi<-function(x){ 
      fun=loglikelihood_weibullgpd(cell,wshape,wscale,tau,beta,x)
      return(fun)
    }
    
    ## logarithmic function proportional to posterior density
    logposterior_xi<-function(x){ 
      fun=log(prior_xi(x))+loglikelihoodfunction_xi(x)
      return(fun)
    }
    
    ## parameter for the proposal density
    m=xi_old
    v=proposal_scale[4] # was assigned to the function
    
    xi_prop=rnorm(1,m,v)  # proposed value for the next iteration
    
    
    ############### comparing the logarithmic acceptance rate with u
    
    acc=acceptance_ratesymm(xi_prop,xi_old,logposterior_xi)
    
    
    if(is.nan(acc)){acc=-Inf            # if NaN occurs, then refuse
                    } # this could happen e.g. for a bad starting value with vlogposterior=-Inf
    
    ## comparing the logarithmic acceptance rate with the logarithm of a uniform distributed rv
    u=log(runif(1))
    if (u<acc){
      xi_old=xi_prop # accept if u<acc, else keep old values
    }
    
    initialval[4]<<-xi_old  # saving the current values in the vector initialval
    Sample_xi[i]<<-xi_old # saving the vector
    return(xi_old) # returns the current value
  }
  
  
  
  #### MH-step 2: sampling tau

    
  mhstep2<-function(vinitialval,cell){ 
    
    wshape=vinitialval[1]  # current parameter values
    wscale=vinitialval[2]
    tau=vinitialval[3]
    xi=vinitialval[4]
    beta=vinitialval[5]
    
    tau_old=tau
    
    loglikelihoodfunction_tau<-function(x){ # log-likelihood function depending on tau
      fun=loglikelihood_weibullgpd(cell,wshape,wscale,x,beta,xi)
      return(fun)
    }
    
    logposterior_tau<-function(x){
      fun=log(prior_tau(x))+loglikelihoodfunction_tau(x)
      return(fun)
    }
    
    # truncated normal distribution as proposal density
    m=tau_old
    v=proposal_scale[3]
    
    # logarithmic proposal density for acceptance rate (since now not symmetric anymore)
    logproposal<-function(x,m){
      fun=log(truncnorm::dtruncnorm(x,min(cell),max(cell),m,v))
      return(fun)
    }
    
    # proposed value
    tau_prop=truncnorm::rtruncnorm(1,min(cell),max(cell),m,v)
    
    # logarithmic acceptance rate
    acc=min(log(1),(logposterior_tau(tau_prop)+logproposal(tau_old,tau_prop)-logposterior_tau(tau_old)-logproposal(tau_prop,tau_old)))
    if(is.nan(acc)){acc=-Inf
    }
    
    u=log(runif(1))
    if (u<acc){
      tau_old=tau_prop

    }
    
    initialval[3]<<-tau_old
    Sample_tau[i]<<-tau_old
    return(tau_old)
  }
  
  
  #### MH-step 3: sampling beta 

  mhstep3<-function(vinitialval,cell){ 
    
    wshape=vinitialval[1]  # current parameter values
    wscale=vinitialval[2]
    tau=vinitialval[3]
    xi=vinitialval[4]
    beta=vinitialval[5]
    
    beta_old=beta
    
    loglikelihoodfunction_beta<-function(x){
      fun=loglikelihood_weibullgpd(cell,wshape,wscale,tau,x,xi)
      return(fun)
    }
    
    logposterior_beta<-function(x){
      fun=log(prior_beta(x))+loglikelihoodfunction_beta(x)
      return(fun)
    }
    
    # log-normal proposal density
    m=log(beta_old)
    v=proposal_scale[5]
    
    logproposal<-function(x,m){
      fun=dlnorm(x,m,v,log=TRUE)
      return(fun)
    }
    
    # proposed value
    beta_prop=rlnorm(1,m,v)
    
    # logarithmic acceptance rate
    acc=min(log(1),(logposterior_beta(beta_prop)+logproposal(beta_old,log(beta_prop))-logposterior_beta(beta_old)-logproposal(beta_prop,log(beta_old))))
    if(is.nan(acc)){acc=-Inf
    }

    
    u=log(runif(1))  
    
    if(u<acc){
      beta_old=beta_prop

    }
    
    initialval[5]<<-beta_old
    Sample_beta[i]<<-beta_old
    return(beta_old)
  }
  
  
  #### MH-step 4: sampling wshape
  
  mhstep4<-function(vinitialval,cell){ 
    
    wshape=vinitialval[1]  # current parameter values
    wscale=vinitialval[2]
    tau=vinitialval[3]
    xi=vinitialval[4]
    beta=vinitialval[5]
    
    wshape_old=wshape 
    
    loglikelihoodfunction_wshape<-function(x){ 
      fun=loglikelihood_weibullgpd(cell,x,wscale,tau,beta,xi)
      return(fun)
    }
    
    logposterior_wshape<-function(x){ 
      fun=log(prior_wshape(x))+loglikelihoodfunction_wshape(x)
      return(fun)
    }
    
    
    m=log(wshape_old)
    v=proposal_scale[1]
    
    logproposal<-function(x,m){
      fun=dlnorm(x,m,v,log=TRUE)
      return(fun)
    }
    
    
    wshape_prop=rlnorm(1,m,v)
    
    acc=min(log(1),(logposterior_wshape(wshape_prop)+logproposal(wshape_old,log(wshape_prop))-logposterior_wshape(wshape_old)-logproposal(wshape_prop,log(wshape_old))))
    if(is.nan(acc)){acc=-Inf
    }
    
    
    u=log(runif(1))
    if (u<acc){
      
      wshape_old=wshape_prop

    }
    
    initialval[1]<<-wshape_old
    Sample_wshape[i]<<-wshape_old
    return(wshape_old)
  }
  
  
  #### MH-step 5: sampling wscale

  mhstep5<-function(vinitialval,cell){ 
    
    wshape=vinitialval[1]  
    wscale=vinitialval[2]
    tau=vinitialval[3]
    xi=vinitialval[4]
    beta=vinitialval[5]
    
    wscale_old=wscale
    
    loglikelihoodfunction_wscale<-function(x){
      fun=loglikelihood_weibullgpd(cell,wshape,x,tau,beta,xi)
      return(fun)
    }
    
    logposterior_wscale<-function(x){
      fun=log(prior_wscale(x))+loglikelihoodfunction_wscale(x)
      return(fun)
    }
    
    m=log(wscale_old)
    v=proposal_scale[2]
    
    logproposal<-function(x,m){ 
      fun=dlnorm(x,m,v,log=TRUE)
      return(fun)
    }
    
    
    wscale_prop=rlnorm(1,m,v)
    
    acc=min(log(1),(logposterior_wscale(wscale_prop)+logproposal(wscale_old,log(wscale_prop))-logposterior_wscale(wscale_old)-logproposal(wscale_prop,log(wscale_old))))
    if(is.nan(acc)){acc=-Inf
    }
    
    
    u=log(runif(1))
    if (u<acc){
      
      gscale_old=wscale_prop

    }
    
    initialval[2]<<-wscale_old
    Sample_wscale[i]<<-wscale_old
    return(wscale_old)
  }
  
  
  ##########################################################################################################
  # ACTUAL ALGORITHM
  ##########################################################################################################
    
    # starting values
    
    
    ########################################################################################################################
    
    initialval=c(start[1],start[2],start[3],start[5],start[4])
    
    ########################################################################################################################

    set.seed(1)
    
    ## actual MH-algorithm with number of iterations = niter + burnin
    # current process of updating of the vector initialval
    # cell is data given to the function
    
    for(i in seq(1:(niter+burnin))){
      mhstep1(initialval,cell)
      mhstep2(initialval,cell)
      mhstep3(initialval,cell)
      mhstep4(initialval,cell)
      mhstep5(initialval,cell)
    }
    
    # determining and saving sample values without burn-in-phase
    xi_estimator=mean(Sample_xi[(burnin+1):(niter+burnin)])
    tau_estimator=mean(Sample_tau[(burnin+1):(niter+burnin)])
    beta_estimator=mean(Sample_beta[(burnin+1):(niter+burnin)])
    wshape_estimator=mean(Sample_wshape[(burnin+1):(niter+burnin)])
    wscale_estimator=mean(Sample_wscale[(burnin+1):(niter+burnin)])
  
  
  # final output
  buildSplicedSevdist("weibull",c(wshape_estimator,wscale_estimator),"gpd",c(tau_estimator,beta_estimator,xi_estimator),tau_estimator,0.5)
}
