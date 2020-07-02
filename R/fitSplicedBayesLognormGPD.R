fitSplicedBayesLognormGPD<-function(cell,prior,burnin=10,niter=100,proposal_scale=evmix::flognormgpd(cell,method="Nelder-Mead")$se,
                          start=evmix::flognormgpd(cell,method="Nelder-Mead")$optim$par){

  # initialization of vectors where values will be assigned later
  # burnin length of the burn-in-phase, niter+burnin iterations -> assign values to function
  Sample_xi=numeric(niter+burnin)
  Sample_tau=numeric(niter+burnin)
  Sample_beta=numeric(niter+burnin)
  Sample_mu=numeric(niter+burnin)
  Sample_sigma=numeric(niter+burnin)
  
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
  
  prior_mu<-function(x){
    fun=dnorm(x,prior$mu[1],prior$mu[2])
    return(fun)
  }
  
  prior_sigma<-function(x){
    fun=truncnorm::dtruncnorm(x,0,Inf,prior$sigma[1],prior$sigma[2])
    return(fun)
  }
  
 
   
  ## log-likelihood function for the distribution model
  loglikelihood_lognormgpd<-function(x,vxi,vtau,vbeta,vmu,vsigma){ 
    likeli=sum(evmix::dlognormgpd(x,vmu,vsigma,vtau,vbeta,vxi,log=TRUE))
    return(likeli)}
 
   
  ## logarithmic acceptance rate for the symmetric proposal densities  
  acceptance_ratesymm<-function(vtheta_prop,vtheta_old,vlogposterior){
    a=min(log(1),(vlogposterior(vtheta_prop)-vlogposterior(vtheta_old)))
    return(a)
  }
  
  #### MH-step 1: sampling xi
  
  mhstep1<-function(vinitialval,cell){ 
    
    xi=vinitialval[1]  # current parameter values
    tau=vinitialval[2]
    beta=vinitialval[3]
    mu=vinitialval[4]
    sigma=vinitialval[5]
    
    xi_old=xi # current value for xi
    
    ## log-likelihood function depending on xi
    loglikelihoodfunction_xi<-function(x){ 
      fun=loglikelihood_lognormgpd(cell,x,tau,beta,mu,sigma)
      return(fun)
    }
    
    ## logarithmic function proportional to posterior density
    logposterior_xi<-function(x){ 
      fun=log(prior_xi(x))+loglikelihoodfunction_xi(x)
      return(fun)
    }
    
    ## parameter for the proposal density
    m=xi_old
    v=proposal_scale[5] # was assigned to the function

    xi_prop=rnorm(1,m,v)  # proposed value for the next iteration
    

    ############### comparing the logarithmic acceptance rate with u
    
    acc=acceptance_ratesymm(xi_prop,xi_old,logposterior_xi) 
    
    if(is.nan(acc)){acc=-Inf            # if NaN occurs, then refuse
                    }    # this could happen e.g. for a bad starting value with vlogposterior=-Inf
    
    ## comparing the logarithmic acceptance rate with the logarithm of a uniform distributed rv    
    u=log(runif(1))
    if (u<acc){
      xi_old=xi_prop # accept if u<acc, else keep old values
    }
    
    return(xi_old)
  }
 
  
   
  #### MH-step 2: sampling tau
 
   
  mhstep2<-function(vinitialval,cell){ 
    
    xi=vinitialval[1]  # current parameter values
    tau=vinitialval[2]
    beta=vinitialval[3]
    mu=vinitialval[4]
    sigma=vinitialval[5]
    
    tau_old=tau
    
    loglikelihoodfunction_tau<-function(x){ # log-likelihood function depending on tau
      fun=loglikelihood_lognormgpd(cell,xi,x,beta,mu,sigma)
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
    
    return(tau_old)
  }
  
  #### MH-step 3: sampling beta 
  
  mhstep3<-function(vinitialval,cell){ 
    
    xi=vinitialval[1]  # current parameter values
    tau=vinitialval[2]
    beta=vinitialval[3]
    mu=vinitialval[4]
    sigma=vinitialval[5]
    
    beta_old=beta
    
    loglikelihoodfunction_beta<-function(x){
      fun=loglikelihood_lognormgpd(cell,xi,tau,x,mu,sigma)
      return(fun)
    }
    
    logposterior_beta<-function(x){
      fun=log(prior_beta(x))+loglikelihoodfunction_beta(x)
      return(fun)
    }
    
    m=log(beta_old)
    v=proposal_scale[4]
    
    # log-normal proposal density
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
    
    return(beta_old)
  }
  
  
  #### MH-step 4: sampling mu
  
  mhstep4<-function(vinitialval,cell){ 
    
    xi=vinitialval[1]  # current parameter values
    tau=vinitialval[2]
    beta=vinitialval[3]
    mu=vinitialval[4]
    sigma=vinitialval[5]
    
    mu_old=mu
    
    
    loglikelihoodfunction_mu<-function(x){
      fun=loglikelihood_lognormgpd(cell,xi,tau,beta,x,sigma)
      return(fun)
    }
    
    logposterior_mu<-function(x){
      fun=log(prior_mu(x))+loglikelihoodfunction_mu(x)
      return(fun)
    }
    
    m=mu_old
    v=proposal_scale[1]
    
    
    logproposal<-function(x){
      fun=dnorm(x,m,v,log=TRUE)
      return(fun)
    }
    
    mu_prop=rnorm(1,m,v)
    
    #symmetric proposal density
    acc=acceptance_ratesymm(mu_prop,mu_old,logposterior_mu) 
    if(is.nan(acc)){acc=-Inf}
    
    u=log(runif(1))
    if (u<acc){
      mu_old=mu_prop
    }
    
    return(mu_old)
  }
  
  #### MH-step 5: sampling sigma
  
  mhstep5<-function(vinitialval,cell){ 
    
    xi=vinitialval[1] 
    tau=vinitialval[2]
    beta=vinitialval[3]
    mu=vinitialval[4]
    sigma=vinitialval[5]
    
    sigma_old=sigma
    
    loglikelihoodfunction_sigma<-function(x){
      fun=loglikelihood_lognormgpd(cell,xi,tau,beta,mu,x)
      return(fun)
    }
    
    logposterior_sigma<-function(x){
      fun=log(prior_sigma(x))+loglikelihoodfunction_sigma(x)
      return(fun)
    }
    
    m=log(sigma_old)
    v=proposal_scale[2]
    
    # log-normal proposal density
    logproposal<-function(x,m){
      fun=dlnorm(x,m,v,log=TRUE)
      return(fun)
    }
    
    sigma_prop=rlnorm(1,m,v)
    
    acc=min(log(1),(logposterior_sigma(sigma_prop)+logproposal(sigma_old,log(sigma_prop))-logposterior_sigma(sigma_old)-logproposal(sigma_prop,log(sigma_old))))
    if(is.nan(acc)){acc=-Inf}
    
    u=log(runif(1))
    if (u<acc){
      sigma_old=sigma_prop
    }
    
    return(sigma_old)
  }
  
  
  ##########################################################################################################
  # ACTUAL ALGORITHM
  ##########################################################################################################
 
    # starting values
    
    
    ########################################################################################################################
    
    initialval=c(start[5],start[3],start[4],start[1],start[2]) 
    
    
    ## actual MH-algorithm with number of iterations = niter + burnin
    # current process of updating of the vector initialval
    # cell is data given to the function
    
    for(i in seq(1:(niter+burnin))){
      initialval[1] = mhstep1(initialval, cell)
      Sample_xi[i] = initialval[1]
      
      initialval[2] = mhstep2(initialval, cell)
      Sample_tau[i] = initialval[2]
      
      initialval[3] = mhstep3(initialval, cell)
      Sample_beta[i] = initialval[3]
      
      initialval[4] = mhstep4(initialval, cell)
      Sample_mu[i] = initialval[4]
      
      initialval[5] = mhstep5(initialval, cell)
      Sample_sigma[i] = initialval[5]
 
    }
    
    
    # determining and saving sample values without burn-in-phase
    xi_estimator=mean(Sample_xi[(burnin+1):(niter+burnin)])
    tau_estimator=mean(Sample_tau[(burnin+1):(niter+burnin)])
    beta_estimator=mean(Sample_beta[(burnin+1):(niter+burnin)])
    mu_estimator=mean(Sample_mu[(burnin+1):(niter+burnin)])
    sigma_estimator=mean(Sample_sigma[(burnin+1):(niter+burnin)])

  
  # final output
  buildSplicedSevdist("lnorm",c(mu_estimator,sigma_estimator),"gpd",c(tau_estimator,beta_estimator,xi_estimator),tau_estimator,0.5)
}