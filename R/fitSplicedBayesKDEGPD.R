fitSplicedBayesKDEGPD<-function(cell,prior,burnin=10,niter=100,proposal_scale=evmix::flognormgpd(cell,method="Nelder-Mead")$se,
                          start=evmix::flognormgpd(cell,method="Nelder-Mead")$optim$par){
  
  # initialization of vectors where values will be assigned later
  # burnin length of the burn-in-phase, niter+burnin iterations -> assign values to function
  Sample_xi=numeric(niter+burnin)
  Sample_tau=numeric(niter+burnin)
  Sample_beta=numeric(niter+burnin)
  
  
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
  
  ## set bandwith based on the density
  bandwidth=density(log(cell),kernel="gaussian")$bw 
  
  ## gaussian kernel density function with boundary correction
  gkernellog<-function(x,vdata,h){ 
    fun=(1/length(vdata))*sum(dlnorm(x,log(vdata),h))
    return(fun)}
  
  ## distribution function
  Fgkernellog<-function(x,vdata,h){
    fun=(1/length(vdata))*sum(plnorm(x,log(vdata),h))
    return(fun)}
  
  ## combound density: kernel density + GPD
  dkerneldensitygpd=function(x,vtau,vbeta,vxi){ 
    if(x<=vtau){
      return(gkernellog(x,cell,bandwidth))
    }else{
      return((1-Fgkernellog(vtau,cell,bandwidth))*evmix::dgpd(x,vtau,vbeta,vxi))
    }
  }
  
  ## vectorized
  vdkerneldensitygpd=Vectorize(dkerneldensitygpd) 
  
 
  ## log-likelihood function for the distribution model
  loglikelihood_kerneldensitygpd<-function(vdata,vtau,vbeta,vxi){
    likeli=sum(log(vdkerneldensitygpd(vdata,vtau,vbeta,vxi)))
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
    
    xi_old=xi # current value for xi
    
    ## log-likelihood function depending on xi
    loglikelihoodfunction_xi<-function(x){ 
      fun=loglikelihood_kerneldensitygpd(cell,tau,beta,x)
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
    

    acc=acceptance_ratesymm(xi_prop,xi_old,logposterior_xi) 
    
    if(is.nan(acc)){acc=-Inf            # if NaN occurs, then refuse
                    }
    
    ## comparing the logarithmic acceptance rate with the logarithm of a uniform distributed rv   
    u=log(runif(1))
    if (u<acc){
      xi_old=xi_prop # accept if u<acc, else keep old values
    }

    return(xi_old) # returns the current value
  }
  
 
   
  #### MH-step 2: sampling tau
  
  
  mhstep2<-function(vinitialval,cell){ 
    
    xi=vinitialval[1]  # current parameter values
    tau=vinitialval[2]
    beta=vinitialval[3]
    
    tau_old=tau
    
    loglikelihoodfunction_tau<-function(x){ # log-likelihood function depending on tau
      fun=loglikelihood_kerneldensitygpd(cell,x,beta,xi)
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
    
    beta_old=beta
    
    loglikelihoodfunction_beta<-function(x){
      fun=loglikelihood_kerneldensitygpd(cell,tau,x,xi)
      return(fun)
    }
    
    logposterior_beta<-function(x){
      fun=log(prior_beta(x))+loglikelihoodfunction_beta(x)
      return(fun)
    }
    
    # log-normal proposal density
    m=log(beta_old)
    v=proposal_scale[4]
    
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
  
  
  ##########################################################################################################
  # ACTUAL ALGORITHM
  ##########################################################################################################

    #initial values determined by flognormgpd

    initialval=c(start[5],start[3],start[4])
    

    ## actual MH-algorithm with number of iterations = niter + burnin
    # current process of updating of the vector initialval
    # cell is data given to the function
    
    for(i in seq(1:(niter+burnin))){
      initialval[1] =  mhstep1(initialval, cell)
      Sample_xi[i] =  initialval[1]
      
      initialval[2] = mhstep2(initialval, cell)
      Sample_tau[i] = initialval[2]
      
      initialval[3] = mhstep3(initialval,cell)
      Sample_beta[i] = initialval[3] 
      
    }
    
    
    # determining and saving sample values without burn-in-phase
    xi_estimator=mean(Sample_xi[(burnin+1):(niter+burnin)])
    tau_estimator=mean(Sample_tau[(burnin+1):(niter+burnin)])
    beta_estimator=mean(Sample_beta[(burnin+1):(niter+burnin)])
  
    
  # final output
  buildSplicedSevdist("lnorm",c(start[1],start[2]),"gpd",c(tau_estimator,beta_estimator,xi_estimator),tau_estimator,0.5) 
}