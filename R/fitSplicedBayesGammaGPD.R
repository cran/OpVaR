fitSplicedBayesGammaGPD<-function(cell,prior,burnin=10,niter=100,proposal_scale=evmix::fgammagpd(cell,method="Nelder-Mead")$se,
                            start=evmix::fgammagpd(cell,method="Nelder-Mead")$optim$par){
   
  # initialization of vectors where values will be assigned later  
  # burnin length of the burn-in-phase, niter+burnin iterations -> assign values to function
  Sample_xi=numeric(niter+burnin)
  Sample_tau=numeric(niter+burnin)
  Sample_beta=numeric(niter+burnin)
  Sample_gscale=numeric(niter+burnin)
  Sample_gshape=numeric(niter+burnin)
  
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
  
  prior_gscale<-function(x){
    fun=truncnorm::dtruncnorm(x,0,Inf,prior$gscale[1],prior$gscale[2])
    return(fun)
  }
  
  prior_gshape<-function(x){
    fun=truncnorm::dtruncnorm(x,0,Inf,prior$gshape[1],prior$gshape[2])
    return(fun)
  }
  
  
  
  ## log-likelihood function for the distribution model
  loglikelihood_gammagpd<-function(x,vgshape,vgscale,vtau,vbeta,vxi){
    likeli=sum(evmix::dgammagpd(x,vgshape,vgscale,vtau,vbeta,vxi,log=TRUE))
    return(likeli)}
  
  ## logarithmic acceptance rate for the symmetric proposal densities
  acceptance_ratesymm<-function(vtheta_prop,vtheta_old,vlogposterior){
    a=min(log(1),(vlogposterior(vtheta_prop)-vlogposterior(vtheta_old)))
    return(a)
  }
  
  
  #### MH-step 1: sampling xi
  
  mhstep1<-function(vinitialval,cell){ 
    
    gshape=vinitialval[1]  # current parameter values
    gscale=vinitialval[2]
    tau=vinitialval[3]
    xi=vinitialval[4]
    beta=vinitialval[5]
    
    xi_old=xi # current value for xi
    
    ## log-likelihood function depending on xi
    loglikelihoodfunction_xi<-function(x){ 
      fun=loglikelihood_gammagpd(cell,gshape,gscale,tau,beta,x)
      return(fun)
    }
    
    ## logarithmic function proportional to posterior density
    logposterior_xi<-function(x){ 
      fun=log(prior_xi(x))+loglikelihoodfunction_xi(x)
      return(fun)
    }
    
    ## parameter for the proposal density
    m=xi_old
    v=proposal_scale[1] # was assigned to the function
  
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
    
    gshape=vinitialval[1]  # current parameter values
    gscale=vinitialval[2]
    tau=vinitialval[3]
    xi=vinitialval[4]
    beta=vinitialval[5]
    
    tau_old=tau
    
    loglikelihoodfunction_tau<-function(x){ # log-likelihood function depending on tau
      fun=loglikelihood_gammagpd(cell,gshape,gscale,x,beta,xi)
      return(fun)
    }
    
    logposterior_tau<-function(x){
      fun=log(prior_tau(x))+loglikelihoodfunction_tau(x)
      return(fun)
    }
    
    # truncated normal distribution as proposal density
    m=tau_old
    v=proposal_scale[2]
    
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
    
    gshape=vinitialval[1]  # current parameter values
    gscale=vinitialval[2]
    tau=vinitialval[3]
    xi=vinitialval[4]
    beta=vinitialval[5]
    
    beta_old=beta
    
    loglikelihoodfunction_beta<-function(x){
      fun=loglikelihood_gammagpd(cell,gshape,gscale,tau,x,xi)
      return(fun)
    }
    
    logposterior_beta<-function(x){
      fun=log(prior_beta(x))+loglikelihoodfunction_beta(x)
      return(fun)
    }
    
    # log-normal proposal density
    m=log(beta_old)
    v=proposal_scale[3]
    
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
  
  
  #### MH-step 4: sampling gshape
  
  
  mhstep4<-function(vinitialval,cell){ 
    
    gshape=vinitialval[1]  # current parameter values
    gscale=vinitialval[2]
    tau=vinitialval[3]
    xi=vinitialval[4]
    beta=vinitialval[5]
    
    gshape_old=gshape 
    
    loglikelihoodfunction_gshape<-function(x){ 
      fun=loglikelihood_gammagpd(cell,x,gscale,tau,beta,xi)
      return(fun)
    }
    
    logposterior_gshape<-function(x){ 
      fun=log(prior_gshape(x))+loglikelihoodfunction_gshape(x)
      return(fun)
    }
    
    
    m=log(gshape_old)
    v=proposal_scale[4]
    
    logproposal<-function(x,m){
      fun=dlnorm(x,m,v,log=TRUE)
      return(fun)
    }
    
    
    gshape_prop=rlnorm(1,m,v)
    
    acc=min(log(1),(logposterior_gshape(gshape_prop)+logproposal(gshape_old,log(gshape_prop))-logposterior_gshape(gshape_old)-logproposal(gshape_prop,log(gshape_old))))
    if(is.nan(acc)){acc=-Inf
                    } 
    
    
    u=log(runif(1))
    if (u<acc){
      
      gshape_old=gshape_prop
    }
    
    return(gshape_old)
  }
  
  
  #### MH-step 5: sampling gscale
  
  mhstep5<-function(vinitialval,cell){ 
    
    gshape=vinitialval[1]  
    gscale=vinitialval[2]
    tau=vinitialval[3]
    xi=vinitialval[4]
    beta=vinitialval[5]
    
    gscale_old=gscale
    
    loglikelihoodfunction_gscale<-function(x){
      fun=loglikelihood_gammagpd(cell,gshape,x,tau,beta,xi)
      return(fun)
    }
    
    logposterior_gscale<-function(x){
      fun=log(prior_gscale(x))+loglikelihoodfunction_gscale(x)
      return(fun)
    }
    
    m=log(gscale_old)
    v=proposal_scale[5]
    logproposal<-function(x,m){
      fun=dlnorm(x,m,v,log=TRUE)
      return(fun)
    }
    
    gscale_prop=rlnorm(1,m,v)
    
    acc=min(log(1),(logposterior_gscale(gscale_prop)+logproposal(gscale_old,log(gscale_prop))-logposterior_gscale(gscale_old)-logproposal(gscale_prop,log(gscale_old))))
    if(is.nan(acc)){acc=-Inf
                    } 
    
    u=log(runif(1))
    if (u<acc){
      
      gscale_old=gscale_prop
    }
  
    return(gscale_old)
  }
  
  
  ##########################################################################################################
  # ACTUAL ALGORITHM
  ##########################################################################################################

    # starting values
    
    
    ########################################################################################################################
    
    initialval=c(start[1],start[2],start[3],start[5],start[4])
    
    
    ## actual MH-algorithm with number of iterations = niter + burnin
    # current process of updating of the vector initialval
    # cell is data given to the function
    
    for(i in seq(1:(niter+burnin))){
      initialval[4] = mhstep1(initialval,cell)
      Sample_xi[i] = initialval[4] 
      
      initialval[3] = mhstep2(initialval,cell)
      Sample_tau[i] = initialval[3] 
      
      initialval[5] = mhstep1(initialval,cell)
      Sample_beta[i] = initialval[5] 

      initialval[1] = mhstep1(initialval,cell)
      Sample_gshape[i] = initialval[1] 
      
      initialval[2] = mhstep1(initialval,cell)
      Sample_gscale[i] = initialval[2] 
    }

    
    # determining and saving sample values without burn-in-phase
    xi_estimator=mean(Sample_xi[(burnin+1):(niter+burnin)])
    tau_estimator=mean(Sample_tau[(burnin+1):(niter+burnin)])
    beta_estimator=mean(Sample_beta[(burnin+1):(niter+burnin)])
    gshape_estimator=mean(Sample_gshape[(burnin+1):(niter+burnin)])
    gscale_estimator=mean(Sample_gscale[(burnin+1):(niter+burnin)])
  
  
  # final output
  buildSplicedSevdist("gamma",c(gshape_estimator,gscale_estimator),"gpd",c(tau_estimator,beta_estimator,xi_estimator),tau_estimator,0.5)
}