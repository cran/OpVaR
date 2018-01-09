fitMixing = function(cell,body,tail,method="L-BFGS-B",c_location0=0.75,c_scale0=2){
  
  dat=cell$Loss
  
  # Initialization Body Distribution
  
  mu=mean(dat)
  sigma_sq=var(dat)
  sd=sd(dat)
  
  if(body=="lnorm"){
    par0_body=c(mean(log(dat)),log(sd(log(dat))))
  }else if(body=="weibull"){
    par0_body=log(c((sd/mu)^(-1.086),mu/(gamma(1+1/((sd/mu)^(-1.086))))))
  }else if(body=="gamma"){
    par0_body=log(c(mu^2/sigma_sq,mu/sigma_sq))
  }
  lower_body=-10*abs(par0_body)
  upper_body=10*abs(par0_body)
  
  # Initializtion Tail Distribution
  if(tail=="lnorm"){
    par0_tail=c(mean(log(dat)),log(sd(log(dat))))
  }else if(tail=="weibull"){
    par0_tail=log(c((sd/mu)^(-1.086),mu/(gamma(1+1/((sd/mu)^(-1.086))))))
  }else if(tail=="gamma"){
    par0_tail=log(c(mu^2/sigma_sq,mu/sigma_sq))
  }else if(tail=="gpd"){
    par0_tail=suppressWarnings(evmix::fgpd(dat))$mle
    par0_tail=c(log(par0_tail[1]),par0_tail[2])
  }
  lower_tail=-5*abs(par0_tail)
  upper_tail=5*abs(par0_tail)
  if(tail=="gpd") lower_tail[1]=0
  
  par0_cauchy=c(quantile(dat,c_location0),log(sd/c_scale0))
  
  
  lower_cauchy=c(.5*abs(par0_cauchy[1]),-5*abs(par0_cauchy[2]))
  upper_cauchy=c(1.5*abs(par0_cauchy[1]),5*abs(par0_cauchy[2]))
  
  nll=function(x){
    p=c(x[1],exp(x[2]))
    if(body=="lnorm"){
      l1=c(x[3],exp(x[4]))
    }else if(body=="weibull"){
      l1=c(exp(x[3]),exp(x[4]))
    }else if(body=="gamma"){
      l1=c(exp(x[3]),exp(x[4]))
    }
    
    if(tail=="lnorm"){
      l2=c(x[5],exp(x[6]))
    }else if(tail=="weibull"){
      l2=c(exp(x[5]),exp(x[6]))
    }else if(tail=="gamma"){
      l2=c(exp(x[5]),exp(x[6]))
    }else if(tail=="gpd"){
      l2=c(0,exp(x[5]),x[6])
    }
    sev_tmp=buildMixingSevdist(body,l1,tail,l2,p)
    val=-sum(dmixing(dat,sev_tmp,log=TRUE))
    return(val)
  }
  
  par0=c(par0_cauchy,par0_body,par0_tail)
  lower0=c(lower_cauchy,lower_body,lower_tail)
  upper0=c(upper_cauchy,upper_body,upper_tail)

  out=suppressWarnings(optim(par0,nll,method=method,lower=lower0,upper=upper0))$par
  
  par_cauchy=as.numeric(c(out[1],exp(out[2])))
  
  if(body=="lnorm"){
    par_body=as.numeric(c(out[3],exp(out[4])))
  }else if(body%in%c("gamma","weibull")){
    par_body=as.numeric(exp(out[3:4]))
  }
  
  if(tail=="lnorm"){
    par_tail=c(out[5],exp(out[6])) 
  }else if(tail%in%c("gamma","weibull")){
    par_tail=exp(out[5:6]) 
  }else if(tail=="gpd"){
    par_tail=c(0,exp(out[5]),out[6])
  }
  
  rslt=buildMixingSevdist(body,par_body,tail,par_tail,par_cauchy)
  
  return(rslt)
}
