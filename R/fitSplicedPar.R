#library(REIns)


fitSplicedPar=function(cell, thresh, body, tail)
  
{
  
dat=cell$Loss
mu=mean(dat)
sigma_sq=var(dat)
sd=sd(dat)
n=length(dat)
body_dat=dat[which(dat<=thresh)]
tail_dat=dat[which(dat>thresh)]
  
if (! body %in% c("lnorm", "weibull", "gamma", "erlang")) stop("Body must be lnorm, weibull, gamma or erlang")
if (! tail %in% c("gpd", "lnorm", "weibull", "gamma", "gh")) stop("Tail must be gpd, lnorm, weibull, gamma or gh")
  
  

if (body=="gamma") {
  nll=function(x){
  y=exp(x)
  -sum(dgamma(dat[which(dat<thresh)],y[1],y[2],log=TRUE)-pgamma(thresh,y[1],y[2],log.p=TRUE))
  }
  par0=c(mu^2/sigma_sq,mu/sigma_sq)
  body_param=exp(optim(log(par0),nll)$par)
}
  

if (body=="lnorm") {
  nll=function(x) {
  y=c(x[1], exp(x[2]))
  -sum(dlnorm(dat[which(dat<thresh)],y[1],y[2],log=TRUE)-plnorm(thresh,y[1],y[2],log.p=TRUE))
  }
 par0=c(log(mu), log(sqrt(sigma_sq)))
 out=optim(c(par0[1], log(par0[2])), nll)$par
 body_param=c(out[1],exp(out[2]))
  }

if (body=="weibull") { 
  nll_weibull_body <- function(x){
  y <- exp(x)
  -sum(dweibull(body_dat, y[1], y[2], log=TRUE) - pweibull(thresh, y[1], y[2], log.p = TRUE))
  }
  med=median(cell$Loss)
  qu=quantile(cell$Loss, c(0.25, 0.75))
  k=log(log(0.25)/log(0.75))/log(qu[2]/qu[1])
  c=med/(log(2)^(1/k))
  par0=as.numeric(c(k, c))
  
  body_param <- as.numeric(exp(optim(log(par0), nll_weibull_body)$par))
}
    
if (body=="erlang") {
      fit=SpliceFitGPD(dat, const=NULL, tsplice=thresh)
      shape_ME=fit$MEfit$shape
      scale_ME=fit$MEfit$theta
      body_param=c(shape_ME, scale_ME)
      
      loc_GPD=thresh
      scale_GPD=fit$EVTfit$sigma
      shape_GPD=1/(fit$EVTfit$gamma)
      tail_param=c(loc_GPD, scale_GPD, shape_GPD)

}
      

if (tail=="gpd" & ! body=="erlang") {
  nll=function(x){
    -sum(dgpd(dat[which(dat > thresh)],thresh,exp(x[1]),x[2],log=T))
  }
  xi_mom=.5*(1-((mean(dat[which(dat > thresh)])-thresh)^2)/var(dat[which(dat > thresh)]))
  if(xi_mom<=0){
    xi_mom=0
  }
  if(xi_mom>=.5){
    xi_mom=.5
  }
  var_mom=(mean(dat[which(dat > thresh)]))*(1-xi_mom)
  fit=optim(c(log(var_mom),xi_mom),nll)
  tail_param=c(thresh,exp(fit$par[1]),fit$par[2])
  }
      
      
if (tail=="lnorm") {
  nll=function(x) {
  y=c(x[1], exp(x[2]))
  -sum(dlnorm(dat[which(dat>thresh)],y[1],y[2],log=TRUE)-log((1-plnorm(thresh,y[1],y[2]))))
  }
  par0=c(log(mu), log(sqrt(sigma_sq)))
  out=exp(optim(c(par0[1], log(par0[2])), nll)$par)
  tail_param=c(out[1],exp(out[2]))
  }
      
      
if (tail=="weibull") {
  nll_weibull_tail <- function(x){
  y <- exp(x)
  - sum(dweibull(tail_dat, y[1], y[2], log=TRUE) - log(1-pweibull(thresh, y[1], y[2])))
  }
  med=median(cell$Loss)
  qu=quantile(cell$Loss, c(0.25, 0.75))
  k=log(log(0.25)/log(0.75))/log(qu[2]/qu[1])
  c=med/(log(2)^(1/k))
  par0 <- c((sd/mu)^(-1.086), mu/(gamma(1+1/((sd/mu)^(-1.086)))))
  #nll_weibull_tail(par0)
  #par0=(c(k, c))
  #nll_weibull_tail(par0)
  tail_param=exp(optim(log(as.numeric(par0)), nll_weibull_tail)$par)
  }
      
if (tail=="gamma") {
  nll_gamma_tail<-function(x) {
      y<-exp(x)
      -sum(dgamma(tail_dat, y[1], y[2], log=TRUE)-log(1-pgamma(thresh, y[1], y[2])))
      }
  par0=c(mu^2/sigma_sq,mu/sigma_sq)
  tail_param=exp(optim(log(par0), nll_gamma_tail)$par)
}
        
if (tail == "gh") {
  q0.1 <- as.numeric(stats::quantile(dat, 0.1))
  q0.9 <- as.numeric(stats::quantile(dat, 0.9))
  A <- stats::median(dat)
  gamma2 <- q0.9 - q0.1
  gamma3 <- (A - q0.1)/(q0.9 - A)
  gamma4 <- (as.numeric(stats::quantile(dat, 0.75)) - as.numeric(stats::quantile(dat, 0.25)))/gamma2
  g <- -log(gamma3)/stats::qnorm(0.9)
  if (g == 0) {
    h <- 2 * log(stats::qnorm(0.75)/(gamma4 * stats::qnorm(0.9)))/(stats::qnorm(0.9)^2 - stats::qnorm(0.75)^2)
  }
  else {
    h <- (2 * log((gamma3^(1 - (stats::qnorm(0.75)/stats::qnorm(0.9))) * (gamma3^(2 * stats::qnorm(0.75)/stats::qnorm(0.9)) -
                                                                            1))/((gamma3^2 - 1) * gamma4)))/(stats::qnorm(0.9)^2 - stats::qnorm(0.75)^2)
  }
  B <- gamma2/(exp(0.5 * h * stats::qnorm(0.9)^2) * (exp(g * stats::qnorm(0.9)) - exp(-g * stats::qnorm(0.9)))/g)
  if(h<0) base::warning("estimated kurtosis parameter is negative")
  pars = c(A, B, g, h)
}
            
return(list(body_param, tail_param))
}
            
            
