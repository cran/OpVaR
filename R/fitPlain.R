fitPlain <-
function(cell,family){
  dat=cell$Loss
  mu=mean(dat)
  sigma=sd(dat)
  if(family=="weibull"){
    nll=function(x){
      -sum(dweibull(dat,exp(x[1]),exp(x[2]),log=TRUE))
    }
    start=c((sigma/mu)^(-1.086),mu/(gamma(1+1/((sigma/mu)^(-1.086)))))
    pars=exp(optim(log(start),nll)$par)
  }
  if(family=="gamma"){
    nll=function(x){
      -sum(dgamma(dat,exp(x[1]),exp(x[2]),log=TRUE))
    }
    start=c(mu^2/sigma^2,mu/sigma^2)
    pars=exp(optim(log(start),nll)$par)
  }
  if(family=="lnorm"){
    pars=c(mean(log(dat)),sd(log(dat)))
  }
  if(family=="gh"){
    q0.1 <- as.numeric(stats::quantile(dat,0.1))
    q0.9 <- as.numeric(stats::quantile(dat,0.9))
    A <- stats::median(dat)
    gamma2 <- q0.9 - q0.1
    gamma3 <- (A-q0.1)/(q0.9-A)
    gamma4 <- (as.numeric(stats::quantile(dat,0.75))-as.numeric(stats::quantile(dat,0.25)))/gamma2
    g <- -log(gamma3)/stats::qnorm(0.9)
    if (g == 0) {
      h <- 2*log(stats::qnorm(0.75) / (gamma4 * stats::qnorm(0.9))) / (stats::qnorm(0.9)^2 - stats::qnorm(0.75)^2)
    } else {
      h <- (2*log( (gamma3^(1-(stats::qnorm(0.75)/stats::qnorm(0.9))) * (gamma3^(2*stats::qnorm(0.75)/stats::qnorm(0.9)) - 1) ) / ((gamma3^2 - 1)*gamma4) )) / (stats::qnorm(0.9)^2 - stats::qnorm(0.75)^2)
    }
    B <- gamma2/(exp(0.5*h*stats::qnorm(0.9)^2)*(exp(g*stats::qnorm(0.9))-exp(-g*stats::qnorm(0.9)))/g)
    pars=c(A,B,g,h)
  }
  return(buildPlainSevdist(family,pars))
}
