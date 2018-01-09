dgh <-
function(x,A,B,g,h,log=FALSE) {
  sapply(x,function(x) {
    x=(x-A)/B
    if(g==0 && h==0) v=stats::dnorm(x)
    else if(h==0) { 
      if((x<=-1/g && g>0) || (x>=-1/g && g<0))  y=NaN
      else y=log(1+g*x)/g
      v=stats::dnorm(y)/(B*exp(g*y))
    }
    else if(h > 0) {
      y=gh_inv(x,g,h)
      v=stats::dnorm(y)/(B*deriv_gh(y,g,h))
    }
    else if(h < 0) {
      y=gh_inv(x,g,h)
      if(x > y[4] && x < y[3]) {
        if(y[2]==Inf) v1=0
        else v1=stats::dnorm(y[2])/(B*deriv_gh(y[2],g,h))
        v=stats::dnorm(y[1])/(B*deriv_gh(y[1],g,h))+v1
      }
      else v=NaN
    }
    if(log==FALSE) v
    else log(v)
  })
}
