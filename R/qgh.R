qgh <-
function(q,A,B,g,h,log.p=FALSE) {
  sapply(q,function(q) {
    if(h >= 0) {
      u=stats::qnorm(q)
      v=A+B*T_g(u,g)*T_h(u,h)
    }
    else {
      z=gh_inv(A,g,h)
      toroot <- function(x,q,A,B,g,h) {
        pgh(x,A,B,g,h)-q
      }
      if(q==0) v=A+B*z[4]
      else if(q==0.5) v=A
      else if(q!=1) v=stats::uniroot(toroot,q=q,A=A,B=B,g=g,h=h,interval=c(z[5],z[6]),extendInt = "yes")$root
      else v=A+B*z[3]
    }
    if(log.p==FALSE) v
    else{
      if (v!=-Inf) log(v)
      else NaN
    }
  })
}
