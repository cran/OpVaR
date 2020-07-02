qgh <-
function(q,A,B,g,h,log.p = FALSE) {
  sapply(q,function(q) {
    if(h >= 0) {
      u=stats::qnorm(q, log.p=log.p)
      v=A+B*T_g(u,g)*T_h(u,h)
    }
    else  {
      base::stop("negative kurtosis parameter")
    }
    v
  })
}
