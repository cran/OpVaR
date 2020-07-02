rgh <-
function(n, A,B,g,h) {
  if(h>=0) {
    z=stats::rnorm(n)
    A+B*T_g(z,g)*T_h(z,h)
  }
  else  {
    base::stop("negative kurtosis parameter")
  }
}
