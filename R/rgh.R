rgh <-
function(n, A,B,g,h) {
  z=stats::rnorm(n)
  A+B*T_g(z,g)*T_h(z,h)
}
