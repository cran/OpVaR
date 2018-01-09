deriv_gh <-
function(z,g,h) {
  if(g==0) exp(h*z^2*0.5)*(1+h*z^2)
  else exp(h*z^2*0.5)*(exp(g*z)+h*z*(exp(g*z)-1)/g)
}
