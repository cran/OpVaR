gh_inv <-
function(x,g,h) {
  sapply(x,function(x) {
    if (h>0)   {
      f1 <- function(z,x,g,h) {
        T_g(z,g)*T_h(z,h)-x
      } 
      tryCatch(stats::uniroot(f1,x=x,g=g,h=h,interval = c(-10,10),extendInt = "yes")$root,
                error=function(err) {x*Inf})
    }
  })
}
