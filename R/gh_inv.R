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
    else {
      if(g!=0) {
        f2 <- function(z,g,h) {
          v=(g+h*z)*exp(g*z)-h*z
          v=max(v,-.Machine$double.xmax)
          v=min(v,.Machine$double.xmax)
          v
        }
        limit1=stats::uniroot(f2,g=g,h=h,interval=c(-10,0),extendInt="yes")$root
        limit2=stats::uniroot(f2,g=g,h=h,interval=c(0,10),extendInt="yes")$root
      }
      else {
        limit1=-sqrt(-1/h)
        limit2=sqrt(-1/h)
      }
      GH_toroot <- function(z,x) {
        (T_g(z,g)*T_h(z,h)) - x
      }
      max=(T_g(limit2,g)*T_h(limit2,h))
      min=(T_g(limit1,g)*T_h(limit1,h))
      if(is.na(max)==TRUE) max=0
      GH_upper_toroot <- function(z,x) {
        if(z>limit2)  (T_g(z,g)*T_h(z,h)) - x
        else  max - x
      }
      GH_lower_toroot <- function(z,x) {
        if(z<limit1)  (T_g(z,g)*T_h(z,h)) - x
        else  min - x
      }
      if(x < min) z1 <- z2 <- 0
      else if(x<0) {
        z1=tryCatch(stats::uniroot(GH_toroot,x=x,interval=c(limit1,limit2))$root,error=function(err) {-Inf})
        z2=tryCatch(stats::uniroot(GH_lower_toroot,x=x,interval=c(-10,limit1),extendInt = "downX")$root,error=function(err) {-Inf})
      }
      else if(x==0) {
        z1=tryCatch(stats::uniroot(GH_toroot,x=x,interval=c(limit1,limit2))$root,error=function(err) {-Inf})
        z2=tryCatch(stats::uniroot(GH_upper_toroot,x=x,interval=c(limit2,100),extendInt = "upX")$root,error=function(err) {Inf})
      }
      else if(x<max) {
        z1=tryCatch(stats::uniroot(GH_toroot,x=x,interval=c(limit1,limit2))$root,error=function(err) {Inf})
        z2=tryCatch(stats::uniroot(GH_upper_toroot,x=x,interval=c(limit2,100),extendInt = "upX")$root,error=function(err) {Inf})
      }
      else z1 <- z2 <- 1
      c(z1,z2,max,min,limit1,limit2)
    }
  })
}
