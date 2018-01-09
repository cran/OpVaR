pgh <-
function(q,A,B,g,h,log.p=FALSE) {
  sapply(q,function(q) {
    if(g==0 && h==0) l=stats::pnorm(q,A,B)
    else if(h==0) { 
      q=(q-A)/B
      if((q<=-1/g && g>0) || (q>=-1/g && g<0))  l=NaN
      else {
        l=stats::pnorm(log(1+g*q)/g)
      }
    }
    else if(h>0) {
      q=(q-A)/B
      z=gh_inv(q,g,h)
      l=stats::pnorm(z)
    }
    else if(h<0) {
      q=(q-A)/B
      z=gh_inv(q,g,h)
      if(q<=z[4]) {
        l <- 0
      }
      else if(q<0) {
        l=min(stats::pnorm(z[1])-stats::pnorm(z[2]),0.5)
      }
      else if(q==0) {
        l=0.5
      }
      else if(q<z[3]) {
        l=1+stats::pnorm(z[1])-stats::pnorm(z[2])
      }
      else l <- 1
    }
    if(log.p==FALSE) l
    else log(l)
  })
}
