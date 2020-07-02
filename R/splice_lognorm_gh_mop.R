splice_lognorm_gh_mop <-
function(data,thresh,distr,weights) {
  if(is.null(thresh)==FALSE) {
    parlist=list(NULL)
    k=0
    data2=list()
    data2[[1]]=data[data<=thresh]
    data2[[2]]=data[data>thresh]
    for(i in 2:1) {
      k=k+1
      if (distr[i]=="lnorm") {
        mu=sum(log(data2[[i]]))/length(data2[[i]])
        sigma2=sum((log(data2[[i]])-mu)^2)/length(data2[[i]])
        lnorm_par=c(mu,sigma2)
        parlist[[k]]=list(NULL)
        if(i==1) { 
          parlist[[k]][[1]]=function(x) ifelse(x<=thresh,stats::dlnorm(x,lnorm_par[1],lnorm_par[2])/stats::plnorm(thresh,lnorm_par[1],lnorm_par[2]),0)
          parlist[[k]][[2]]=function(q) ifelse(q<=thresh,stats::plnorm(q,lnorm_par[1],lnorm_par[2])/stats::plnorm(thresh,lnorm_par[1],lnorm_par[2]),1)
          parlist[[k]][[3]]=function(p) stats::qlnorm(p*stats::plnorm(thresh,lnorm_par[1],lnorm_par[2]),lnorm_par[1],lnorm_par[2])
        }
        else { 
          parlist[[k]][[1]]=function(x) ifelse(x>thresh,stats::dlnorm(x,lnorm_par[1],lnorm_par[2])/(1-stats::plnorm(thresh,lnorm_par[1],lnorm_par[2])),0)
          parlist[[k]][[2]]=function(q) ifelse(q>thresh,(stats::plnorm(q,lnorm_par[1],lnorm_par[2])-stats::plnorm(thresh,lnorm_par[1],lnorm_par[2]))/(1-stats::plnorm(thresh,lnorm_par[1],lnorm_par[2])),0)
          parlist[[k]][[3]]=function(p) stats::qlnorm(stats::plnorm(thresh,lnorm_par[1],lnorm_par[2])+p*(1-stats::plnorm(thresh,lnorm_par[1],lnorm_par[2])),lnorm_par[1],lnorm_par[2])
        }
        parlist[[k]][[4]]=lnorm_par
      }
      else if (distr[i]=="gh") {
        q0.1 <- as.numeric(stats::quantile(data2[[i]],0.1))
        q0.9 <- as.numeric(stats::quantile(data2[[i]],0.9))
        A <- stats::median(data2[[i]])
        gamma2 <- q0.9 - q0.1
        gamma3 <- (A-q0.1)/(q0.9-A)
        gamma4 <- (as.numeric(stats::quantile(data2[[i]],0.75))-as.numeric(stats::quantile(data2[[i]],0.25)))/gamma2
        g <- -log(gamma3)/stats::qnorm(0.9)
        if (g == 0) {
          h <- 2*log(stats::qnorm(0.75) / (gamma4 * stats::qnorm(0.9))) / (stats::qnorm(0.9)^2 - stats::qnorm(0.75)^2)
        } 
        else {
          h <- (2*log( (gamma3^(1-(stats::qnorm(0.75)/stats::qnorm(0.9))) * (gamma3^(2*stats::qnorm(0.75)/stats::qnorm(0.9)) - 1) ) / ((gamma3^2 - 1)*gamma4) )) / (stats::qnorm(0.9)^2 - stats::qnorm(0.75)^2)
        }
        B <- gamma2/(exp(0.5*h*stats::qnorm(0.9)^2)*(exp(g*stats::qnorm(0.9))-exp(-g*stats::qnorm(0.9)))/g)
        gh_par=c(A,B,g,h)
        parlist[[k]]=list(NULL)
        if(i==1) {  
          parlist[[k]][[1]]=function(x) ifelse(x<=thresh,dgh(x,gh_par[1],gh_par[2],gh_par[3],gh_par[4])/pgh(thresh,gh_par[1],gh_par[2],gh_par[3],gh_par[4]),0)
          parlist[[k]][[2]]=function(q) ifelse(q<=thresh,pgh(q,gh_par[1],gh_par[2],gh_par[3],gh_par[4])/pgh(thresh,gh_par[1],gh_par[2],gh_par[3],gh_par[4]),1)
          parlist[[k]][[3]]=function(p) qgh(p*pgh(thresh,gh_par[1],gh_par[2],gh_par[3],gh_par[4]),gh_par[1],gh_par[2],gh_par[3],gh_par[4])
        }
        else { ###truncated to the left
          parlist[[k]][[1]]=function(x) ifelse(x>thresh,dgh(x,gh_par[1],gh_par[2],gh_par[3],gh_par[4])/(1-pgh(thresh,gh_par[1],gh_par[2],gh_par[3],gh_par[4])),0)
          parlist[[k]][[2]]=function(q) ifelse(q>thresh,(pgh(q,gh_par[1],gh_par[2],gh_par[3],gh_par[4])-pgh(thresh,gh_par[1],gh_par[2],gh_par[3],gh_par[4]))/(1-pgh(thresh,gh_par[1],gh_par[2],gh_par[3],gh_par[4])),0)
          parlist[[k]][[3]]=function(p) qgh(pgh(thresh,gh_par[1],gh_par[2],gh_par[3],gh_par[4])+p*(1-pgh(thresh,gh_par[1],gh_par[2],gh_par[3],gh_par[4])),gh_par[1],gh_par[2],gh_par[3],gh_par[4])
        }
        parlist[[k]][[4]]=gh_par
      }
    }
    list(type="spliced",par=parlist,weights=weights,thresh=thresh)
  }
  else {
    par=list()
    if (distr=="lnorm") {
      mu=sum(log(data))/length(data)
      sigma2=sum((log(data)-mu)^2)/length(data)
      par=c(mu,sigma2)
      parlist=list(NULL)
      parlist[[1]]=function(x) stats::dlnorm(x,par[1],par[2])
      parlist[[2]]=function(q) stats::plnorm(q,par[1],par[2])
      parlist[[3]]=function(p) stats::qlnorm(p,par[1],par[2])
      parlist[[4]]=function(n) stats::rlnorm(n,par[1],par[2])
      parlist[[5]]=par
    }
    else if (distr=="gh") {
      q0.1 <- as.numeric(stats::quantile(data,0.1))
      q0.9 <- as.numeric(stats::quantile(data,0.9))
      A <- stats::median(data)
      gamma2 <- q0.9 - q0.1
      gamma3 <- (A-q0.1)/(q0.9-A)
      gamma4 <- (as.numeric(stats::quantile(data,0.75))-as.numeric(stats::quantile(data,0.25)))/gamma2
      g <- -log(gamma3)/stats::qnorm(0.9)
      if (g == 0) {
        h <- 2*log(stats::qnorm(0.75) / (gamma4 * stats::qnorm(0.9))) / (stats::qnorm(0.9)^2 - stats::qnorm(0.75)^2)
      } 
      else {
        h <- (2*log( (gamma3^(1-(stats::qnorm(0.75)/stats::qnorm(0.9))) * (gamma3^(2*stats::qnorm(0.75)/stats::qnorm(0.9)) - 1) ) / ((gamma3^2 - 1)*gamma4) )) / (stats::qnorm(0.9)^2 - stats::qnorm(0.75)^2)
      }
      B <- gamma2/(exp(0.5*h*stats::qnorm(0.9)^2)*(exp(g*stats::qnorm(0.9))-exp(-g*stats::qnorm(0.9)))/g)
      par=c(A,B,g,h)
      parlist=list(NULL)
      parlist[[1]]=function(x) dgh(x,par[1],par[2],par[3],par[4])
      parlist[[2]]=function(q) pgh(q,par[1],par[2],par[3],par[4])
      parlist[[3]]=function(p) qgh(p,par[1],par[2],par[3],par[4])
      parlist[[4]]=function(n) rgh(n,par[1],par[2],par[3],par[4])
      parlist[[5]]=par
    }
    list(type="plain",par=parlist)
  }
}
