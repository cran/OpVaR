print.sevdist<-function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\nSeverity Distribution Object:\n---\n")
  cat("Distribution Type: ")
  print(x$type)
  if(x$type=="plain") {
    cat("---\nDistribution Family: ")
    print(x$par[[6]])
    cat("Parameters: ")
    print(x$par[[5]])
  }
  if(x$type=="mixing") {
    cat("---\nBody Distribution Family: ")
    print(x$par[[2]][[6]])
    cat("Parameters: ")
    print(x$par[[2]][[5]])
    cat("Tail Distribution Family: ")
    print(x$par[[1]][[6]])
    cat("Parameters: ")
    print(x$par[[1]][[5]])
    cat("---\nMixing CDF: ")
    print("cauchy")
    cat("Parameters: ")
    print(x$par[[3]][[5]])
    
  }
  if(x$type=="spliced") {
    cat("---\nBody Distribution Family: ")
    print(x$par[[2]][[6]])
    cat("Parameters: ")
    print(x$par[[2]][[5]])
    cat("Tail Distribution Family: ")
    print(x$par[[1]][[6]])
    cat("Parameters: ")
    print(x$par[[1]][[5]])
    cat("---\nThreshold: ")
    print(x$thresh)
  }
}

plot.sevdist=function (x, xmax = NULL, npoints = 500, main = "Severity Distribution", ...)
{
  sevdist = x
  if (x$type == "plain") {
    if (is.null(xmax))
      xmax = qsevdist(0.99, x)
    xseq = seq(0, xmax, length.out = 1000)
    plot(xseq, dsevdist(xseq, sevdist), t = "l", col = "royalblue",
         lwd = 2, main = main, xlab = "Loss", ylab = "Density")
  }
  if (x$type == "spliced") {
    if (is.null(xmax))
      xmax = qsevdist(0.99, x)
    thresh = sevdist$thresh
    xseq1 = seq(0, floor(thresh - 0.1), length.out = ceiling(npoints *
                                                               thresh/xmax))
    xseq2 = seq(ceiling(thresh + 0.1), xmax, length.out = ceiling(npoints *
                                                                    (xmax - thresh)/xmax))
    f = dsevdist(c(xseq1, xseq2), sevdist)
    plot(xseq1, f[1:length(xseq1)], t = "l", col = "royalblue",
         lwd = 2, main = main, xlab = "Loss",
         ylab = "Density", xlim = c(0, xmax))
    lines(xseq2, f[(length(xseq1) + 1):(length(xseq1) + length(xseq2))],
          col = "royalblue", lwd = 2)
    par(xpd = FALSE)
    abline(v = sevdist$thresh, lty = 2, col = "darkorange4")
  }
  if (sevdist$type == "mixing") {
    if (is.null(xmax))
      xmax = max(rsevdist(200, sevdist))
    curve(dsevdist(x, sevdist), 0, xmax, col = "royalblue",
          lwd = 2, main = main, xlab = "Loss",
          ylab = "Density")
    abline(v = sevdist$par[[3]][[5]][1], lty = 2, col = "darkorange4")
  }
}