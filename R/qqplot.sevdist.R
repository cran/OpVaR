qqplot.sevdist <- function(cell, sevdist, probabilities = NULL){
  
  if (is.null(probabilities)){
    probabilities = seq(0.01, 0.99, 0.001)
  }
  if (sevdist$type == "plain"){
    quantiles = qsevdist(x = probabilities, sevdist = sevdist)
    stats::qqplot(cell$Loss, quantiles, xlab = "Experimental Loss Quantiles", 
           ylab = "Theoretical Loss Quantiles", main = 
             paste("Plain:", paste(toupper(substr(sevdist$par[[6]][[1]][[1]], 1, 1)),
                                  substr(sevdist$par[[6]][[1]][[1]],
                                         2, nchar(sevdist$par[[6]][[1]][[1]])),
                                  sep = "") ))
    abline(a = 0, b = 1, col = "royalblue", lwd = 2) 
  }
  
  
  if (sevdist$type == "spliced") {
    quantiles = qsevdist(x = probabilities, sevdist = sevdist)
    stats::qqplot(cell$Loss, quantiles, xlab = "Experimental Loss Quantiles", 
           ylab = "Theoretical Loss Quantiles", main = 
             paste("Spliced:", paste(toupper(substr(sevdist$par[[1]][[6]][[1]], 1, 1)),
                                    substr(sevdist$par[[1]][[6]][[1]], 2,
                                           nchar(sevdist$par[[1]][[6]][[1]])),sep = ""),
                   "(Body)", paste(toupper(substr(sevdist$par[[2]][[6]][[1]], 1, 1)),
                                   substr(sevdist$par[[2]][[6]][[1]], 2,
                                          nchar(sevdist$par[[2]][[6]][[1]])),
                                   sep = ""), "(Tail)"))
    abline(a=0, b=1, col = "royalblue", lwd = 2) 
  }
  
  if (sevdist$type == "mixing") {
    quantiles = qmixing(p = probabilities, sevdist = sevdist)
    stats::qqplot(cell$Loss, quantiles, xlab = "Experimental Loss Quantiles",
           ylab = "Theoretical Loss Quantiles", main = 
             paste("Mixing:", paste(toupper(substr(sevdist$par[[2]][[6]][[1]], 1, 1)),
                                   substr(sevdist$par[[2]][[6]][[1]], 2, 
                                          nchar(sevdist$par[[2]][[6]][[1]])), sep = ""),
                   "(Body)", paste(toupper(substr(sevdist$par[[1]][[6]][[1]], 1, 1)),
                                   substr(sevdist$par[[1]][[6]][[1]], 2, 
                                          nchar(sevdist$par[[1]][[6]][[1]])),
                                   sep = ""), "(Tail)" ))
    abline(a=0, b=1, col = "royalblue", lwd = 2)
    
  }
}