print.freqdist<-function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\nFrequency Distribution Object:\n---\n")
    cat("Distribution: ")
    print(x[[1]])
    cat("Parameters: ")
    print(x[[2]])
  }
