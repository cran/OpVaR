\name{fitSplicedBestFit}
\alias{fitSplicedBestFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fitting a spliced distribution over a given data set 
}

\description{
  A spliced distribution is fitted based on given loss data and required distributions for the body and the tail using a best-fit method. The function returns a sevdist object with type "spliced".
}

\usage{
  fitSplicedBestFit(cell, body, tail, thresh0 = 0.7, thresh.max = 0.98)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cell}{
    List containing the loss data in cell$Loss.}
  \item{body}{
    Character string giving the name of the distribution in the body. Can be chosen between "gamma", "lnorm", "weibull" or "erlang".
  }
   \item{tail}{
  Character string giving the name of the distribution in the tail. Can be chosen between "gpd", "gamma", "lnorm", "weibull" or "gh".
  }
  \item{thresh0}{Initial value for the quantile of the threshold between body and tail}
  \item{thresh.max}{Terminal value for the quantile of the threshold between body and tail}

}

\details{
  Starting with a spliced distribution with threshold given by the quantile thresh0, the respective parameters of the given body and tail distribution are fitted, as well as the weights of both parts. A list of possible thresholds is given by all loss values between all quantiles between thresh0 and thresh.max. For each of these thresholds a spliced distribution will be fitted and then the best fitting threshold with the corresponding parameters will be chosen with a Kolmogorow-Smirnov test and the result is returned in a sevdist object with type "spliced". 
  }
  
  \author{
    Christina Zou
  }
  
  
  \examples{
    data(lossdat)
    \donttest{
    param<-fitSplicedBestFit(lossdat[[1]], "lnorm", "gpd")
    param[[1]] ## Parameters of the body distribution
    param[[2]] ## Parameters of the tail distribution
    }
  }
  
  \seealso{
    fitSplicedPar, fitThreshold, fitSpliced
  }
  
  % Add one or more standard keywords, see file 'KEYWORDS' in the
  % R documentation directory.