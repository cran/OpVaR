\name{fitFreqdist}
\alias{fitFreqdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fitting the frequency distribution
}
\description{
  Fitting the frequency distribution of operational losses
}
\usage{
  fitFreqdist(cell, distr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cell}{
    A data frame giving the losses (cell$Loss) and the user-defined period (cell$Period) in which the loss occured.
  }
  \item{distr}{
    The frequency distribution: either "pois" or "nbinom"
  }
  
}
\details{
  First, the number of losses per period is counted. With the results a Poisson or Negative Binomial distribution is fitted for the number of losses per period using a Maximum Likelihood estimator.
}
\value{Returns a freqdist object with the given distribution fitted to the loss frequency.}

\author{
  Christina Zou
}


\examples{
  data(lossdat)
  fitFreqdist(lossdat[[1]],"pois")
  fitFreqdist(lossdat[[2]],"nbinom")
}

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
    The frequency distribution is fitted using \code{\link{goodfit}} with method 'ML' from the package 'vcd'.
}

