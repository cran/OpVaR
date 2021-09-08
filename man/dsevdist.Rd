\name{dsevdist}
\alias{dsevdist}
\alias{psevdist}
\alias{qsevdist}
\alias{rsevdist}
\alias{print.sevdist}
\title{
Evaluating Plain, Spliced or Mixing Severity Distributions
}
\description{
Functions for evaluating pdf, cdf, quantile function and random number generation given a loss severity model (sevdist object)
}
\usage{
dsevdist(x, sevdist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{sevdist}{
%%     ~~Describe \code{sevdist} here~~
}
}

\author{
Marius Pfeuffer
}

\examples{

## The function is currently defined as
function (x, sevdist) 
{
    if (sevdist$type == "plain") 
        return(sevdist$par[[1]](x))
    if (sevdist$type == "spliced") 
        return(dspliced(x, sevdist))
    if (sevdist$type == "mixing") 
        return(dmixing(x, sevdist))
  }
}

