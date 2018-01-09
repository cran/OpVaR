\name{buildFreqdist}
\alias{buildFreqdist}
\alias{print.freqdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Building a freqdist object
}
\description{
  Building a freqdist object
}
\usage{
  buildFreqdist(freq.distr, freq.param)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{freq.distr}{
    A string giving the name of the distribution.
    
    Distributions "pois" and "nbinom" are recognised.
  }
  \item{freq.param}{
    List of numerical values giving the parameters for the given distribution.
  }
  
}
\details{
  A freqdist object with type 'plain' is generated, i.e. for a given cell the losses occur with frequency given by distribution 'distr' with parameters given by 'param'.
}

\author{
  Christina Zou
}


\examples{
  # Poisson distributed frequency with lambda = 50
  
  freqdist1 = buildFreqdist("pois", 50)
}

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
    The sevdist objects according to different types are built via \code{\link{buildPlainSevdist}}, \code{\link{buildSplicedSevdist}} and \code{\link{buildMixingSevdist}}.
}

