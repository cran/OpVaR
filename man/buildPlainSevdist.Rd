\name{buildPlainSevdist}
\alias{buildPlainSevdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Building a sevdist object with a plain distribution
}
\description{
  Building a sevdist object with a plain distribution
}
\usage{
  buildPlainSevdist(distr, param)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distr}{
    A string giving the name of the distribution.
  }
  \item{param}{
    The parameters for the given distribution
  }
  
}
\details{
  A sevdist object with type 'plain' is generated, i.e. the severity has distribution 'distr' with parameters given by 'param'.
}

\author{
  Christina Zou
}


\examples{
  # Log-gamma distributed severity with shape = 2.2 and rate = 1.2
  
  sevdist1 = buildPlainSevdist("lgamma", c(2.2, 1.2))
  plot(sevdist1)
}

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
    Other sevdist objects with type 'spliced' and 'mixing' are built via \code{\link{buildSplicedSevdist}} and \code{\link{buildMixingSevdist}}.
}

