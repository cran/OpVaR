\name{buildSplicedSevdist}
\alias{buildSplicedSevdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Building a sevdist object with a spliced distribution
}
\description{
  Building a sevdist object with a spliced distribution
}
\usage{
  buildSplicedSevdist(body.distr, body.param, tail.distr, tail.param, thresh, weight)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{body.distr}{
    A string giving the name of the body distribution.
  }
  \item{body.param}{
    Vector of the parameters for the given body distribution.
  }
  \item{tail.distr}{
    A string giving the name of the tail distribution.
    
    Distributions "lgamma", "weibull", "gpd", "gh" are recognised.
  }
  \item{tail.param}{
    Vector of the parameters for the given tail distribution.
  }
  \item{thresh}{
    Numeric value giving the threshold of the distribution, the severity follows a truncated 'body.distr' distribution below the threshold and a truncated 'tail.distr' distribution above the threshold.
  }
  \item{weight}{
    Numeric value in [0,1] giving the probability that the severity distribution takes values below the threshold.
  }
  
}
\details{
  A sevdist object with type 'spliced' is generated, i.e. the severity has distribution 'body.distr' in the body with parameters given by 'body.param' and distribution 'tail.distr' in the tail with parameters given by 'tail.param'.
  
  The density \eqn{f(x)} for a spliced distribution with threshold \eqn{\tau}, weights \eqn{w}, body distribution with density \eqn{g} and cumulative distribution function \eqn{G} and tail distribution with density \eqn{h} and cumulative distribution function \eqn{H} is given as below:
    
    If \eqn{x\le \tau}: \eqn{f(x) = w* g(x)/G(\tau)}
  
  If \eqn{x>\tau}: \eqn{f(x) = (1-w)* h(x)/(1-H(\tau))}
}

\author{
  Christina Zou
}


\examples{
  # Spliced distributed severity with gamma distributed body with shape = 1.23, rate = 0.12 
  # and GPD distributed tail with shape = 716 and scale = 0.1 and threshold = 2000. 
  # The weight for the body is 0.8.
  
  sevdist1 = buildSplicedSevdist("gamma", c(2.5, 0.012), "gpd", c(2000, 716, 0.1), 2000, 0.8)
  plot(sevdist1)
}

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
    Other sevdist objects with type 'plain' and 'mixing' are built via \code{\link{buildPlainSevdist}} and \code{\link{buildMixingSevdist}}.
}

