\name{buildMixingSevdist}
\alias{buildMixingSevdist}\alias{evaluate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Building a dynamic mixture model as a sevdist object
}
\description{
  Building a dynamic mixture model as a sevdist object
}
\usage{
  buildMixingSevdist(body.distr, body.param, tail.distr, tail.param, mixing.param)
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
  }
  \item{tail.param}{
    Vector of the parameters for the given tail distribution.
    
    Distributions "lgamma", "weibull", "gpd", "gh" are recognised.
  }
  \item{mixing.param}{
    Vector of the parameters for the given mixing distribution.
  }
  
}
\details{
  A sevdist object with type 'mixing' is generated, i.e. the distribution of the severity is a dynamic mixture of the distributions 'body.distr' with parameters given by 'body.param' and 'tail.distr' with parameters given by 'tail.param' using the cdf of a mixing distribution 'mixing.distr' with parameters 'mixing.param' according to Frigessi et al. (2002).
  
  The density \eqn{f(x)} for a dynamic mixture model with body distribution having density \eqn{g} and tail distribution having density \eqn{h} and cumulative distribution function \eqn{p} of the mixing distribution is given below:
    
    \deqn{f(x) = C*[p(x)*g(x) +(1-p(x))*h(x)]}
}
\author{
  Christina Zou
}

\examples{
  # Create mixing sevdist object
  sevdist = buildMixingSevdist("weibull", c(2,6), "gpd", c(0,4,-0.5), c(8,.8))
  # Plot pdf
  plot(sevdist)
}



\references{
  Frigessi, A., Haug, O., & Rue, H. (2002). A dynamic mixture model for unsupervised tail estimation without threshold selection. Extremes, 5(3), 219-235.
}

