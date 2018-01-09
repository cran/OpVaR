\name{Mixing}
\alias{dmixing}
\alias{pmixing}
\alias{qmixing}
\alias{rmixing}

\title{
  Extended Dynamic Weighted Mixture Model 
}
\description{
  Implementation of the dynamic weighted mixture model of Frigessi et al., 2001
}
\usage{
  dmixing(x,sevdist,log=FALSE)
  pmixing(q,sevdist,lower.tail=TRUE)
  qmixing(p,sevdist,lower.tail=TRUE)
  rmixing(n,sevdist)
}

\arguments{
  \item{x,q}{
    vector of quantiles
  }
  \item{p}{
    vector of probabilities.
  }
  \item{n}{
    positive integer, number of observations
  }
  \item{sevdist}{
    sevdist object
  }
  \item{log}{
    logical; if TRUE, probabilities p are given as log(p).
  }
  \item{lower.tail}{
    logical; if TRUE (default), probabilities are P[X <= x] otherwise, P[X > x].
  }
}
\details{
  The dynamic mixture pdf is given by l(x)=1/r*((1-W(x))*f_1(x)+W(x)*f_2(x)), where W(x) denotes a Cauchy cdf, f_1(x) and f_2(x) body and tail pdfs, r is a normalizing constant.
}
\author{
Marius Pfeuffer
}

\references{
  A. Frigessi et al.: \emph{A Dynamic Mixture Model for Unsupervised Tail Estimation without Threshold Selection}, Extremes 5(3):219-235, 2002
  Y. Hu: \emph{User's Guide for evmix Package in R.}, Working Paper, 2013
}

\examples{
  ### Create mixing sevdist object
  sevdist=buildMixingSevdist("weibull", c(2,6), "gpd", c(0,4,-0.5), c(8,.8))
  ### Evaluate pdf
  dmixing(1,sevdist)
}
