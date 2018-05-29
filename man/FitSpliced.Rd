\name{fitSpliced}
\alias{fitSpliced}
\alias{fitSplicedPar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimation of the threshold, the body and the tail parameters for a spliced distribution
}

\description{
  Given a dataset with a chosen distribution for the data in the body and another distribution in the tail, the threshold, the parameters of the body and the tail distributions and the weights are estimated.
}

\usage{
  fitSpliced(cell, body, tail, method, thresh = NULL)

  fitSplicedPar(cell, thresh, body, tail)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cell}{
   List containing the data in the component cell$Loss
  }
  \item{body}{
   Distribution in the body. Can be chosen between "gamma", "lnorm", "weibull" or "erlang"
  }
  \item{tail}{
  Distribution in the tail. Can be chosen between "gpd", "gamma", "lnorm", "weibull" or "gh"
  }
  \item{method}{
  Method for the threshold estimation. In case of a GPD tail, it can be chosen between "BestFit", "Fixed", "dAMSE", "danielsson", "DK", "hall", "Himp", "HW", "mindist" or "RT"
}
\item{thresh}{
 Predetermined threshold quantile (if method="Fixed")
}
}
\value{Returns a sevdist object of type 'spliced' with the given body and tail distributions fitted to the loss data.}
\details{
  In the spliced model, the distribution of the data is spliced into two distributions, one in the body and another in the tail.
  Given the threshold t, a density function \eqn{d1} for the body distribution and a density function \eqn{d2} for the tail distribution, the incidental spliced density function
  has the formula \eqn{d(x)=w*d1(x)} if \eqn{x<=t} and \eqn{d(x)=(1-w)*d2(x)} of \eqn{x>t}.
  The weight w is needed to normalize the density function.
  The estimation of the spliced distribution consists of three steps.
  In the first step, the threshold given chosen body and tail distribution is estimated.
  For further details on the methods for the threshold estimation, see fitThreshold.
  In the second step, the parameters of the body and the tail distribution are obtained by maximum likelihood estimation. 
  For spliced distributions, truncated body distributions are fitted to truncated data. For the estimation of the parameters of the tail distribution, only the data points above the threshold are used.
  In the last step, the weight w is obtained by maximum likelihood estimation.
  
}


\author{
  Christina Zou, Leonie Wicht
}

\examples{
 data(lossdat)
 fitSpliced(lossdat[[3]],"gamma","gpd",method="Fixed",thresh=2000)
}

  \seealso{
    fitSplicedBestFit, fitThreshold
  }
