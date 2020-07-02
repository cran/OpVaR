\name{fitSpliced}
\alias{fitSpliced}
\alias{fitSplicedPar}
\alias{dspliced}
\alias{qspliced}
\alias{pspliced}
\alias{rspliced}
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
  The density \eqn{f(x)} for a spliced distribution with threshold \eqn{\tau}, weights \eqn{w}, body distribution with density \eqn{g} and cumulative distribution function \eqn{G} and tail distribution with density \eqn{h} and cumulative distribution function \eqn{H} is given as below:
    
    If \eqn{x\le \tau}: \eqn{f(x) = w* g(x)/G(\tau)}
  
  If \eqn{x>\tau}: \eqn{f(x) = (1-w)* h(x)/(1-H(\tau))}

  The weight \eqn{w} is required to normalise the density function.
  The estimation of the spliced distribution consists of three steps.
  In the first step, the threshold for the selected body and tail distribution is estimated.
  For further details on the methods for the threshold estimation, see fitThreshold.
  In the second step, the parameters of the body and the tail distribution are obtained by maximum likelihood estimation. 
  For spliced distributions, truncated body distributions are fitted to truncated data. For the estimation of the parameters of the tail distribution, only the data points above the threshold are used.
  In the last step, the weight \eqn{w} is obtained.
  
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
