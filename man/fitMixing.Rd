\name{fitMixing}
\alias{fitMixing}

\title{
Maximum Likelihood Estimation 
}
\description{
Maximum Likelihood Estimation for the dynamic weighted mixture model of Frigessi et al., 2002
}
\usage{
fitMixing(cell, body, tail, method="L-BFGS-B", c_location0=0.75, c_scale0=2)
}

\arguments{
 \item{cell}{
  lossdat cell
  }
  \item{body}{
  body distribution, either gamma, lnorm or weibull
  }
  \item{tail}{
  tail distribution, either gamma, lnorm, weibull or gpd
  }
  \item{method}{
  optimization method, default is "L-BFGS-B"
  }
  \item{c_location0}{
  empirical quantile of loss severity data used for initialization of Cauchy location parameter in optimization: quantile(cell$Loss,c_location0)
  }
  \item{c_scale0}{
  scaling factor for empirical standard deviation used for initialization of Cauchy scale parameter in optimization: sd(cell$Loss)/c_scale0 
  }
}
\details{
Body and tail parameters are initialized by method of moments estimators. Cauchy location is initialized by empirical 70% quantile, Cauchy scale by 1/10 of the empirical standard deviation.
}
\value{
Returns a sevdist object of type 'mixing' with the given body and tail distributions fitted to the loss data.
}
\author{
Marius Pfeuffer
}
\references{
Frigessi et al. \emph{A Dynamic Mixture Model for Unsupervised Tail Estimation without Threshold Selection}, Extremes 5(3):219-235, 2003
}

\examples{
\dontrun{
data(lossdat)
sev=fitMixing(lossdat[[1]],"weibull","gpd")
sev
plot(sev,5000)
}
}

