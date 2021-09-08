\name{qqplot.sevdist}
\alias{plot.sevdist}
\alias{qqplot.sevdist}


\title{
Density plot and Q-Q plot for plain, mixing, and spliced distributions
}
\description{
Given a dataset and a fitted sevdist object of type plain, mixing, or spliced, a density plot or qqplot is created to help determine the best distribution for modelling loss severity.
}

\usage{
qqplot.sevdist(cell, sevdist, probabilities = NULL)
}

\arguments{
  \item{sevdist}{
  Fitted sevdist object
}
  \item{cell}{
  List containing the data in the component cell$Loss
  }
  \item{probabilities}{
  Vector containing probabilities for desired quantiles. If no vector
  is given,
  default probabilities of 0.001, 0.002, ...., 0.99 will be used.
  }
}
\value{
Returns a density plot of the fitted sevdist object or a qqplot of the theoretical loss quantiles (sevdist object) vs. the experimental loss quantiles (loss data).
}

\author{
Nina Buoni, Marius Pfeuffer, Christina Zou
}

\examples{
data(lossdat)
sev1 <- fitPlain(lossdat[[1]], "weibull")
plot(sev1)

sev2 <- fitSpliced(lossdat[[2]], "gamma", "gpd", method = "dAMSE")
qqplot.sevdist(lossdat[[2]], sev2)

}
