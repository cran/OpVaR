\name{fitDependency}
\alias{fitDependency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for fitting bivariate Copulas
}
\description{
Function for fitting bivariate dependency between loss frequencies and severities.
}
\usage{
fitDependency(cell, copula_family, method = "mle")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cell}{
a lossdat cell
}
  \item{copula_family}{
0: independence, 1: Gaussian, 2: t, 3: Clayton, 4: Gumbel, 5: Frank, 6: Joe, see Vine Copula Package ?BiCop encoding
}
\item{method}{
as in fitCopula of the copula package. The method "itau" should not be used as dependencies between a continuous and a discrete random variable are modeled.
}
}
\value{
A BiCop object
}
\references{
Schepsmeier et al.: VineCopula package, CRAN
}
\author{
Marius Pfeuffer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
### Fit Joe Copula
data(lossdat)
fitDependency(lossdat[[1]],6)
}
