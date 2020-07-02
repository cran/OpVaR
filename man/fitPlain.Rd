\name{fitPlain}
\alias{fitPlain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit plain distribution models
}
\description{
Function for estimating the parameters of plain (non-spliced, non-mixing) distribution models
}
\usage{
fitPlain(cell, family)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cell}{
list containing the data in the component cell$Loss
}
  \item{family}{
distribution family
}
}

\value{
Returns a sevdist object of type 'plain' with the given distribution fitted to the loss data.
}
\author{
Marius Pfeuffer
}

\examples{
data(lossdat)
fitPlain(lossdat[[1]],"lnorm")
}