\name{gh}
\alias{dgh}
\alias{qgh}
\alias{pgh}
\alias{rgh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tukey's gh Distribution
}
\description{
Density, distribution function, quantile function, and random generation for Tukey's gh distribution.
}
\usage{
dgh(x, A, B, g, h, log = FALSE)
pgh(q, A, B, g, h, log.p = FALSE)
qgh(q, A, B, g, h, log.p =FALSE)
rgh(n, A, B, g, h)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{
vector of quantiles.
}
 \item{p}{
vector of probabilities.
}
 \item{n}{
  number of observations.
}
  \item{A}{
 location parameter.
}
  \item{B}{
 scale parameter, has to be positive.
}
  \item{g}{
 skewness parameter.
}
  \item{h}{
 kurtosis parameter.
}
  \item{log, log.p}{
logical; if TRUE, probabilities p are given as log(p)
}

}
\details{
Tukey's gh distribution with location parameter A, scale parameter B, skewness parameter g, and kurtosis parameter h is obtained by transforming a standard normal variable X by

\deqn{T(X)=A+B exp(h/2 X^2)(exp(gX)-1)/g}

if g is not equal to zero, and else by 

\deqn{T(X)=A+B exp(h/2 X^2) X.}
}
\value{
dgh gives the density, pgh gives the distribution function, qgh gives the quantile function, and rgh generates random deviates.

The length of the result is determined by n for rgh, and is the length of the numerical arguments for the other function.
}
\references{
Tukey, J. W. (1960): The Practical Relationship between the Common Transformations of Counts of Amounts.
Technical Report 36, Princeton University Statistical Techniques Research Group, Princeton. 

Klein, I. and Fischer, M. (2002): Symmetrical gh-transformed Distributions. in S. Mittnek and I. Klein: \emph{Contributions to Modern Econometrics, Kluwer Academic Publishers.}

Pfaelzner, F. (2017): Einsatz von Tukey-type Verteilungen bei der Quantifizierung von operationellen Risiken. MMasterthesis Friedrich-Alexander-University Erlangen-Nueremberg.
}
\author{
Linda Moestel
}


\examples{
##Parameters  for a gh distribution
  A=500 
  B=3
  g=0.2
  h=0.5 
  
  hist(rgh(n=1000,A,B,g,h))
  curve(dgh(x,A,B,g,h),480,520)
  curve(pgh(x,A,B,g,h),480,520)
  curve(qgh(x,A,B,g,h),0,1)
}