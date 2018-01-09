\name{gpd}
\alias{dgpd}
\alias{qgpd}
\alias{pgpd}
\alias{rgpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized Pareto Distribution
}
\description{
Density, distribution function, quantile function, and random generation for generalized Pareto distribution.
}
\usage{
dgpd(x, loc, scale, shape)
pgpd(q, loc, scale, shape) 
qgpd(p, loc, scale, shape) 
rgpd(n, loc, scale, shape) 
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
 \item{loc}{
  location parameter.
}
  \item{scale}{
 scale parameter.
}
  \item{shape}{
shape parameter.
}

}
\details{
Generalized Pareto Distribution with density

\deqn{f(x)=1/scale * (1 + shape * (x - loc)/scale)^(-1/shape - 1)}
}
\value{
dgpd gives the density, pgpd gives the distribution function, qgpd gives the quantile function, and rgpd generates random deviates.
}

\examples{
dgpd(110,100,2,4) # 0.01112233
}