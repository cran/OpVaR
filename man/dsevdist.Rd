\name{dsevdist}
\alias{dsevdist}
\alias{psevdist}
\alias{qsevdist}
\alias{rsevdist}
\alias{print.sevdist}
\alias{plot.sevdist}
\title{
Evaluating Plain, Spliced or Mixing Severity Distributions
}
\description{
Functions for evaluating pdf, cdf, quantile function and random number generation given a loss severity model (sevdist object)
}
\usage{
dsevdist(x, sevdist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{sevdist}{
%%     ~~Describe \code{sevdist} here~~
}
}

\author{
Marius Pfeuffer
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, sevdist) 
{
    if (sevdist$type == "plain") 
        return(sevdist$par[[1]](x))
    if (sevdist$type == "spliced") 
        return(dspliced(x, sevdist))
    if (sevdist$type == "mixing") 
        return(dmixing(x, sevdist))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
