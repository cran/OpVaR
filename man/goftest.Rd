\name{goftest}
\alias{goftest}
\title{
  Goodness of fit tests for severity distributions
}
\description{
  Testing the goodness-of-fit of the distributions given in a sevdist or a freqdist object to the loss data.
}
\usage{
  goftest(cell, object)
}
\arguments{
  \item{cell}{
    A data frame giving the losses (cell$Loss) and the user-defined period (cell$Period) in which the loss occured.
  }
    \item{object}{
    A loss severity model (sevdist object) or a loss frequency model (freqdist object)
  } 
}
\details{
  If object is of type sevdist, then Anderson-Darling, Cramer-von-Mises and Kolmogorov-Smirnov test is performed on the loss data of a single cell and the respected fitted distribution given by the sevdist object.
  
  If object is of type freqdist, then a \eqn{\chi^2} goodness-of-fit test will be performed.
}
\value{
  If object is of type sevdist, then the return is a
  list with the test results from all three tests:
  [[1]] Anderson-Darling, 
  [[2]] Cramer-von-Mises, 
  [[3]] Kolmogorow-Smirnov. 
  
  If object is of type freqdist, then the functions returns the result of the \eqn{\chi^2} test.
  
  Refer with $p.value to the p-values and $statistic to the test statistic of each of the tests.
}
\author{
 Marius Pfeuffer, Christina Zou
}

\examples{
data("lossdat")
opriskmodel = list()
opriskmodel[[1]] = list()
opriskmodel[[1]]$freqdist = fitFreqdist(lossdat[[1]],"pois")
opriskmodel[[1]]$sevdist = fitPlain(lossdat[[1]],"gamma")

# perform test on the sevdist object
goftest(lossdat[[1]], opriskmodel[[1]]$sevdist)

# show result for e.g. only the Kolmogorow-Smirnov test
test = goftest(lossdat[[1]], opriskmodel[[1]]$sevdist)
test[[3]]$p.value
test[[3]]$statistic

# perform test on the freqdist object
goftest(lossdat[[1]], opriskmodel[[1]]$freqdist)

# the p-value is given by 
goftest(lossdat[[1]], opriskmodel[[1]]$freqdist)$p.value
}

\references{
Anderson, Theodore W., and Donald A. Darling. Asymptotic theory of certain "goodness of fit" criteria based on stochastic processes. The annals of mathematical statistics (1952): 193-212.

Anderson, Theodore W., and Donald A. Darling. "A test of goodness of fit". Journal of the American statistical association 49.268 (1954): 765-769.

Birnbaum, Z. W., and Fred H. Tingey. "One-sided confidence contours for probability distribution functions". The Annals of Mathematical Statistics 22.4 (1951): 592-596.

Conover, William Jay, and William Jay Conover. "Practical nonparametric statistics" (1980).

Csorgo, Sandor, and Julian J. Faraway. "The exact and asymptotic distributions of Cramer-von Mises statistics". Journal of the Royal Statistical Society. Series B (Methodological) (1996): 221-234.

Durbin, James. "Distribution theory for tests based on the sample distribution function". Society for Industrial and Applied Mathematics, 1973.
}

\seealso{
We use the goodness-of-fit tests \code{\link{ad.test}}, \code{\link{cvm.test}} from the package \code{\link{goftest}} and \code{\link{ks.test}} and \code{\link{chisq.test}} from \code{\link{stats}}.
}