\name{fitSplicedBayes}
\alias{fitSplicedBayesWeibullGPD}
\alias{fitSplicedBayesGammaGPD}
\alias{fitSplicedBayesLognormGPD}
\alias{fitSplicedBayesKDEGPD}

\title{
Parameter Estimation for Spliced Distributions
}

\description{
MCMC based parameter estimation for spliced distribution functions, where the body is modeled by a Weibull, gamma or log-normal distribution or Kernel density estimation for a log-normal distribution and the tail is fitted with a GPD.
}

\usage{
fitSplicedBayesWeibullGPD(cell, prior, burnin=10, niter=100, 
proposal_scale=evmix::fweibullgpd(cell,method="Nelder-Mead")$se, 
start=evmix::fweibullgpd(cell,method="Nelder-Mead")$optim$par)

fitSplicedBayesGammaGPD(cell, prior, burnin=10, niter=100, 
proposal_scale=evmix::fgammagpd(cell,method="Nelder-Mead")$se, 
start=evmix::fgammagpd(cell,method="Nelder-Mead")$optim$par)

fitSplicedBayesLognormGPD(cell, prior, burnin=10, niter=100, 
proposal_scale=evmix::flognormgpd(cell,method="Nelder-Mead")$se, 
start=evmix::flognormgpd(cell,method="Nelder-Mead")$optim$par)

fitSplicedBayesKDEGPD(cell, prior, burnin=10, niter=100, 
proposal_scale=evmix::flognormgpd(cell,method="Nelder-Mead")$se, 
start=evmix::flognormgpd(cell,method="Nelder-Mead")$optim$par)
}

\arguments{
\item{cell}{
  database for the parameter estimation
}
\item{prior}{
  parameters of the particular prior distribution
}
\item{burnin}{
  length of the burn-in-phase, standard value = 10
}
\item{niter}{
  MCMC sample size, standard value = 100
}
\item{proposal_scale}{
  respective scale values used for determining the prior distribution for each parameter
}
\item{start}{
  starting values (shapes) for the MH-algorithm steps 
}
}

\details{
The input contains the dataset, whose compound distribution should be fitted with the function. The second part of the input is a list with the parameters of the prior distribution. It is necessary to indicate the \code{prior} as a list containing the full parameter vector in the order \code{bodypar_1=c(),bodypar_2=c(), threshold=c(), beta=c(), xi=c()}, where \code{bodypar_1} and \code{bodypar_2} are the parameters of the body (Weibull: wshape and wscale as shape and scale parameter, log-norm: mu and sigma, gamma: gshape and gscale as shape and scale parameter).\code{threshold} is the threshold seperating the distributions.\code{beta} and \code{xi} are the scale and shape parameters of the generalized Pareto distribution.
Using the MH-algorithm the true parameters of the body and tail distribution should be approximated. 
}
\value{
Sevdist object containing the parameters of the target distribution.
}
\references{
Behrens, C. N.; Lopes, H. F.; Gamerman, D. (2004). Bayesian analysis of extreme events with threshold estimation. Statistical Modelling 4, 227-244.

Dehler, K. (2017). Bayesianische Methoden im operationellen Risiko. Master's thesis.

Ergashev, B.; Mittnik, S.; Sekeris, E. (2013). A Bayesian Approach to Extreme Value Estimation in Operational Risk Modeling. The Journal of Operational Risk 8(4), 55-81.

MacDonald, A.; Scarrott, C.J.; Lee, D.; Darlow, B.; Reale, M.; Russell, G. (2011). A flexible extreme value mixture model. Computational Statistics and Data Analysis 55(6), 2137-2157.
}
\author{
Kristina Dehler, Nicole Derfuss
}

\seealso{
Read \code{\link{buildSplicedSevdist}} for further information about building sevdist objects with type "spliced". The initial parameter values are fitted using the methods from the package 'evmix', i.e. for the WeibullGPD the function \code{\link{fweibullgpd}}. For the algorithm also the package 'truncnorm' is required.
}
\examples{
### Example for estimating the parameters of the spliced Weibull-GPD

\dontrun{

data("lossdat")
data=lossdat[[1]]$Loss

## starting values - method of moments

# Weibull distribution (shape, scale)
model <- function(x) c(F1 = x[2]*gamma(1+1/x[1])-mean(data),
	F2 = x[2]^2*(gamma(1+2/x[1])-(gamma(1+1/x[1]))^2)-sd(data)^2)
weibullpar <- rootSolve::multiroot(f = model,start = c((sd(data)/mean(data))^(-1.086),
	mean(data)/(gamma(1+1/((sd(data)/mean(data))^(-1.086))))))$root

# generalized Pareto distribution (tresh, beta, xi)
thresh=quantile(data,.9) # threshold
data2=data[which(data>thresh)]
model=function(x) c(F1=thresh+x[1]/(1-x[2])-mean(data2),
	F2=x[1]^2/(1-x[2])^2/(1-2*x[2])-sd(data2)^2)
gpdpar=c(thresh,rootSolve::multiroot(f=model,start=c(mean(data2),1/mean(data2)))$root)

## parameters of the prior distribution

prior <- list(xi = c(gpdpar[3],sd(data)),tau = c(quantile(data,.9),sd(data)/10),
	beta = c(gpdpar[2],sd(data)),wscale = c(weibullpar[2],sd(data)),
	wshape = c(weibullpar[1],sd(data)))

## estimation of (shape, scale, thresh, beta, xi)

fitSplicedBayesWeibullGPD(data, prior = prior, proposal_scale = evmix::fweibullgpd(data, 
	method = "Nelder-Mead", pvector = c(weibullpar[1], weibullpar[2], gpdpar[1], 
	gpdpar[2], gpdpar[3]))$se, start = evmix::fweibullgpd(data, 
	method = "Nelder-Mead", pvector = c(weibullpar[1], weibullpar[2], gpdpar[1], 
	gpdpar[2], gpdpar[3]))$optim$par)
}
}