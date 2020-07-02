\name{mcSim}
\alias{mcSim}
\alias{VaR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo Simulation from opriskmodel objects for total loss estimation
}
\description{
Function for conducting Monte Carlo Simulation of complete opriskmodel objects (list of cells with (1) frequency model, (2) severity model and (3) dependencymodel)
}
\usage{
mcSim(opriskmodel, n_sim, verbose=TRUE)
VaR(mc_out, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{opriskmodel}{
an opriskmodel object
}
  \item{n_sim}{
number of simulations
}
 \item{mc_out}{
Monte Carlo simulation output
}
  \item{alpha}{
significance level for quantile/value-at-risk
}
 \item{verbose}{
verbose mode
}
}
\value{
A mcsim object, which can be further processed by the VaR function to estimate empirical quantiles as value-at-risk measure
}
\author{
Marius Pfeuffer
}

\seealso{
\code{\link{sla}}
}
\examples{
\donttest{
### Load Example Data Set
data(lossdat)

### Estimation of Complete Risk Model
opriskmodel1=list()
for(i in 1:length(lossdat)){
  opriskmodel1[[i]]=list()
  opriskmodel1[[i]]$freqdist=fitFreqdist(lossdat[[i]],"pois")
  opriskmodel1[[i]]$sevdist=fitPlain(lossdat[[i]],"lnorm")
}

### Cell 1: Gumbel Copula, Cell 2: Independence, Cell 3: Frank Copula, Cell 4: Independence
opriskmodel1[[1]]$dependency=fitDependency(lossdat[[1]],6)
opriskmodel1[[3]]$dependency=fitDependency(lossdat[[3]],4)

### Monte Carlo Simulation
mc_out=mcSim(opriskmodel1,100)

### Evaluation of 95
VaR(mc_out,.95)
sla(opriskmodel1,.95)

### Monte Carlo Simulation
mc_out=mcSim(opriskmodel1,100)

### Evaluation of 95% Value-at-Risk
VaR(mc_out,.95)
}
}