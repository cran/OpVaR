\name{sla}
\alias{sla}
\title{
  Single-Loss Approximation for Operational Value at Risk
}

\description{
  Given the object opriskmodel using the single-loss approximation a list of alpha-quantiles is created for every cell in opriskmodel.
}

\usage{
  sla(opriskmodel, alpha, xi_low = 0.8, xi_high = 1.2, plot = FALSE)
}
\arguments{
  \item{opriskmodel}{
    Object containing the parameters for the severities and the frequencies 
  }
  \item{alpha}{
    Level of the quantile for the total loss-process of the single-loss-approximation
  }
  \item{xi_low}{
    Lower interpolation-point for the spline-function, standard value = 0.8
  }
  \item{xi_high}{
    Upper interpolation-point for the spline-function, standard value = 1.2
  }
  \item{plot}{
    Plot of the interpolated correction term if xi is between xi_low and xi_high, standard value = FALSE
  }
}

\details{
  In the first step the tailindex of the severity distribution is determined. 
  If it does not lie in the critical zone, i.e., not between xi_low and xi_high, the closed-form single-loss approximation from Degen is used.
  If the tailindex lies in the critical zone, then the R-function splinefun with method "hyman" is used to create a spline with given data points xi lying in [xi_low - 0.2, xi_low] and [xi_high, xi_high + 0.2] using the value at xi = 1 as an anchor. 
  A plot of the interpolation is provided such that user-defined adjustment of xi_low and xi_high is possible.
}


\references{
  Bocker, Klaus, and Claudia Kluppelberg. "Operational VaR: a closed-form approximation." Risk-London-Rsik Magazine Limited- 18.12 (2005): 90.
  
  Degen, Matthias. "The calculation of minimum regulatory capital using single-loss approximations." The Journal of Operational Risk 5.4 (2010): 3.
  
}

\author{
  Benedikt Graswald, Jixuan Wang, Christina Zou
}


\examples{
  \donttest{
  #Example: SLA for the spliced log-gamma gpd model (tail-index = 0.014, no interpolation required)
  opriskmodel = list()
  opriskmodel[[1]] = list()
  opriskmodel[[1]]$sevdist = buildSplicedSevdist("lgamma", c(1.23, 0.012), 
                                                 "gpd", c(200, 716, 0.014), 2000, 0.8)
  opriskmodel[[1]]$freqdist = buildFreqdist("pois", 50)
  
  #Example: SLA for the spliced log-gamma gpd model (tail-index = 0.9, interpolation performed)
  opriskmodel[[2]] = list()
  opriskmodel[[2]]$sevdist = buildSplicedSevdist("lgamma", c(1.23, 0.012), 
                                                 "gpd", c(200, 716, 0.9), 2000, 0.8)
  opriskmodel[[2]]$freqdist = buildFreqdist("pois", 50)
  
  sla(opriskmodel, alpha = 0.95) 
  
  #generate plot if interpolation was performed
  sla(opriskmodel, alpha = 0.95, plot = TRUE) 
  
  }
  
}
