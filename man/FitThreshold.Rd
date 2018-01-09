\name{fitThreshold}
\alias{fitThreshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Threshold estimation for spliced distribution
}

\description{
  Given a dataset with a chosen distribution for the data in the body and another distribution in the tail, the threshold is estimated. 
}

\usage{
  fitThreshold(cell, body, tail, method)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cell}{
    list containing the data in the component Loss
  }
 \item{body}{
    list containing the data in the component Loss
  }
 \item{tail}{
    list containing the data in the component Loss
  }
  \item{method}{
    Method for the threshold estimation. In case of a GPD tail, it can be chosen between "dAMSE", "danielsson", "DK", "GH", "gomes", "hall", "Himp", "HW", "mindist", "RT"
  }  
}

\details{
  %%  ~~ If necessary, more details than the description above ~~
## GPD tail
In case of a Generalized Pareto Function in the tail, there are several methods to estimate the threshold which are based on extreme value statistics.
Most of them are based on the Asymptotic Mean Squared Error (AMSE) approach.
where the asymptotic squared error \deqn{Asy E((Hill(N,k)-gamma)^2)} between the true inverse tail index gamma and the Hill-estimator is minimized: "dAMSE"", "danielsson"", "DK", "GH"", gomes", "hall", "Himp" and "HW". For further details, see tea.
"RT" is a heuristic approach where the smallest order statistic k that minimizes the expression \eqn{1/k sum_i=1^k i^beta |gamma_i-median(gamma_1,...,gamma_k)|} is chosen for the threshold. 
"mindist" minimizes the distance between the largest upper order statistics (the empirical tail) and the theoretical tail of a Pareto distribution. For further details, see tea.

}

\value{
Returns a numeric value giving threshold fitted according the given method.
}
\references{
  %% ~put references to the literature/web site here ~
  # GPD Tail:
    
    Caeiro, F. and Gomes, M.I. (2014) On the bootstrap methodology for the estimation of the tail
  sample fraction. Proceedings of COMPSTAT, 545-552.
  
    Caeiro, J. and Gomes, M.I. (2016) Threshold selection in extreme value analysis. Extreme Value
  Modeling and Risk Analysis:Methids and Applications, 69-86
  
    Danielsson, J. and Ergun, L.M. and de Haan, L. and de Vries, C.G. (2016) Tail Index Estimation:
  Quantile Driven Threshold Selection
  
    Danielsson, J. and Haan, L. and Peng, L. and Vries, C.G. (2001) Using a bootstrap method to choose
  the sample fraction in tail index estimation. Journal of Multivariate analysis, 2, 226-248.
  
    Drees, H. and Kaufmann, E. (1998) Selecting the optimal sample fraction in univariate extreme
  value estimation. Stochastic Processes and their Applications, 75(2), 149-172.
  
    Gomes, M.I. and Figueiredo, F. and Neves, M.M. (2012) Adaptive estimation of heavy right tails:
  resampling-based methods in action. Extremes, 15, 463-489
  
    Guillou, A. and Hall, P. (2001) A Diagnostic for Selecting the Threshold in Extreme Value Analysis
  Journal of the Royal Statistical Society, 63(2), 293-305
  
    Hall, P. (1990) Using the Bootstrap to Estimate Mean Squared Error and Select Smoothing Parameter
  in Nonparametric Problems. Journal of Multivariate Analysis, 32, 177-203

    Hall, P. and Welsh, A.H. (1985) Adaptive estimates of parameters of regular variation. The Annals
  of Statistics, 13(1), 331-341.
  
    Reiss, R.-D. and Thomas, M. (2007) Statistical Analysis of Extreme Values: With Applications to
  Insurance, Finance, Hydrology and Other Fields. Birkhauser, Boston.
}

\author{
  Christina Zou, Leonie Wicht
}


\examples{
data(lossdat)
fitThreshold(lossdat[[1]], method="dAMSE")
  
}
