\name{fitWeights}
\alias{fitWeights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fitting the weights of the body and the tail for a spliced distribution
}

\description{
  Fitting the weights of the spliced distribution.
}

\usage{
  fitWeights(cell, thresh)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cell}{
    List containing the data in the component Loss
  }
  \item{thresh}{
    Threshold between the body and the tail.
  }
}

\details{
  %%  ~~ If necessary, more details than the description above ~~
    Given the threshold t, this function fits the weight on the body distribution as part of fitting the parameters of the spliced distribution.
  The weights depends only on the number of data points below and above the threshold, not on the distribution functions in the body and the tail.
  }
  
  \author{
    Christina Zou, Leonie Wicht
  }
  
  
  \examples{
    data(lossdat)
    thresh<-quantile(lossdat[[1]]$Loss, 0.9)
    w<-fitWeights(lossdat, thresh)
    
  }
  
  \seealso{
    fitThreshold, fitSpliced
  }
  
  % Add one or more standard keywords, see file 'KEYWORDS' in the
  % R documentation directory.