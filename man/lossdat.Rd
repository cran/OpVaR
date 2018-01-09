\name{lossdat}
\alias{lossdat}
\docType{data}
\title{
Example loss data set
}
\description{
A minimal loss data set comprising 4 cells (referring to a 2x2 business line/event type operational risk setting).
}
\usage{data("lossdat")}
\format{
  The format is:
List of 4
 $ :'data.frame':	1965 obs. of  3 variables:
  ..$ Loss  : num [1:1965] 1877 1807 918 1480 1218 ...
  ..$ Period: num [1:1965] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ Date  : Date[1:1965], format: "2016-12-30" "2016-12-25" "2016-12-22" ...
 $ :'data.frame':	2025 obs. of  3 variables:
  ..$ Loss  : num [1:2025] 181 610 961 1312 259 ...
  ..$ Period: num [1:2025] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ Date  : Date[1:2025], format: "2016-12-30" "2016-12-30" "2016-12-28" ...
 $ :'data.frame':	1995 obs. of  3 variables:
  ..$ Loss  : num [1:1995] 1334 1067 1068 714 1590 ...
  ..$ Period: num [1:1995] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ Date  : Date[1:1995], format: "2016-12-31" "2016-12-28" "2016-12-26" ...
 $ :'data.frame':	1941 obs. of  3 variables:
  ..$ Loss  : num [1:1941] 650 847 242 513 264 ...
  ..$ Period: num [1:1941] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ Date  : Date[1:1941], format: "2016-12-31" "2016-12-28" "2016-12-22" ...
}
\details{
4 cells with Losses, Time Stamp and Discrete-Time Period (referring to e.g., a quarterly or annual basis). 
}
\source{
Hypothetical, simulated example with bivariate dependencies between loss frequencies and severities
}
\examples{
data(lossdat)
summary(lossdat)
}
\keyword{datasets}
