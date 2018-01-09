#library(tea)
#library(h2o) 

fitThreshold <- function(cell, body, tail, method) 
{  
  
if (!is.character(method)) stop("The method must be a character")
if (! method %in% c("dAMSE", "danielsson", "DK", "GH", "gomes", "hall", "Himp", "HW", "mindist", "RT")) stop("The method must be dAMSE, danielsson, DK, GH, gomes, hall, Himp, HW, mindist or RT")
dat=cell$Loss
dat<-dat[which(dat>=0)]
  
if (method=="dAMSE") {thresh<-tea::dAMSE(dat)$threshold}
if (method=="danielsson") {thresh<-tea::danielsson(dat)$threshold} 
if (method=="DK") {thresh<-tea::DK(dat)$threshold} 
if (method=="GH") {thresh<-tea::GH(dat)$threshold}
if (method=="gomes") {thresh<-tea::gomes(dat)$threshold}
if (method=="hall") {thresh<-tea::hall(dat)$threshold}
if (method=="Himp") {thresh<-tea::Himp(dat)$threshold}
if (method=="HW") {thresh<-tea::HW(dat)$threshold}
if (method=="mindist") {thresh<-tea::mindist(dat)$threshold}
if (method=="RT") {thresh<-tea::RT(dat)$threshold[1]}
return(thresh)
}

