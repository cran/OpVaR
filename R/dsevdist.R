dsevdist <-
function(x,sevdist){
  if(sevdist$type=="plain") return(sevdist$par[[1]](x))
  if(sevdist$type=="spliced") return(dspliced(x,sevdist))
  if(sevdist$type=="mixing") return(dmixing(x,sevdist))
}
