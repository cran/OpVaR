psevdist <-
function(x,sevdist){
  if(sevdist$type=="plain") return(sevdist$par[[2]](x))
  if(sevdist$type=="spliced") return(pspliced(x,sevdist))
  if(sevdist$type=="mixing") return(pmixing(x,sevdist))
}
