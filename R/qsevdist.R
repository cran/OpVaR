qsevdist <-
function(x,sevdist){
  if(sevdist$type=="plain") return(sevdist$par[[3]](x))
  if(sevdist$type=="spliced") return(qspliced(x,sevdist))
  if(sevdist$type=="mixing") return(qmixing(x,sevdist))
}
