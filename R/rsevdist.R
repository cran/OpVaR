rsevdist <-
function(n,sevdist){
  if(sevdist$type=="plain") return(sevdist$par[[4]](n))
  if(sevdist$type=="spliced") return(rspliced(n,sevdist))
  if(sevdist$type=="mixing") return(rmixing(n,sevdist))
}
