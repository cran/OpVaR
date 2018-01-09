rsevdist <-
function(x,sevdist){
  if(sevdist$type=="plain") return(sevdist$par[[4]](x))
  if(sevdist$type=="spliced") return(rspliced(x,sevdist))
  if(sevdist$type=="mixing") return(rmixing(x,sevdist))
}
