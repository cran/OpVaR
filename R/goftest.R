
goftest=function(cell, object){
  if(class(object) == "sevdist"){
    save_sevdist=FALSE
    if(exists("object")) save_sevdist=TRUE; s2=object
    object<<-object
    ddistn <- function(x) dsevdist(x,object)
    pdistn <- function(x) psevdist(x,object)
    qdistn <- function(x) qsevdist(x,object)
    
    Standard1 = goftest::ad.test(cell$Loss, null = pdistn, nullname = "psevdist")
    Standard2 = goftest::cvm.test(cell$Loss, null = pdistn, nullname = "psevdist")
    Standard3 = suppressWarnings(stats::ks.test(cell$Loss, y = pdistn, alternative = "two.sided"))
    
    rm(object,pos=.GlobalEnv)
    if(save_sevdist) object<<-s2 
    
    return(list(Standard1, Standard2, Standard3))
    
  }else if(class(object) == "freqdist"){
    periods = unique(cell$Period)
    # count losses per period
    count = stats::ftable(c(periods, cell$Period), col.vars = 1, 
                          row.vars = NULL)
    count = as.vector(count)
    
    # all possible numbers of losses occurred
    values = seq(0, max(count), by = 1)
    
    # the probability that the respective values are attained according to fitted distribution
    probabilities = evaluate("d", object[[1]], object[[2]], values)
    probabilities = c(probabilities, 1-sum(probabilities))
    
    # count the occurances for each value
    frequencies = rep(0, max(count)+2)
    for(i in 1:length(values)){
      frequencies[i] = length(which(count == values[i]))
    }
    result = suppressWarnings(stats::chisq.test(x = frequencies, p = probabilities))
    return(result)
  }

}
