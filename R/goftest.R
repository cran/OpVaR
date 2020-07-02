
goftest=function(cell, object){
  if(class(object) == "sevdist"){
    N = length(cell$Loss)
    sim = rsevdist(N, object)
    res <- suppressWarnings(stats::ks.test(cell$Loss, sim, alternative = "two.sided"))
    
    return(list(res))
    
  }else if(class(object) == "freqdist"){
    periods = unique(cell$Period)
    # count losses per period
    count = table(cell$Period)
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
