
#################################        Auxiliary functions       #################################

## Retrieve the tail index from a sevdist-object
tailindex_plain <- function(sevdist){
  return(switch(sevdist$par[[6]], 
                "gpd" = sevdist$par[[5]][3],
                "lgamma" = 1 / sevdist$par[[5]][2],
                "gh" = 1 / sevdist$par[[5]][2], 
                "lnorm" = 0,
                "weibull" = 0))
}

tailindex_spliced <- function(sevdist){
  return(switch(sevdist$par[[1]][[6]], 
                "gpd" = sevdist$par[[1]][[5]][3], 
                "lgamma" = 1 / sevdist$par[[1]][[5]][2], 
                "gh" = 1 / sevdist$par[[1]][[5]][2], 
                "lnorm" = 0, 
                "weibull" = 0))
}

tailindex_mixing <- function(sevdist){
  return(switch(sevdist$par[[1]][[6]], 
                "gpd" = sevdist$par[[1]][[5]][3], 
                "lgamma" = 1 / sevdist$par[[1]][[5]][2], 
                "gh" = 1 / sevdist$par[[1]][[5]][2], 
                "lnorm" = 0, 
                "weibull" = 0))
}


## Calculate a temporary sevdist from xi
sevdist_xi_plain <- function(xi, sevdist){
  sevdist_tmp <- sevdist
  # new parameters
  if(sevdist$par[[6]] == "gpd"){
    sevdist_tmp$par[[5]][3] <- xi
  } else if(sevdist$par[[6]] %in% c("lgamma", "gh")){
    sevdist_tmp$par[[5]][2] <- 1/xi
  }
  
  # new cdf
  sevdist_tmp$par[[2]] <- function(q) {evaluate("p", sevdist_tmp$par[[6]], sevdist_tmp$par[[5]], q)}
  # new quantile function
  sevdist_tmp$par[[3]] <- function(p) {evaluate("q", sevdist_tmp$par[[6]], sevdist_tmp$par[[5]], p)}
  return(sevdist_tmp)
}

sevdist_xi_spliced <- function(xi, sevdist){
  sevdist_tmp <- sevdist
  # new parameters
  if(sevdist$par[[1]][[6]] == "gpd"){
    sevdist_tmp$par[[1]][[5]][3] <- xi
  } else if(sevdist$par[[1]][[6]] %in% c("lgamma", "gh")){
    sevdist_tmp$par[[1]][[5]][2] <- 1/xi
  }
  
  # d/p of underlying plain distribution with new tail index:
  pdisttail <- function(q) {do.call(paste0("p", sevdist_tmp$par[[1]][[6]]), c(list(q), sevdist_tmp$par[[1]][[5]]))}
  qdisttail <- function(p) {do.call(paste0("q", sevdist_tmp$par[[1]][[6]]), c(list(p), sevdist_tmp$par[[1]][[5]]))}
  # new cdf in tail
  sevdist_tmp$par[[1]][[2]] <- function(q) {ifelse(q>sevdist_tmp$thresh, (pdisttail(q) - pdisttail(sevdist_tmp$thresh))/(1-pdisttail(sevdist_tmp$thresh)), 0)}
  # new quantile function in tail:
  sevdist_tmp$par[[1]][[3]] <- function(p) {qdisttail(pdisttail(sevdist_tmp$thresh) + p * (1 - pdisttail(sevdist_tmp$thresh)))}
  return(sevdist_tmp)
}

sevdist_xi_mixing <- function(xi, sevdist){
  sevdist_tmp <- sevdist
  if(sevdist$par[[1]][[6]] == "gpd"){
    sevdist_tmp$par[[1]][[5]][3] <- xi
  } else if(sevdist$par[[1]][[6]] %in% c("lgamma", "gh")){
    sevdist_tmp$par[[1]][[5]][2] <- 1/xi
  }
  
  # d/p of underlying plain distribution with new tail index:
  pdisttail <- function(q) {do.call(paste0("p", sevdist_tmp$par[[1]][[6]]), c(list(q), sevdist_tmp$par[[1]][[5]]))}
  qdisttail <- function(p) {do.call(paste0("q", sevdist_tmp$par[[1]][[6]]), c(list(p), sevdist_tmp$par[[1]][[5]]))}
  # new cdf in tail
  sevdist_tmp$par[[1]][[2]] <- function(q) {evaluate("p", sevdist$par[[1]][[6]], sevdist_tmp$par[[1]][[4]], q)}
  # new quantile function in tail:
  sevdist_tmp$par[[1]][[3]] <- function(p) {evaluate("q", sevdist$par[[1]][[6]], sevdist_tmp$par[[1]][[4]], p)}
  return(sevdist_tmp)
}




correction_high <- function(xi, sevdist, type, mean_freq, disp_freq, alpha){
  sevdist_xi <- do.call(paste("sevdist_xi", as.character(type), sep = "_"), list(xi, sevdist))
  basic_SLA <- qsevdist(1 - (1-alpha)/mean_freq, sevdist_xi)
  c_xi <- (1-xi) /2 /gamma(1-2/xi) * (gamma(1-1/xi)^2)
  return((mean_freq + disp_freq -1) * (1-alpha) / mean_freq * basic_SLA * c_xi / (1-1/xi))
}

correction_low <- function(xi, sevdist, type, mean_freq, disp_freq){
  sevdist_xi <- do.call(paste("sevdist_xi", as.character(type), sep = "_"), list(xi, sevdist))
  mean_sev <- try(pracma::quadinf(function(x){1-psevdist(x, sevdist_xi)}, xa = 0, xb = Inf)$Q, silent=TRUE)
  return((mean_freq + disp_freq -1) * mean_sev)
}

correction_1 <- function(xi, sevdist, type, mean_freq, disp_freq, alpha){
  sevdist_xi <- do.call(paste("sevdist_xi", as.character(type), sep = "_"), list(xi, sevdist))
  basic_SLA <- qsevdist(1 - (1-alpha)/mean_freq, sevdist_xi)
  return((mean_freq + disp_freq -1) * pracma::integral(function(x){1-psevdist(x, sevdist_xi)}, xmin = 0, xmax = basic_SLA, method = "Kronrod"))
}


# Spline interpolation
spline_SLA <- function(sevdist, xi_low, xi_high, type, mean_freq, disp_freq, alpha){
  # define the points used for calculating the spline
  grid_xi_high <- seq(xi_high, xi_high + 0.2, 0.05)
  grid_xi_low <- seq(xi_low - 0.2, xi_low, 0.05)
  grid_xi <- c(grid_xi_low, 1, grid_xi_high)
  
  # evaluate the correction term outside the singularity area
  spline_high <- sapply(grid_xi_high, correction_high, sevdist, type, mean_freq, disp_freq, alpha)
  spline_low <- sapply(grid_xi_low, correction_low, sevdist, type, mean_freq, disp_freq)
  spline_1 <- correction_1(1, sevdist, type, mean_freq, disp_freq, alpha)
  spline_values <- c(spline_low, spline_1, spline_high)
  
  # calculate the spline function
  correct_spline <- try(splinefun(grid_xi, spline_values, method = "hyman"), silent=TRUE)
  return(correct_spline)
}



#############################         SLA with spline        ###############################
sla <- function(opriskmodel, alpha, xi_low = 0.8, xi_high = 1.2, plot = FALSE){
  # initialize the list containing the approximated VaR for each cell
  approx <- list()
  # iterate over the cells
  for(cell in 1:length(opriskmodel)){
    approx[[cell]] <- list()
    approx[[cell]]$interpolation <- FALSE
    
    # track the progress
    #print(paste("Cell ", cell))
    
    # type of current cell
    cell_type <- (opriskmodel[[cell]]$sevdist)$type
    
    # mean loss frequency of current cell
    mean_freq <- switch(opriskmodel[[cell]]$freqdist[[1]], 
                        "pois" = opriskmodel[[cell]]$freqdist[[2]], 
                        "nbinom" =  opriskmodel[[cell]]$freqdist[[2]][1]*(1-opriskmodel[[cell]]$freqdist[[2]][2])/opriskmodel[[cell]]$freqdist[[2]][2])
    #print(paste("mean_freq ", mean_freq))
    
    # dispersion factor of the frequency distribution of current cell
    disp_freq <- switch(opriskmodel[[cell]]$freqdist[[1]], 
                        "pois" = 1, 
                        "nbinom" = 1/opriskmodel[[cell]]$freqdist[[2]][2])
    #print(paste("disp_freq ", disp_freq))
    
    # original SLA 
    simple_VaR <- qsevdist(1 - (1-alpha)/mean_freq,  opriskmodel[[cell]]$sevdist)
    
    # test if mean severity is finite
    mean_sev <- try(pracma::quadinf(function(x){x*dsevdist(x, opriskmodel[[cell]]$sevdist)}, xa = 0, xb = Inf)$Q, silent=FALSE)
    
    # retrieve the tail index from the tail severity distribution
    tailindex <- do.call(paste("tailindex", as.character(cell_type), sep = "_"), list(opriskmodel[[cell]]$sevdist))
    
    if(!is.null(tailindex)){
      # calculate the correction term for the quantile of the total-loss distribution
      if(tailindex <= xi_low){    # case of finite severity mean?!?!
        correction <- (mean_freq + disp_freq -1) * mean_sev
      } else if(tailindex >= xi_high){
        c_xi <- (1-tailindex) /2 /gamma(1-2/tailindex) * (gamma(1-1/tailindex)^2)
        correction <- (mean_freq + disp_freq -1) * (1-alpha) / mean_freq * simple_VaR * c_xi / (1-1/tailindex)
      } else if(tailindex == 1){
        correction <- (mean_freq + disp_freq -1) * (pracma::integral(function(x){1-psevdist(x, opriskmodel[[cell]]$sevdist)}, xmin = 0, xmax = simple_VaR, method = "Kronrod"))
      } else {   # xi lies in the sigularity area ==> need monotone spline interpolation
        correction_spline <- spline_SLA(opriskmodel[[cell]]$sevdist, xi_low, xi_high, cell_type, mean_freq, disp_freq, alpha)
        xi_high_tmp <- xi_high
        xi_low_tmp <- xi_low
        # track the progress
        #print(paste("current sigularity zone: ", xi_low_tmp, "-", xi_high_tmp))
        
        while(class(correction_spline) == "try-error"){
          xi_high_tmp <- xi_high_tmp + 0.05
          xi_low_tmp <- xi_low_tmp - 0.05
          correction_spline <- spline_SLA(opriskmodel[[cell]]$sevdist, xi_low_tmp, xi_high_tmp, cell_type, mean_freq, disp_freq, alpha)
          # track the progress
          #print(paste("current sigularity zone: ", xi_low_tmp, "-", xi_high_tmp))
        }
        
        # evaluate the spline at xi of interest
        correction <- correction_spline(tailindex)
        approx[[cell]]$interpolation <- TRUE
      }
      
      # Plot:
      if(approx[[cell]]$interpolation == TRUE && plot == TRUE) {
        # Plot the divergence
        grid_high <- seq(1.005, 1.5, 0.01)
        grid_low <- seq(0.5, 0.995, 0.01)
        grid <- c(grid_low, 1, grid_high)
        
        correction_high_values <- sapply(grid_high, correction_high, opriskmodel[[cell]]$sevdist, cell_type, mean_freq, disp_freq, alpha)
        correction_low_values <- sapply(grid_low, correction_low, opriskmodel[[cell]]$sevdist, cell_type, mean_freq, disp_freq)
        
        plot(grid_low, correction_low_values, 
             xlim=c(0.5, 1.5), 
             ylim=c(min(correction_low_values, correction_high_values), max(correction_low_values, correction_high_values)),
             type="l", col="royalblue", lwd = 2,
             xlab = expression(xi), ylab = "SLA correction term")
        lines(grid_high, correction_high_values, col="royalblue", lwd = 2)
        abline(v = 1, lty = "dashed")
        
        # Plot the spline
        if(!exists("xi_low_tmp")){
          xi_low_tmp <- xi_low
          xi_high_tmp <- xi_high
        }
        curve(correction_spline, from = xi_low_tmp, to = xi_high_tmp, lty = "dashed", 
              col="darkorange4", lwd = 2, add = TRUE)
        points(1, correction_spline(1), col = "darkorange", pch = 19)
      }
      
    }
    else{
      correction <- (mean_freq + disp_freq -1) * mean_sev
      if(!is.numeric(correction)) correction <- 0
      approx[[cell]]$interpolation <- FALSE
      
    }
    
    
    # catch other unexpected errors :correction <- 0
    
    # append the SLA of current cell to the output-list
    
    approx[[cell]]$value <- simple_VaR + correction
  }
  output = NULL
  for(i in 1:length(approx)){
    output = rbind(output, c(approx[[i]]$value, as.character(approx[[i]]$interpolation)))
  }
  cat("OpRiskmodel - SLA: ", alpha*100, "% \n")
  cat("----------------------------------------------\n")
  colnames(output) =  c("VaR", "Interpolation")
  rownames(output) = paste("Cell", c(seq(1, length(approx), by=1)), ": ")
  prmatrix(output, quote = FALSE)
}
