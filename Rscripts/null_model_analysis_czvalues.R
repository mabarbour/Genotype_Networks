null_model_analysis_czvalues <- function(web, N_null_webs = 10, null_model = "swap.web", cz_level = "higher") {
  require('bipartite')
  
  # generate null models
  null_web_list <- nullmodel(web, N = N_null_webs, method = null_model)
  
  # add row and column names to each null model
  t <- sapply(null_web_list, provideDimnames)
  
  ### ON THE RIGHT TRACK. NEED TO ADD DIMNAMES TO ALL NULL MODELS, THEN I SHOULD BE ABLE TO CALCULATE THE CZ VALUES.
  browser()
  # need to extract cz values
  #nulls <- sapply(null_web_list, computeModules)

  null_cz.values <- list()
  for(i in 1:length(null_web_list)){
    tmp <- computeModules(null_web_list[[i]])
    browser()
    nulls[[i]] <- czvalues(tmp, weighted = TRUE, level = cz_level)
  }
  #browser()

  
  #null_cz.values <- sapply(nulls, FUN = function(x) czvalues(x, weighted = TRUE, level = cz_level))
  
 
  
  null_values <- vector()
  for(i in 1:dim(nulls)[2]){
    null_values[i] <- nulls[ ,i]$statistic[3] 
  }
  
  mean_null_values <- mean(null_values)
  sd_null_values <- sd(null_values)
  
  z_score <- (observed_value - mean_null_values)/sd_null_values
  names(z_score) <- "z-score"
  names(mean_null_values) <- "mean of null values"
  names(sd_null_values) <- "SD of null values"
  return(c(z_score, mean_null_values, sd_null_values))
}

