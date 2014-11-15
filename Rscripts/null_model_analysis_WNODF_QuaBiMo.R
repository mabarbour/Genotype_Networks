null_model_analysis_WNODF_QuaBiMo <- function(web, type, observed_value, N_null_webs = 10000, null_model = "swap.web") {
  require('bipartite')
  source('~/Documents/vegan/R/nestednodf.R')
  
  # generate null models
  null_web_list <- nullmodel(web, N = N_null_webs, method = null_model)
  
  # WNODF analysis
  if(type == "WNODF"){
    
    nulls <- sapply(null_web_list, FUN = function(x) nestednodf(x, order = TRUE, weighted = TRUE))
    
    null_values <- vector()
    for(i in 1:dim(nulls)[2]){
      null_values[i] <- nulls[ ,i]$statistic[3]
    }
  }
  
  # QuaBiMo analysis
  if(type == "QuaBiMo"){
    
    nulls <- sapply(null_web_list, computeModules)
    
    null_values <- sapply(nulls, function(x) x@likelihood)  
  }
  
  mean_null_values <- mean(null_values)
  sd_null_values <- sd(null_values)
  
  z_score <- (observed_value - mean_null_values)/sd_null_values
  names(z_score) <- "z-score"
  names(mean_null_values) <- "mean of null values"
  names(sd_null_values) <- "SD of null values"
  return(c(z_score, mean_null_values, sd_null_values))
}