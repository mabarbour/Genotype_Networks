null_model_analysis_czvalues <- function(web, N_null_webs = 10, null_model = "swap.web", ColNames = colnames(web), RowNames = rownames(web)) {
  require('bipartite')
  require('plyr')
  
  # generate null models
  null_web_list <- nullmodel(web, N = N_null_webs, method = null_model)
  
  # add row and column names to each null model. Necessary to run czvalues function on moduleWebObject
  for(i in 1:N_null_webs){
    colnames(null_web_list[[i]]) <- ColNames
    rownames(null_web_list[[i]]) <- RowNames
  }
  
  # compute modules for each null model
  nulls <- sapply(null_web_list, computeModules)

  # compute "higher" trophic level cz values for each null model
  high.nulls_cz.values <- sapply(nulls, function(x) module.contribution(x, weighted = TRUE, level = "higher"))
  
  
 
  # extract c values
  high.null_c.values <- list()
  for(i in 1:N_null_webs){
    high.null_c.values[[i]] <- high.nulls_cz.values[ ,i]$c 
  }
  
  # extract z values
  high.null_z.values <- list()
  for(i in 1:N_null_webs){
    high.null_z.values[[i]] <- high.nulls_cz.values[ ,i]$module.percent.contribution
  }
  
  # summary stats for c-values
  high.null_c.values.df <- ldply(high.null_c.values)
  
  # summary stats for z-values
  high.null_z.values.df <- ldply(high.null_z.values)
  
  
  # compute "lower" trophic level cz values for each null model
  low.nulls_cz.values <- sapply(nulls, function(x) module.contribution(x, weighted = TRUE, level = "lower"))
  
  # extract c values
  low.null_c.values <- list()
  for(i in 1:N_null_webs){
    low.null_c.values[[i]] <- low.nulls_cz.values[ ,i]$c 
  }
  
  # extract z values
  low.null_z.values <- list()
  for(i in 1:N_null_webs){
    low.null_z.values[[i]] <- low.nulls_cz.values[ ,i]$module.percent.contribution
  }
  
  # summary stats for c-values
  low.null_c.values.df <- ldply(low.null_c.values)
  
  # summary stats for z-values
  low.null_z.values.df <- ldply(low.null_z.values)
  
  
  # make a list
  list_cz.values <- list(high.null_c.values.df, high.null_z.values.df, low.null_c.values.df, low.null_z.values.df)
  names(list_cz.values) <- c("high.null_c.values","high.null_z.values", "low.null_c.values","low.null_z.values")
  return(list_cz.values)
}

