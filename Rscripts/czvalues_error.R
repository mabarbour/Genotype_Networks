czvalues_error <- function(web, QuaBiMo_reps = 100) {
  require('bipartite')
  
  observed_web_repeats <- list()
  for(i in 1:QuaBiMo_reps){ 
    observed_web_repeats[[i]] <- web
  }
  
  QuaBiMo_repeats <- sapply(observed_web_repeats, computeModules)
  Q_values <- sapply(QuaBiMo_repeats, function(x) x@likelihood)
  cz_higher_repeats <- sapply(QuaBiMo_repeats, function(x) czvalues(x, weighted = TRUE, level = "higher"))
  cz_lower_repeats <- sapply(QuaBiMo_repeats, function(x) czvalues(x, weighted = TRUE, level = "lower"))
  
  c_higher_df <- ldply(cz_higher_repeats[1, ])
  z_higher_df <- ldply(cz_higher_repeats[2, ])
  
  c_lower_df <- ldply(cz_lower_repeats[1, ])
  z_lower_df <- ldply(cz_lower_repeats[2, ])
  
  # make a list
  list_Qcz.values <- list(Q_values,
                         c_higher_df, 
                         z_higher_df,
                         c_lower_df,
                         z_lower_df)
  names(list_Qcz.values) <- c("Q_values", "c_higher", "z_higher", "c_lower", "z_lower")
  return(list_Qcz.values)
}