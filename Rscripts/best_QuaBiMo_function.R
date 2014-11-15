best_QuaBiMo <- function(web, QuaBiMo_reps = 100) {
  require('bipartite')
  
  observed_web_repeats <- list()
  for(i in 1:QuaBiMo_reps){ 
    observed_web_repeats[[i]] <- web
  }
  
  QuaBiMo_repeats <- sapply(observed_web_repeats, computeModules)
  best_observed_value <- max(sapply(QuaBiMo_repeats, function(x) x@likelihood)) # takes the maximum out of a user specified number of replications, since this is an optimization algorithm and we want to avoid reaching a local optimum.
  
  best_reps <- which.max(sapply(QuaBiMo_repeats, function(x) x@likelihood)) # identify best partition of modularity algorithm
  
  output <- list(best_observed_value = best_observed_value,
                 best_reps = best_reps,
                 best_module_info = QuaBiMo_repeats[best_reps])
}

