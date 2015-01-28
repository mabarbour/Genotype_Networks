source('~/Documents/Genotype_Networks/Rscripts/teammates.R')

teammates_table <- function(web, focal = "vLG_Platy", node.location = "column", QuaBiMo_reps = 100) { # focal
  require('bipartite')
  
  # run modularity algorithm and create a list of the module information
  observed_web_repeats <- list()
  for(i in 1:QuaBiMo_reps){ 
    observed_web_repeats[[i]] <- web
  }
  
  QuaBiMo_repeats <- sapply(observed_web_repeats, computeModules)  
  modlist <- lapply(QuaBiMo_repeats, listModuleInformation)
  
  run_teammates <- lapply(modlist, function(x) teammates(x, focal = focal, node.location = node.location))
  
  total <- do.call(c, run_teammates)
  
  output <- table(total)
}
  