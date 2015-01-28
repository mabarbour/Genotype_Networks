
# need to calculate z-score of removing each interaction. It is not sufficient, to just remove a node and recalculate modularity to see how it changes, since it is sensitive to the number of nodes, number of links in the web, and the total weights in the web.

node.Q.contribution <- function(web, n_null_webs = 1, Null_model = "swap.web"){
  require('bipartite')
  require('plyr')
  require('dplyr')
  
  upper.node.exclusion.webs <- list()
  for(i in 1:dim(web)[2]){
    upper.node.exclusion.webs[[i]] <- web[ ,-i]
  }
  
  upper.node.excl.modules <- sapply(upper.node.exclusion.webs, computeModules)
  upper.node.Q.values <- sapply(upper.node.excl.modules, function(x) x@likelihood)
  
  upper.node.mod.contribs <- list()
  for(i in 1:dim(web)[2]){
    upper.node.mod.contribs[[i]] <- null_model_analysis_WNODF_QuaBiMo(web = upper.node.exclusion.webs[[i]], type = "QuaBiMo", observed_value = upper.node.Q.values[i], N_null_webs = n_null_webs, null_model = Null_model)
  }
  
  upper.node.mod.contribs.df <- ldply(upper.node.mod.contribs)
  return(upper.node.mod.contribs.df)
}

