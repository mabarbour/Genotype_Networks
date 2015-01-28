#### Andrew MacDonald wrote this code for me.

teammates <- function(modelrun, focal = "vLG_Platy", node.location = "column"){
  node.set <- ifelse(node.location == "column", 2, 1)
  
  interact <- sapply(modelrun[[2]], function(x) x[[node.set]])
  
  has_focal <- sapply(interact, function(x) focal %in% x)
  
  in_mod <- interact[[which(has_focal)]]
  
  in_mod[!(in_mod %in% focal)]
}

