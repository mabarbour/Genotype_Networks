module.contribution <- function (moduleWebObject, weighted = FALSE, level = "higher") 
{
  # this function was modified from the czvalues function in the bipartite package. Anything with comments represents a modification from the czvalues function. I kept all naming the same as in the czvalues function just to avoid errors. I did not modify the c-value calculation, but I changed the z-value calculation to estimate the percent contribution of a node to the total number of links in the module. I then changed the output to 'module.percent.contribution' to reflect what it is actually calculating.
  
  require(bipartite) # added because I don't know how to put this into the 'bipartite' environment
  
  #if (!isCorrectModuleWebObject(moduleWebObject)) # I omitted this because it seems that I need to be in the 'bipartite' environment to reference this function
   # stop("This function cannot be applied to this type of object!") 
  zvalues <- function(web, weighted = FALSE) {
    if (weighted) {
      depL <- web/matrix(rowSums(web), nrow = NROW(web), 
                         ncol = NCOL(web), byrow = FALSE)
      k.is <- colSums(depL)
    }
    else {
      k.is <- colSums(web > 0)
    }
    out <- k.is/sum(k.is) # altered this to calculate the % contribution of each node to the number of links in the module. Note that this was only changed for the weighted version.
    out
  }
  modInfo <- listModuleInformation(moduleWebObject)
  modules <- modInfo[[2]]
  nModules <- length(modInfo[[2]])
  if (nModules < 2) 
    stop("This web has no modules.")
  web <- moduleWebObject@originalWeb
  if (level == "lower") 
    web <- t(web)
  if (weighted) {
    depL <- web/matrix(rowSums(web), nrow = NROW(web), ncol = NCOL(web), 
                       byrow = FALSE)
    k.i <- colSums(depL)
  }
  else {
    k.i <- colSums(web > 0)
  }
  z.values.all <- web[1, ]
  k.it.for.t <- matrix(NA, ncol = NCOL(web), nrow = nModules)
  for (t in 1:nModules) {
    h.level.species.names <- if (level == "higher") 
      modules[[t]][[2]]
    else modules[[t]][[1]]
    l.level.species.names <- if (level == "higher") 
      modules[[t]][[1]]
    else modules[[t]][[2]]
    if (weighted) {
      depL.mod1 <- web[l.level.species.names, , drop = FALSE]/matrix(rowSums(web[l.level.species.names, 
                                                                                 , drop = FALSE]), nrow = NROW(web[l.level.species.names, 
                                                                                                                   , drop = FALSE]), ncol = NCOL(web), byrow = FALSE)
      k.it <- colSums(depL.mod1)
    }
    else {
      k.it <- colSums(web[l.level.species.names, , drop = FALSE] > 
                        0)
    }
    k.it.for.t[t, ] <- (k.it/k.i)^2
    z.values.all[h.level.species.names] <- zvalues(web[l.level.species.names, 
                                                       h.level.species.names, drop = FALSE], weighted = weighted)
  }
  c.values.all <- 1 - colSums(k.it.for.t) # note that the calculation for c-values never changed.
  names(c.values.all) <- names(z.values.all)
  z.values.all[is.nan(z.values.all)] <- 0
  return(list(c = c.values.all, module.percent.contribution = z.values.all))
}