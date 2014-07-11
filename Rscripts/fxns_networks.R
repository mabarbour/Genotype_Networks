### Functions for Network Analysis

# empty the columns of a food web matrix, but not the rows. This is biologically reasonable when, for example, a lower trophic level is collected (e.g. a gall or leaf mine) but no upper trophic levels were associated with it.
empty.col = function(web){
  # based on "empty" function for bipartite
  cempty <- which(colSums(web) == 0)
  rempty <- which(rowSums(web) == 0)
  cind <- if (length(cempty) == 0) 
    1:NCOL(web)
  else (1:NCOL(web))[-cempty]
  rind <- if (length(rempty) == 0) 
    1:NROW(web)
  else (1:NROW(web))#[-rempty] # attempting to prevent it from removing empty rows.
  out <- web[rind, cind, drop = FALSE]
  return(out)
}

network_metrics <- function(web){
  # this code was adapted from the "bipartite" package's "networklevel" function. I had to create my own code, because some of my food web matrices are too small for the "networklevel" function to handle, and since I would be calculating a lot of these, I wanted it to be repeatable.
  preytot.mat <- matrix(rep(colSums(web), NROW(web)), 
                        NROW(web), byrow = TRUE)
  preyprop.mat <- web/preytot.mat
  predtot.mat <- matrix(rep(rowSums(web), NCOL(web)), 
                        NROW(web), byrow = FALSE)
  predprop.mat <- web/predtot.mat
  H_Nk <- apply(preyprop.mat, 2, function(x) -sum(x * log(x), na.rm = TRUE))
  H_Pk <- apply(predprop.mat, 1, function(x) -sum(x * log(x), na.rm = TRUE))
  n_Nk <- ifelse(colSums(web) != 0, exp(H_Nk),   0)
  n_Pk <- ifelse(rowSums(web) != 0, exp(H_Pk),     0)
  V <- sum(rowSums(web)/sum(web) * n_Pk) # vulnerability
  G <- sum(colSums(web)/sum(web) * n_Nk) # generality
  LD_q <- 0.5 * (V + G) # linkage density
  connect_wt <- LD_q/sum(dim(web)) # weighted connectance
  preyRich <- NROW(web) # prey richness
  predRich <- NCOL(web) # predator richness
  out <- data.frame(V = V, G = G, LD_q = LD_q, connect_wt = connect_wt, preyRich = preyRich, predRich = predRich, totRich = preyRich + predRich)
  return(out)
}

# summarySE from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
