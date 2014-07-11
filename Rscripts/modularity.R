# These modularity metrics are only useful for qualitative data (presence/absence).  The "undirected" function is for a unipartite, undirected network.  The "bipartite" function is for a bipartite, undirected network. This code uses the Newman's spectral algorithm to optimize the modularity value.

# The input for these functions requires a matrix of values of 0 or 1

# I received this code from Rudolf Rohr while I was visiting the Bascompte Lab in the fall of 2012.

split <- function(B){
  
  r <- eigen(B,symmetric=TRUE)
  
  index <- order(r$values,decreasing=TRUE)[1]
  
  v <- r$vectors[,index]
  
  s <- (v>0)*2 -1 
  
  Q <- t(s) %*% B %*% s # t() transposes a matrix
    
  out <- list(s = s, Q = Q)  

  out 
  
}

####################################################################

# allows you to conduct a null model comparing the observed to expected values under random chance.
rel_modularity_bipartite <- function(a, Nsamp=100){
  
  modul <- modularity_bipartite(a)
  modul <- modul$Q
  
  modul.samp <- rep(NA,Nsamp)
  
  p1 <- rep(rowSums(a),times = dim(a)[2])
  dim(p1) <- dim(a)
  
  p2 <- rep(colSums(a),each = dim(a)[1])
  dim(p2) <- dim(a)
  
  p <- 0.5 * (p1 / dim(a)[2] + p2 / dim(a)[1])
  
  for (i in 1:Nsamp){
    
    rand <- runif(prod(dim(a)))
    dim(rand) <- dim(a)
    
    y <- rand < p

    temp <- modularity_bipartite(y)
    
    modul.samp[i] <- temp$Q
    
  }
  
  relmodul <- (modul - mean(modul.samp))/mean(modul.samp) # relative modularity value, if positive (more modularity than you would expect), if negative (less modularity than you would expect). Kind of like an effect size. Note that this is "rQ" in the output.
  
  zscore <- (modul - mean(modul.samp))/sd(modul.samp)  # z-score of comparison to null model
  
  pvalue <- 1-sum(modul>modul.samp)/Nsamp # p-value, relative to the null model.  According to this null model.  These parasitoid communities are not more modular than you would expect by random chance.
  
  out <- list(Q = modul, rQ = relmodul, zscore = zscore, pvalue = pvalue, Qsamp = modul.samp)
  
  out
  
}





####################################################################

# note that a square matrix (A) is needed for this function
modularity_undirected <- function(A){
  
  k <- rowSums(A)
  m <- sum(k)/2
  B <- A - k %*% t(k) / (2*m)  
  
  #first split
  
  r <- split(B)
  Q <- r$Q / (4*m)
  print(Q)
  
  if (Q>1e-4){
  
    groups <- r$s
    groups[groups == -1] <- 2
  
    nb_groups <- 2
    sp <- c(1,1)
    
    #sub split
    
    while (sum(sp)>0){
    
      index <- (1:nb_groups)[sp == 1]
    
      for (i in index){
      
        Bg_temp <- B[groups==i,groups==i]
        Bg <- Bg_temp - diag(rowSums(Bg_temp))
        r <- split(Bg)
        deltaQ <- r$Q / (4*m)
        print(deltaQ)
      
        if (deltaQ>1e-4){
        
          Q <- Q + deltaQ
          nb_groups <- nb_groups + 1
          groups[groups == i] <- r$s -2        
          groups[groups == -1] <- i
          groups[groups == -3] <- nb_groups
          sp[i] <- 1
          sp <- c(sp,1)
        
        } else {sp[i] <- 0}
        
      }
      
    }
    
  } else { groups <- r$S *0
           nb_groups <- 1}
  
  print(sp)
  out <- list(Q = Q, groups=groups, nb_groups=nb_groups)
  out
  
}


####################################################################

modularity_bipartite <- function(A){
  
  S1 <- dim(A)[1]
  S2 <- dim(A)[2]
  
  k1 <- rowSums(A)
  k2 <- colSums(A)
  m <- sum(k1)
  Bs <- A - k1 %*% t(k2) / m
  
  B <- matrix(0,nrow=S1+S2,ncol=S1+S2)
  B[1:S1,(1+S1):(S1+S2)] <- Bs
  B[(S1+1):(S1+S2),1:S1] <- t(Bs)
  
  
  #first split
  
  r <- split(B)
  Q <- r$Q / (4*m)
  
  if (Q>1e-4){
    
    groups <- r$s
    groups[groups == -1] <- 2
    
    nb_groups <- 2
    sp <- c(1,1)
    
    #sub split
    
    while (sum(sp)>0){
      
      index <- (1:nb_groups)[sp == 1]
      
      for (i in index){
        
        Bg_temp <- B[groups==i,groups==i]
        Bg <- Bg_temp - diag(rowSums(Bg_temp))
        r <- split(Bg)
        deltaQ <- r$Q / (4*m)
        
        if (deltaQ>1e-4){
          
          Q <- Q + deltaQ
          nb_groups <- nb_groups + 1
          groups[groups == i] <- r$s -2        
          groups[groups == -1] <- i
          groups[groups == -3] <- nb_groups
          
          if (sum(r$s == 1) > 1) { sp[i] <- 1 } else {sp[i] <- 0}
          if (sum(r$s == -1) > 1) { sp <- c(sp,1) } else {sp <- c(sp,0)}
                    
        } else {sp[i] <- 0}
        
      }
      
    }
    
  } else { groups <- r$S *0
           nb_groups <- 1}
  
  out <- list(Q = as.vector(Q), groups1=groups[1:S1], groups2=groups[(S1+1):(S1+S2)], nb_groups=as.vector(nb_groups))
  out
  
}
