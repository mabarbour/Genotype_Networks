library(devtools)

install_github("tpoisot/betalink")

library(betalink)

w1 = matrix(c(1, 0, 0, 1), ncol=2)
w2 = matrix(c(1, 1, 0, 1), ncol=2)
w3 = matrix(c(1, 0, 1, 1), ncol = 2)
colnames(w1) = c('p1', 'p2')
rownames(w1) = c('h1', 'h2')
colnames(w2) = colnames(w1)
rownames(w2) = rownames(w1)
colnames(w3) = colnames(w1)
rownames(w3) = rownames(w1)

realizations = list(w1, w2, w3)

metaweb_1 = metaweb(W = realizations)

betalink(w1, w2, bf = B01)
betalink.plot(m1 = w1, m2 = w2, by.unique = TRUE)
beta

betalink.dist(realizations, bf = B01)
beta.os_prime(W = realizations)


rabie.test = function(w1,w2,bf="jaccard"){
  sp1 = list(top=rowSums(w1),bottom=colSums(w1),all=c(rowSums(w1),colSums(w1)[!names(colSums(w1))%in%names(rowSums(w1))]))
  sp2 = list(top=rowSums(w2),bottom=colSums(w2),all=c(rowSums(w2),colSums(w2)[!names(colSums(w2))%in%names(rowSums(w2))]))
  topCom=vec2data.frame(sp1$top,sp2$top,all=T)
  bottomCom=vec2data.frame(sp1$bottom,sp2$bottom,all=T)
  allCom=vec2data.frame(sp1$all,sp2$all,all=T)  
  topCom[is.na(topCom)]=0
  bottomCom[is.na(bottomCom)]=0
  allCom[is.na(allCom)]=0
  beta_U = vegdist(topCom,method=bf,diag=T,upper=T)
  beta_L = vegdist(bottomCom,method=bf,diag=T,upper=T)
  beta_S = vegdist(allCom,method=bf,diag=T,upper=T)
  # Common species
  Csp = names(sp1$all)[names(sp1$all) %in% names(sp2$all)]
  CUsp = names(sp1$top)[names(sp1$top) %in% names(sp2$top)]
  CLsp =names(sp1$bottom)[names(sp1$bottom) %in% names(sp2$bottom)]
  return(list(U = beta_U[1], L = beta_L[1], S = beta_S[1]))
}
