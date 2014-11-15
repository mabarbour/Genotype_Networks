#library(devtools)

#install_github("tpoisot/betalink")

library(betalink)
library(reshape2)
library(plyr)
library(dplyr)

# Comparing 2 networks whose dissimilarity is entirely driven by "rewiring" of species interactions
w1 = matrix(c(1, 0, 0, 
              0, 1, 0, 
              0, 0, 1), ncol=3)
w2 = matrix(c(0, 1, 1, 
              1, 0, 1,
              1, 1, 0), ncol=3)

colnames(w1) = c('p1', 'p2', 'p3')
rownames(w1) = c('h1', 'h2', 'h3')
colnames(w2) = colnames(w1)
rownames(w2) = rownames(w1)

realizations = list(w1, w2)

metaweb_1 = metaweb(W = realizations)

betalink(w1, w2, bf = B01)
betalink.plot(m1 = w1, m2 = w2, by.unique = TRUE)

w1_trans <- melt(w1)
w1_trans <- mutate(w1_trans, site = rep("w1", dim(w1_trans)[1]))

w2_trans <- melt(w2)
w2_trans <- mutate(w2_trans, site = rep("w2", dim(w2_trans)[1]))

w1_w2_comp <- rbind.data.frame(w1_trans, w2_trans)
w1_w2_comp_cast <- dcast(w1_w2_comp, site ~ Var1 + Var2, sum)

library(vegan)
betadiver(w1_w2_comp_cast[ ,-1], method = "w")

# Comparing 2 networks that differ solely due to species turnover. Note that there has to be at least 1 overlapping link for betalink to calculate species turnover...does that make sense!?!?!
w3 = matrix(c(1, 0, 0, 
              0, 1, 0, 
              0, 0, 1), ncol=3)

colnames(w3) = c('p4', 'p5', 'p6') # play around with the number of overlapping species
rownames(w3) = c('h3', 'h4', 'h5') # play around with the number of overlapping species

betalink(w1, w3, B01)
betalink.plot(w1, w3)

realization_2 = list(w1, w3)
betalink.dist(realization_2, bf = B01)
beta.os_prime(W = realizations)

w3_trans <- melt(w3)
w3_trans <- mutate(w3_trans, site = rep("w3", dim(w3_trans)[1]))

w1_w3_comp <- rbind.data.frame(w1_trans, w3_trans)
w1_w3_comp_cast <- dcast(w1_w3_comp, site ~ Var1 + Var2, sum)

betadiver(w1_w3_comp_cast[ ,-1], method = "w") # note that calculating dissimilarity this way actually reflects a realistic dissimilarity in interactions between sites!

# Comparing 2 networks with a small number of interactions
w4 = matrix(c(1), ncol = 1)
colnames(w4) <- "p1"
rownames(w4) <- "h1"

betalink(w1,w4, bf=B01)
betalink.plot(w1, w4)

realizations_2 = list(w1, w4)
metaweb_2 = metaweb(realizations_2)

betalink.dist(realizations_2)
