source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')
source('~/Documents/Genotype_Networks/Rscripts/network_management_genotype_level.R')

library(ggplot2)
library(ggdendro)

geno.level.net <- genotype_level_interaxn_traits_size[ ,interaxns_noPont]

geno.members <- table(tree_level_interaxn_all_plants_traits_size$Genotype)

plot(hclust(d = dist(geno.level.net)^2, method = "cen", members = geno.members), 
     labels = genotype_level_interaxn_traits_size$Genotype, hang = -1)

plot(hclust(d = dist(geno.level.net), method = "ave", members = geno.members), 
     labels = genotype_level_interaxn_traits_size$Genotype, hang = -1)

t <- hclust(d = dist(geno.level.net), method = "ward", members = geno.members)#, 
     #labels = genotype_level_interaxn_traits_size$Genotype, hang = -1, leaflab = "perpendicular")
t <- as.dendrogram(t)

hc <- hclust(d = dist(geno.level.net), method = "ward", members = geno.members)
dhc <- as.dendrogram(hc)
ddata <- dendro_data(dhc)

ggplot(data = segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, y = yend)) #+ 
  #coord_flip() + 
  #scale_y_reverse(expand = c(0.2, 0))
