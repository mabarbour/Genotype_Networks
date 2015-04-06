source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')
source('~/Documents/Genotype_Networks/Rscripts/network_management_genotype_level.R')

library(ggplot2)
library(ggdendro)

geno.level.net <- genotype_level_interaxn_traits_size[-21,interaxns_noPont] # remove Genotype U and focus on interactions
genotype.labels <- genotype_level_interaxn_traits_size$Genotype[-21]
levels(genotype.labels)[1] <- "C" # change for plotting aesthetics
rownames(geno.level.net) <- genotype.labels

hc.net <- hclust(d = vegdist(geno.level.net, method = "bray")*100, # multiplied by 100 to interpret as percent dissimilarity
                 method = "ave")

## Currently the best dendrogram for plotting the networks and their relationship to each other.
dhc.net <- as.dendrogram(hc.net, hang = 0.05)
#require(ape)
#plot(as.phylo(hc.net), label.offset = 0.01)
#plot(dhc.net, main = "", ylab = "Bray-Curtis Dissimilarity", xlab = "", ylim = c(0,100),  leaflab = "perpendicular")
plot(dhc.net, horiz = TRUE, xlab = "Percent dissimilarity", xlim = c(100,0), dLeaf = -1)


## everything below this is old
dhc.net.df <- dendro_data(dhc.net)

label.positions <- which(segment(dhc.net.df)$x %in% 1:25)
dhc.net.df$segments$label.positions <- ""
dhc.net.df$segments$label.positions[label.positions] <- as.character(dhc.net.df$labels$label)

ggplot(data = segment(dhc.net.df), aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_segment()  +
  geom_text(aes(x = x, y = y - 0.06, label = label.positions)) +
  ylab("Bray-Curtis Dissimilarity") +
  xlab("")
ggdendrogram(dhc.net.df)
#geno.members <- table(tree_level_interaxn_all_plants_traits_size$Genotype)

#plot(hclust(d = dist(geno.level.net)^2, method = "cen", members = geno.members), 
  #   labels = genotype_level_interaxn_traits_size$Genotype, hang = -1)

#plot(hclust(d = dist(geno.level.net), method = "ave", members = geno.members), 
 #    labels = genotype_level_interaxn_traits_size$Genotype, hang = -1)

t <- hclust(d = dist(geno.level.net), method = "ward.D", members = geno.members)#, 
     #labels = genotype_level_interaxn_traits_size$Genotype, hang = -1, leaflab = "perpendicular")
t <- as.dendrogram(t)

hc <- hclust(d = dist(geno.level.net), method = "ward", members = geno.members)
dhc <- as.dendrogram(hc)
ddata <- dendro_data(dhc)

ggdendrogram(ddat)

ggplot(data = segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, y = yend)) #+ 
  #coord_flip() + 
  #scale_y_reverse(expand = c(0.2, 0))
