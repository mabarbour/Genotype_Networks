
library(reshape)
library(reshape2)
library(plyr)
library(dplyr)
library(vegan)


# upload molten gall network data
gall_net_melt <- read.csv("~/Documents/Genotype_Networks/data/gall_network_data.csv")
gall_net_melt <- tbl_df(gall_net_melt)

##### Does the number of shoots sampled vary among genotypes?

# data frame
shootEsts_df <- gall_net_melt %>%
  mutate(plant.position = as.factor(plant.position)) %>%
  group_by(Gender, Genotype, plant.position) %>%
  summarise(shootEst.no18 = mean(shootEst.no18), shootEst.all = mean(shootEst.all), galls.found = mean(Galls.found)) # needed to take the mean because these shoot estimates are duplicated throughout the molten dataframe
shootEsts_df <- mutate(shootEsts_df, plant.position = as.character(plant.position))

shootEsts_genotype_df <- shootEsts_df %>%
  group_by(Genotype) %>%
  summarise(n = n(), shootEst.no18 = sum(shootEst.no18), shootEst.all = sum(shootEst.all))


##### How dissimilar are gall-parasitoid interaction networks among genotypes?

### Genotype-level networks
levels(gall_net_melt$gall_contents)

genotype_net_filter <- gall_net_melt %>%
  filter(gall_contents %in% c("Eulo.fem", "Eulo.mal", "Eury.fem", "Eury.mal", "Lathro.fem", "Lathro.mal", "Lestodip", "Mesopol", "Mymarid", "Platy", "Ptero.2", "Tory.fem", "Tory.mal")) # removed survival data for gall species

genotype_net_add <- mutate(genotype_net_filter, gall_contents_collapse = revalue(gall_contents, c("Eulo.fem" = "Eulo", "Eulo.mal" = "Eulo", "Eury.fem" = "Eury", "Eury.mal" = "Eury", "Lathro.fem" = "Lathro", "Lathro.mal" = "Lathro", "Ptero.2" = "Mesopol", "Tory.fem" = "Tory", "Tory.mal" = "Tory") ))


genotype_net <- cast(genotype_net_add, gall.sp ~ gall_contents_collapse | Genotype, sum)
genotype_net <- genotype_net[-21] # never any interactions for Genotype U

shootEsts_genotype_df_no_U <- shootEsts_genotype_df %>% filter(Genotype != "U") # update corresponding shoot estimate data.

# this loop prepares each network for analysis in betalink
for(i in 1:length(genotype_net)) {
  rownames(genotype_net[[i]]) <- as.character(genotype_net[[i]]$gall.sp) # change row names to gall specie names
  genotype_net[[i]] <- as.data.frame(select(genotype_net[[i]], -gall.sp)) # remove gall.sp as a column of data. Turning this into a matrix is important for easier analysis with betalink package
  genotype_net[[i]] <- round(genotype_net[[i]]/shootEsts_genotype_df_no_U$shootEst.no18[i]*1000) # scale it to the density of interactions per 1000 shoots, rounded to the nearest whole integer.
}


### Beta-diversity of interaction network analysis
source('~/Documents/betalink/R/betalink.R')
source('~/Documents/betalink/R/vec2data.frame.R')
source('~/Documents/betalink/R/measures.r')
source('~/Documents/betalink/R/betalink.b.r')
source('~/Documents/betalink/R/betalink.dist.r')
source('~/Documents/betalink/R/betalink.q.R')

betalink_genotype_net_qualitative <- betalink.dist(genotype_net, bf=B01, triangular = T) # tried using Sorenson's index (B11), but it was giving very wonky results. Also giving wonky results when I try to use Paul Rabie's qualitative code...

mean(betalink_genotype_net_qualitative$WN) # average dissimilarity is 65% among genotypes
mean(betalink_genotype_net_qualitative$OS) # 0.06
mean(betalink_genotype_net_qualitative$S_all.taxa) # 0.48
mean(betalink_genotype_net_qualitative$S_upper) # 0.40
mean(betalink_genotype_net_qualitative$S_lower) # 0.54
mean(betalink_genotype_net_qualitative$ST) # 0.59
mean(betalink_genotype_net_qualitative$contrib) # 87% of the dissimilarity is due to species turnover.
hist(betalink_genotype_net_qualitative$contrib) # wow, a lot of 100% dissimilarity based on species composition. Apparent, bimodal distribution.

# even after removing combinations with absolutely no species in common "S_all.taxa = 1", most of the dissimilarity is due to species turnover
qualitative_contrib_less_than_1 <- betalink_genotype_net_qualitative$contrib[which(betalink_genotype_net_qualitative$contrib < 1)]
hist(qualitative_contrib_less_than_1)
mean(qualitative_contrib_less_than_1) # 48% of the dissimilarity.

# Betalink Quantitative
genotype_net_matrices <- list()
for(i in 1:length(genotype_net)) {
  genotype_net_matrices[[i]] <- as.matrix(genotype_net[[i]])
} # need to be converted into matrices for "as.table" function to work for quantitative data.

betalink_genotype_net_quantitative <- betalink.dist(genotype_net_matrices, bf = "euclidean", triangular = T) #

mean(betalink_genotype_net_quantitative$WN) # 67% dissimilarity among genotypes with "horn"
hist(betalink_genotype_net_quantitative$WN)

mean(betalink_genotype_net_quantitative$OS)
hist(betalink_genotype_net_quantitative$OS)

mean(betalink_genotype_net_quantitative$S_all.taxa)
mean(betalink_genotype_net_quantitative$S_upper)
mean(betalink_genotype_net_quantitative$S_lower)
mean(betalink_genotype_net_quantitative$ST)

mean(betalink_genotype_net_quantitative$contrib) # species turnover contributes to 66% of the dissimilarity in interaction networks, according to euclidean distance.
hist(betalink_genotype_net_quantitative$contrib)

########### Everything below this may need to be adjusted if I use it.
# even after removing combinations with absolutely no species in common "S_all.taxa = 1", most of the dissimilarity is due to species turnover
table(betalink_genotype_net_quantitative$contrib < 1)
quantitative_contrib_less_than_1 <- betalink_genotype_net_quantitative$contrib[which(betalink_genotype_net_quantitative$contrib < 1)]
hist(quantitative_contrib_less_than_1)
mean(quantitative_contrib_less_than_1) # species turnover only contributes to 44% of the dissimilarity



##### Let's examine how modularity corresponds with the dissimilarity of species interaction networks.

# upload module data
module.id.info <- read.csv("~/Documents/Genotype_Networks/data/module_info_genotype_gall_parasitoid_network.csv")
genotype.module.info <- module.id.info %>%
  filter(Trophic == "Genotype") %>%
  mutate(Module.ID = factor(Module.ID))

# create hierarchical clustering to compare modularity to
clust.betalink.q <- hclust(betalink_genotype_net_quantitative$WN, method = "ward")
plot(clust.betalink.q, labels = names(genotype_net))
cluster.info <- factor(cutree(clust.betalink.q, k = 3))

cluster.df <- data.frame(Genotypes = names(genotype_net), ClusterID = cluster.info)

comp.betalink.q <- hclust(betalink_genotype_net_quantitative$WN, method = "complete")
plot(comp.betalink.q, labels = names(genotype_net))

# test for differences among modules in some network properties
library(bipartite)

networklevel_geno_net <- function(x) networklevel(x, intereven = "sum", logbase = 2) # changes from default
geno_net_indices <- ldply(.data = genotype_net_matrices, .fun = networklevel_geno_net) # this function is sweet.
geno_net_indices <- mutate(geno_net_indices, Genotypes = names(genotype_net))
geno_net_indices <- left_join(geno_net_indices, select(genotype.module.info, Genotypes = Vertex, Module.ID), by = "Genotypes")

geno_net_indices_cluster <- left_join(geno_net_indices, cluster.df, by = "Genotypes")

geno_net_indices_cluster_vLG_ptism <- left_join(geno_net_indices_cluster, select(vLG_geno_parasitism, Genotypes = Genotype, vLG_ptoid_attack, vLG_vLG.pupa, vLG_prop_ptoid_attack), by = "Genotypes")

# test for an effect of Module ID on weighted linkage density (or complexity)
plot(geno_net_indices_cluster_vLG_ptism$"linkage density" ~ Module.ID, geno_net_indices_cluster_vLG_ptism)

wLD_lm <- lm(geno_net_indices_cluster_vLG_ptism$"linkage density" ~ Module.ID, geno_net_indices_cluster_vLG_ptism) # qualitatively the same results even if Module #2 (row #24) is removed from the analysis.
summary(wLD_lm)
anova(wLD_lm)

summary(lm(geno_net_indices_cluster_vLG_ptism$"linkage density" ~ ClusterID, geno_net_indices_cluster_vLG_ptism)) # same pattern revealed even when I use an independent clustering algorithm

# test for an effect of Module ID on attack rates of vLG. Only used those modules for which vLG participated in.
plot(vLG_prop_ptoid_attack ~ Module.ID, filter(geno_net_indices_cluster_vLG_ptism, Module.ID %in% c(1,3,5)))

vLG_ptism_module <- glm(cbind(vLG_ptoid_attack, vLG_vLG.pupa) ~ Module.ID, filter(geno_net_indices_cluster_vLG_ptism, Module.ID %in% c(1,3,5)), family = "binomial") # still significant even if data point #7 isn't removed.
summary(vLG_ptism_module)
anova(vLG_ptism_module, test = "Chi")
plot(vLG_ptism_module) # 7 is an outlier. Results qualitatively the same when #7 is removed

plot(vLG_prop_ptoid_attack ~ ClusterID, geno_net_indices_cluster_vLG_ptism) 
summary(glm(cbind(vLG_ptoid_attack, vLG_vLG.pupa) ~ ClusterID, geno_net_indices_cluster_vLG_ptism, family = "binomial")) # much stronger effect when I use the cluster assignments.


# use PERMANOVA to test whether module ID explains the dissimilarity in gall-parasitoid interaction networks among genotypes
ST_mod <- betalink_genotype_net_quantitative$ST
ST_mod[c(22,292)] <- c(0,0)

contrib_mod <- betalink_genotype_net_quantitative$contrib
contrib_mod[c(22,292)] <- c(0,0)

adonis(betalink_genotype_net_quantitative$WN ~ genotype.module.info$Module.ID, permutations = 1000) # R2 = 66%, P < 0.001
adonis(betalink_genotype_net_quantitative$ST ~ genotype.module.info$Module.ID, permutations = 1000) # won't calculate because 2 dissimilarities are negative...
adonis(ST_mod ~ genotype.module.info$Module.ID, permutations = 1000)
adonis(betalink_genotype_net_quantitative$S_upper ~ genotype.module.info$Module.ID, permutations = 1000)
adonis(betalink_genotype_net_quantitative$S_lower ~ genotype.module.info$Module.ID, permutations = 1000)
adonis(betalink_genotype_net_quantitative$OS ~ genotype.module.info$Module.ID, permutations = 1000)
adonis(contrib_mod ~ genotype.module.info$Module.ID, permutations = 1000) # no differences in the contribution of network dissimilarity, suggesting that a similar process is occuring both between and among modules.

adonis(betalink_genotype_net_quantitative$WN ~ cluster.info, permutations = 1000) # explains 55% with only 3 cluster (74% with 5 cluster though), so it beats modularity... but this may be an unfair comparison since the clustering is based directly off these dissimilarities.

module.betalink.betadisper <- betadisper(betalink_genotype_net_quantitative$WN, genotype.module.info$Module.ID, bias.adjust = TRUE) # did the bias adjustment for small sample sizes.
boxplot(module.betalink.betadisper)
anova(module.betalink.betadisper) # I guess the variance is different between groups... May be driven by Module 3.

module.betalinkOS.betadisper <- betadisper(betalink_genotype_net_quantitative$OS, genotype.module.info$Module.ID, bias.adjust = TRUE) # did the bias adjustment for small sample sizes.
boxplot(module.betalinkOS.betadisper)
anova(module.betalinkOS.betadisper)

module.betalinkST.betadisper <- betadisper(ST_mod, genotype.module.info$Module.ID, bias.adjust = TRUE) # did the bias adjustment for small sample sizes.
boxplot(module.betalinkST.betadisper)
anova(module.betalinkST.betadisper)

# identify "how" different modules differ in their interaction networks
betalink.meandist <- meandist(betalink_genotype_net_quantitative$WN, module.id.info$Module.ID)
low <- lower.tri(betalink.meandist)
betalink.meandist.lower <- betalink.meandist[which(low == TRUE)]
mean(betalink.meandist.lower) # 80% mean dissimilarity among modules (horn)
sd(betalink.meandist.lower)
mean(diag(betalink.meandist), na.rm = T) # 44% dissimilarity within modules
sd(diag(betalink.meandist), na.rm = T)

mean(betalink.meandist.lower)/mean(diag(betalink.meandist), na.rm = T)

contrib.meandist <- meandist(contrib_mod, module.id.info$Module.ID)
low.contrib <- lower.tri(contrib.meandist)
contrib.meandist.lower <- contrib.meandist[which(low.contrib == TRUE)]
mean(contrib.meandist.lower) # 82% 
mean(diag(contrib.meandist), na.rm = T) 

# parasitoids (LOWER)
ptoid.meandist <- meandist(betalink_genotype_net_quantitative$S_lower, module.id.info$Module.ID)
low.ptoid <- lower.tri(ptoid.meandist)
ptoid.meandist.lower <- ptoid.meandist[which(low.ptoid == TRUE)]
mean(ptoid.meandist.lower[-5]) # 0.58 to 0.63
mean(diag(ptoid.meandist), na.rm = T) # 0.41

# galls (UPPER)
gall.meandist <- meandist(betalink_genotype_net_quantitative$S_upper, module.id.info$Module.ID)
low.gall <- lower.tri(gall.meandist)
gall.meandist.lower <- gall.meandist[which(low.gall == TRUE)]
mean(gall.meandist.lower) # 0.52
mean(diag(gall.meandist), na.rm = T) # 0.2

plot(ptoid.meandist.lower ~ gall.meandist.lower) # parasitoid dissimilarity is usually high among compartments, whereas gall dissimilarity varies from low to high.

# plot dissimilarity among genotypes with module ID overlain
library(ggplot2)
library(ggvegan)
source('~/Documents/miscellaneous_R/ggplot_themes.R')

# testing different ones
metaMDS_qWN <- metaMDS(betalink_genotype_net_quantitative$WN) # stress = 0.1168
#fortify_qWN <- fortify(metaMDS_qWN)

metaMDS_qWN_hack <- data.frame(metaMDS_qWN$points)
metaMDS_qWN_hack <- mutate(metaMDS_qWN_hack, Genotype = genotype.module.info$Vertex, Module.ID = genotype.module.info$Module.ID) # cluster.info

# color schemes
col.mod1 <- "red3"
col.mod2 <- "grey"
col.mod3 <- "steelblue"
col.mod4 <- "gold"
col.mod5 <- "darkorchid4"

net_dissimWN <- ggplot(metaMDS_qWN_hack, aes(x = MDS1, y = MDS2, shape = Module.ID, fill = Module.ID)) + 
  geom_point(size = 20) + 
  geom_text(data = metaMDS_qWN_hack, aes(x = MDS1, y = MDS2, label = Genotype),
            size = 7) +
  coord_fixed() + 
  scale_shape_manual(values = c(21,23,22,25,24), guide = "none") + 
  scale_fill_manual(values = c(col.mod1, col.mod2, col.mod3, col.mod4, col.mod5), guide = "none") + 
  xlab("NMDS 1") + 
  ylab("NMDS 2") +
  theme_ordination_black_white + 
  theme(axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text = element_text(size = 20)) 

#Code to override clipping
net_dissimWN
gt <- ggplot_gtable(ggplot_build(net_dissimWN))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

#summary(capscale(betalink_genotype_net_quantitative$WN ~ genotype.module.info$Module.ID))


# use betalink plot to visualize differences in interaction networks among genotypes
#library(bipartite)
#plotweb(genotype_net_matrices[[19]])
#plotweb(genotype_net_matrices[[25]])
source('~/Documents/betalink/R/betalink.plot.r')
source('~/Documents/betalink/R/metaweb.r')

S_qual <- genotype_net_matrices[[19]]
S_qual[S_qual > 0] <- 1

Xqual <- genotype_net_matrices[[23]]
Xqual[Xqual > 0] <- 1

Tqual <- genotype_net_matrices[[20]]
Tqual[Tqual > 0] <- 1

betalink.plot(S_qual, Xqual) # comparing genotypes S & X
betalink.q(genotype_net_matrices[[19]], genotype_net_matrices[[23]], bf = "bray") # 0.58 dissimilarity between genotypes, species turnover contributes to 2.5% of this

betalink.plot(S_qual, Tqual)
betalink.q(genotype_net_matrices[[19]], genotype_net_matrices[[20]], bf = "bray") # horn doesn't count these as dissimilar...Although they are different....What gives? 

betalink.plot(Xqual, Tqual)


#########  Mechanisms leading to dissimilarity in gall-parasitoid interaction networks
galldens <- read.csv("~/Documents/Genotype_Networks/data/gall_density_genotype_summary_df.csv")
galldens <- select(galldens, -X)

#gallsize <- read.csv("~/Documents/Genotype_Networks/data/gall_size_genotype_summary_df.csv")
#gallsize <- select(gallsize, -X)

vLG_gallsize <- read.csv("~/Documents/Genotype_Networks/data/vLG_gall_size_geno.df.csv")
vLG_gallsize <- select(vLG_gallsize, -X)

galldens.size <- left_join(galldens, vLG_gallsize)
galldens.size.noU <- filter(galldens.size, Genotype != "U")
galldens.size.noU <- mutate(galldens.size.noU, prop.small.galls = small.vLG.count/(large.vLG.count + small.vLG.count))

# vLG_Platy
plot(vLG_Platy ~ vLG_abund, galldens.size.noU)

#library(car)
#scatterplotMatrix(galldens.size.noU[ ,2:6])
#plot(log(rG_abund) ~ log(vLG_abund), data = galldens.size.noU)
#cor.test(log(galldens.size.noU$rG_abund+1), log(galldens.size.noU$vLG_abund+1))

cap_qWN <- capscale(betalink_genotype_net_quantitative$WN ~ galldens.size.noU$vLG_abund*galldens.size.noU$g.height.vLG.mean, na.action = na.omit)

summary(cap_qWN)
plot(cap_qWN)
#text(cap_qWN, labels = names(genotype_net))
anova(cap_qWN, by = "margin", step = 10000) # marginally significant effect of vLG_size

cap_qWN <- capscale(betalink_genotype_net_quantitative$WN ~ genotype.module.info$Module.ID)
plot(cap_qWN, type = "n")
points(cap_qWN, col = genotype.module.info$Module.ID)
text(cap_qWN, labels = genotype.module.info$Vertex)

metaMDS_qWN <- metaMDS(betalink_genotype_net_quantitative$WN)
plot(metaMDS_qWN, type = "n")
#text(metaMDS_qWN, labels = genotype.module.info$Vertex)
points(metaMDS_qWN, col = genotype.module.info$Module.ID)

cb_palette <- c("#56B4E9", "#999999", "#E69F00", "#009E73", "#CC79A7")
fortify_qWN <- fortify(metaMDS_qWN)

metaMDS_qWN_hack <- data.frame(metaMDS_qWN$points)
metaMDS_qWN_hack <- mutate(metaMDS_qWN_hack, Genotype = genotype.module.info$Vertex, Module.ID = genotype.module.info$Module.ID)
levels(metaMDS_qWN_hack$Genotype)[1] <- "C" # for plotting aesthestics.

source('~/Documents/miscellaneous_R/ggplot_themes.R')
source('~/Documents/ggvegan/R/fortify.metaMDS.R')
library(ggvegan)

ggplot(metaMDS_qWN_hack, aes(x = MDS1, y = MDS2, shape = Module.ID, fill = Module.ID)) +
  geom_point(size = 20) + 
  coord_fixed() + 
  scale_shape_manual(values = 21:25, guide = "none") + 
  scale_fill_manual(values = c("blue","grey","red","yellow","white"), guide = "none") + 
  theme_ordination_black_white + 
  xlab("NMDS 1") + 
  ylab("NMDS 2") +
  geom_text(data = metaMDS_qWN_hack, aes(x = MDS1, y = MDS2, label = Genotype))
ggsave("~/Documents/Genotype_Networks/gall_net_dissimilarity_NMDS_ordination.png", height = 9, width = 9, units = "in", dpi = 300)

fortify_qWN_genotypes <- filter(fortify_qWN, Score == "sites")
fortify_qWN_genotypes <- mutate(fortify_qWN_genotypes, Genotype = genotype.module.info$Vertex, Module.ID = genotype.module.info$Module.ID)

ggplot(fortify_qWN_genotypes, aes(x = Dim1, y = Dim2, shape = Module.ID, fill = Module.ID)) + geom_point(size = 10) + coord_fixed() + scale_shape_manual(values = 21:25, guide = "none") + scale_fill_manual(values = cb_palette, guide = "none") + theme_ordination_black_white


autoplot(cap_qWN, geom = "point", layers = "sites", size = 12)

adonis(betalink_genotype_net_quantitative$WN ~ genotype.module.info$Module.ID, permutations = 10000)
adonis(betalink_genotype_net_quantitative$WN ~ cluster.info) # interesting, clustering provides a better partition of the community, but this may be an unfair comparison since the clustering is based directly off these dissimilarities.
module.betalink.betadisper <- betadisper(betalink_genotype_net_quantitative$WN, genotype.module.info$Module.ID, bias.adjust = TRUE) # did the bias adjustment for small sample sizes.
boxplot(module.betalink.betadisper)
anova(module.betalink.betadisper)

test.gall.abunds <- read.csv('~/Documents/Genotype_Networks/data/gall_density_genotype_summary_df.csv')
test.gall.abunds <- select(test.gall.abunds, aSG_abund:rsLG_abund)
test.gall.abunds.dist <- vegdist(test.gall.abunds[-21, ], method = "bray")

mantel(test.gall.abunds.dist, betalink_genotype_net_quantitative$WN, method = "spearman", permutations = 10000) # results are robust to the method. Need to try and incorporate gall size now...


# create

