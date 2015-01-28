# NOTE THAT I HAVE NOT PERFORMED THESE ANALYSES AFTER REMOVING PP 373 (GENOTYPE S: OUTLYING DATAPOINT FOR VLG AND RG)

## upload source code for managing network data at the tree-level
source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')

## upload source code for analyzing beta-diversity of interaction networks. Note that this code is necessary for partitioning the components of network beta-diversity, but unnecessary for calculating the total dissimilarity of interaction networks (species turnover + interaction turnover)
source('~/Documents/betalink/R/betalink.R')
source('~/Documents/betalink/R/vec2data.frame.R')
source('~/Documents/betalink/R/measures.r')
source('~/Documents/betalink/R/betalink.b.r')
source('~/Documents/betalink/R/betalink.dist.r')
source('~/Documents/betalink/R/betalink.q.R')

## upload required libraries for analysis
library(vegan)

## main dataset is 'tree_level_interaxn_all_plants'

#### Analysis of dissimilarity of gall-parasitoid interaction networks among genotypes
abund_interaxns_gallsurv_noPont <- rowSums(tree_level_interaxn_all_plants_traits_size[ ,interaxns_gallsurv_noPont] 

abund_interaxns_noPont <- rowSums(tree_level_interaxn_all_plants_traits_size[ ,interaxns_noPont])

table(filter(tree_level_interaxn_all_plants_traits_size, abund_interaxns_noPont > 0)$Genotype)
## when 'net_df <- filter(tree_level_interaxn_all_plants_traits_size, abund_interaxns_noPont > 0'. The replicate samples of many of the Genotypes is reduced to a very low level. To see whether removing these genotypes affects the outcome of the analysis, I have created the following subsets. 
geno_2plus <- c("*","A","B","D","E","F","I","K","L","Q","S","T","V","W","X","Y","Z") # 17
geno_3plus <- c("*","B","D","I","K","L","Q","S","V","W","X","Y","Z") # 13
geno_4plus <- c("*","B","I","K","L","S","V","X","Z") # 9 (if Pontania galls were included, then I would have retained Genotype D as well.)

## dataset to use for PERMANOVA
net_df <- tree_level_interaxn_all_plants_traits_size#filter(tree_level_interaxn_all_plants_traits_size, abund_interaxns_noPont > 0, Genotype %in% geno_3plus)

## Permutational multivariate analysis of variance (PERMANOVA)
net_qual <- ifelse(net_df[ ,interaxns_noPont] > 0, 1, 0)

guild_adonis <- adonis(sqrt(net_df[ ,interaxns_guild]/shootEst.no18) ~ Genotype,
                        data = net_df,
                        #strata = net_df$Genotype,
                        method = "euclidean")
sum(colSums(net_df[ ,interaxns_guild])[c(2,4,5)])/sum(colSums(net_df[ ,interaxns_guild]))
plot(rda(sqrt(net_df[ ,interaxns_guild]/net_df$shootEst.no18) ~ Genotype, data = net_df), display = c("sp"), scaling = 2) 

net_df_adonis <- adonis(sqrt(net_df[ ,interaxns_noPont]/shootEst.no18) ~ Genotype,
       data = net_df,
       #strata = net_df$Genotype,
       method = "euclidean")


C.replace <- levels(net_df$Genotype)
C.replace[1] <- "C"
net_df$Genotype <- factor(net_df$Genotype, levels = levels(net_df$Genotype), labels = C.replace)
net_df_meandist <- meandist(vegdist(sqrt(net_df[ ,interaxns_noPont]/net_df$shootEst.no18), 
                                    method = "euclidean"),
                            grouping = net_df$Genotype)
table.geno <- table(net_df$Genotype)
names(table.geno)[1] <- "C" # switch to C for plotting aesthetics
net_df_meandist <- as.dist(net_df_meandist)
net_df_hclust <- hclust(d = net_df_meandist, method = "ave", members = table.geno)
plot(net_df_hclust, hang = -1)

library(ggplot2)
library(ggdendro)
net_df_dend <- as.dendrogram(net_df_hclust)
ddata <- dendro_data(net_df_dend, type = "rectangle")

ggdendrogram(net_df_hclust, leaf_labels = FALSE, theme_dendro = FALSE) + 
  scale_y_reverse() + 
  geom_text(data = label(ddata), aes(x = x, y = y-0.003, label = label)) + 
  ylab("Euclidean distance") +
  #theme_dendro + 
  theme(axis.text.y = element_text(size = 15, angle = 0, color = "black"),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), axis.title.y = element_text(size = 18, vjust = 1),
        axis.line.y = element_line(color = "black"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())


plot(rda(sqrt(net_df[ ,ptoids.noPont]/net_df$shootEst.no18) ~ Genotype, data = net_df), display = c("sp"), scaling = 2) # look at relationships among species. In line with what I expected. All of the ectoparasitoids are aligned on the same axis.

adonis(net_qual ~ Genotype, data = net_df, method = "euclidean") # results are qualitatively the same whether or not I included shoot estimates as a covariate in the model

rc_dis <- raupcrick(net_qual) # calculate raupcrick dissimilarity for teasing apart differences in alpha vs. beta-diversity
adonis(rc_dis ~ Genotype, data = net_df)

anova(betadisper(d = vegdist(sqrt(net_df[ ,ptoids.noPont]/net_df$shootEst.no18), method = "euclidean"),
           group = net_df$Genotype, bias.adjust = TRUE)) # significant overdispersion. However, when I use subsets of the same dataset, Genotype still remains a significant predictor of network dissimilarity and I no longer have overdispersion.
anova(betadisper(d = vegdist(net_qual, method = "euclidean"), group = net_df$Genotype, bias.adjust = TRUE)) # still significant overdispersion with the qualitative data and the full dataset.

## Code for calculating heritability of network dissimilarity. Triple check this code.
# First, calculate the average sample size among all the genotypes.
# Then calculate the mean Sum of Squares due to Genotype and the Residual
# Subtract Genotype MS - Residual MS, and divide that by the sample size to get the variance due to Genotype.
# Divide Genotypic variance by the Genotype variance plus residual mean square (which is equal to the residual variance) to obtain heritability estimate.
gall.adonis.ss <- mean(table(net_df$Genotype)) # 5.6 average sample size for genotypes
gall.gen.ms <- net_df_adonis$aov.tab$MeanSqs[1]
gall.res.ms <- net_df_adonis$aov.tab$MeanSqs[2] # also equal to residual variance
gall.gen.var <- (gall.gen.ms - gall.res.ms)/gall.adonis.ss

gall.h2 <- gall.gen.var/(gall.gen.var + gall.res.ms)
gall.h2

guild.adonis.ss <- mean(table(net_df$Genotype)) # 5.6 average sample size for genotypes
guild.gen.ms <- guild_adonis$aov.tab$MeanSqs[1]
guild.res.ms <- guild_adonis$aov.tab$MeanSqs[2] # also equal to residual variance
guild.gen.var <- (guild.gen.ms - guild.res.ms)/guild.adonis.ss

guild.h2 <- guild.gen.var/(guild.gen.var + guild.res.ms)
guild.h2
#### Partition beta-diversity of interaction networks into species turnover and interaction components.

tree_net_filter <- gall_net_melt %>%
  filter(gall_contents %in% c("Eulo.fem", "Eulo.mal", "Lestodip", "Mesopol", "Mymarid", "Platy", "Ptero.2", "Tory.fem", "Tory.mal")) # only retaining gall-parasitoid interactions. Removed rsLG parasitoids "Eury.fem", "Eury.mal", "Lathro.fem", "Lathro.mal",

tree_net_add <- mutate(tree_net_filter, gall_contents_collapse = revalue(gall_contents, c("Eulo.fem" = "Eulo", "Eulo.mal" = "Eulo",  "Ptero.2" = "Mesopol", "Tory.fem" = "Tory", "Tory.mal" = "Tory") )) # removed rsLG parasitoids: "Eury.fem" = "Eury", "Eury.mal" = "Eury", "Lathro.fem" = "Lathro", "Lathro.mal" = "Lathro",

## this only retains plants with known gall-parasitoid interactions (77 plants). 
tree_net <- cast(tree_net_add, gall.sp ~ gall_contents_collapse | plant.position, sum) 

## this loop prepares each network for analysis in betalink. 
for(i in 1:length(tree_net)) {
  rownames(tree_net[[i]]) <- as.character(tree_net[[i]]$gall.sp) # change row names to gall specie names
  tree_net[[i]] <- as.data.frame(select(tree_net[[i]], -gall.sp)) # remove gall.sp as a column of data. Turning this into a matrix is important for easier analysis with betalink package
}

## beta-diversity partitioning of qualitative data
tree_net_qual <- list()
for(i in 1:length(tree_net)) {
  tree_net_qual[[i]] <- ifelse(tree_net[[i]] > 0, 1, 0)
} 

betalink_tree_net_qualitative <- betalink.dist(tree_net_qual, bf = "euclidean", triangular = T) 

mean(betalink_tree_net_qualitative$WN) 
mean(betalink_tree_net_qualitative$OS) 
mean(betalink_tree_net_qualitative$S_all.taxa)
mean(betalink_tree_net_qualitative$S_upper) 
mean(betalink_tree_net_qualitative$S_lower)
mean(betalink_tree_net_qualitative$ST) 
mean(betalink_tree_net_qualitative$contrib) 
hist(betalink_tree_net_qualitative$contrib)
#genotype_net <- cast(genotype_net_add, gall.sp ~ gall_contents_collapse | Genotype, sum)

## Partition beta-diversity using quantitative data 
tree_net_trans <- list()
for(i in 1:length(tree_net)) {
  tree_net_trans[[i]] <- log(tree_net[[i]]+1)
} 

betalink_tree_net_quantitative <- betalink.dist(tree_net_trans, bf = "euclidean", triangular = T) #

mean(betalink_tree_net_quantitative$WN) 
hist(betalink_tree_net_quantitative$WN)

mean(betalink_tree_net_quantitative$OS)
hist(betalink_tree_net_quantitative$OS)

mean(betalink_tree_net_quantitative$S_all.taxa)
mean(betalink_tree_net_quantitative$S_upper)
mean(betalink_tree_net_quantitative$S_lower)
mean(betalink_tree_net_quantitative$ST)
hist(betalink_tree_net_quantitative$ST)

mean(betalink_tree_net_quantitative$contrib) 
hist(betalink_tree_net_quantitative$contrib)

#### PERMANOVA

## align tree info with tree numbering for 'tree_net' for PERMANOVA
pp.df <- data.frame(plant.position = as.character(names(tree_net)))
tree_info_df_dis <- join(x = pp.df, 
                         y = select(tree_level_interaxn_all_plants, plant.position, 
                                    shootEst.no18, Gender, Genotype), 
                         by = "plant.position")

## Qualitative. Dissimilarity in networks is primarily due to species turnover.
adonis(betalink_tree_net_qualitative$WN ~ log(shootEst.no18) + Genotype, tree_info_df_dis) # WN
adonis(betalink_tree_net_qualitative$OS ~ log(shootEst.no18) + Genotype, tree_info_df_dis) # OS
adonis(betalink_tree_net_qualitative$ST ~ log(shootEst.no18) + Genotype, tree_info_df_dis) # ST

## Quantitative.
adonis(betalink_tree_net_quantitative$WN ~ log(shootEst.no18) + Genotype, tree_info_df_dis) # WN
adonis(betalink_tree_net_quantitative$OS ~ log(shootEst.no18) + Genotype, tree_info_df_dis) # OS
adonis(betalink_tree_net_quantitative$ST ~ log(shootEst.no18) + Genotype, tree_info_df_dis) # ST

#### Visualizing differences in network structure among genotypes. 
# Play around with a lot of the variations. This should also help guide which dissimilarity metric will be most appropriate.
# Qualitative or quantitative? Mean value of each interaction, corrected for the number of shoots sampled? May first be a good idea to do an ordination of the data. 
# For visualizing the total dissimilarity, I'm going to chose the ones setup in the "community" format so I can visualize which interactions are driving the differences in the community.

# right now, the raup crick visualization is looking pretty reasonable
# When plotting the full dataset (N = 146), I'm getting this: Error in chol.default(cov) : the leading minor of order 2 is not positive definite. May explain why some of the genotype standard errors aren't plotting.

net_df_plot <- filter(tree_level_interaxn_all_plants_traits_size, abund_interaxns_noPont > 0, Genotype %in% geno_3plus) # need at least 3 data points for each genotype to plot ellipses
net_qual_plot <- ifelse(net_df_plot[ ,interaxns_noPont] > 0, 1, 0)
rc_dis_plot <- raupcrick(net_df_plot[ ,interaxns_noPont])

rda_sub <- rda(sqrt(net_df_plot[ ,interaxns_noPont]/net_df_plot$shootEst.no18) ~ Genotype, 
               data = net_df_plot)
RsquareAdj(rda_sub)
anova(betadisper(d = vegdist(sqrt(net_df_plot[ ,interaxns_noPont]/net_df_plot$shootEst.no18), 
                       method = "euclidean"),
                 group = net_df_plot$Genotype, bias.adjust = TRUE))
plot(rda_sub, display = "sp") # species correlations patterns aren't the same as with the full dataset.
plot(rda_sub, display = c("cn","sp"), scaling = 3)
plot(rda_sub, display = "cn")
ordiellipse(rda_sub, net_df_plot$Genotype, kind = "se")

## RDA of presence/absence network
rda_qual <- rda(net_qual ~ Genotype, 
                data = net_df_plot)
RsquareAdj(rda_qual)

plot(rda_qual, display = "sp")
plot(rda_qual, display = c("cn","sp"), scaling = 3)
plot(rda_qual, display = "cn")
ordiellipse(rda_qual, net_df_plot$Genotype, kind = "se") 

## Capscale of Raup-crick dissimilarity
cap_rc <- capscale(rc_dis ~ Genotype, 
                   data = net_df_plot)
RsquareAdj(cap_rc)
plot(cap_rc, display = c("cn"))
ordiellipse(cap_rc, net_df_plot$Genotype, kind = "se") # suggests that Genotype Q is driving much of the beta-diversity in interactions among the subsetted dataset.

## Capscale for dissimilarity using relative species composition (e.g. bray curtis or horn)
cap_sub <- capscale(sqrt(net_df_plot[ ,interaxns_noPont]/net_df_plot$shootEst.no18) ~ Genotype, 
                   data = net_df_plot, distance = "bray")
RsquareAdj(cap_sub)
plot(cap_sub, display = "sp")
plot(cap_sub, display = "cn")
ordiellipse(cap_sub, net_df_plot$Genotype, kind = "se")

#### Plot for manuscript
ord <- rda_sub

## extract centroids for plotting
centroids <- scores(ord, choices = c(1,2), display = "cn")
row.names(centroids) <- c("C","B","D","I","K","L","Q","S","V","W","X","Y","Z") # replaced "*" with "C" for plotting aesthetics.

## generate plot. Coded out the manner in which I saved the plot

pdf(file="~/Documents/Genotype_Networks/gall qualitative network 13 genotypes.pdf",width=8, height=8)
plot.new()
plot.window( xlim=c(-0.5,0.5), ylim = c(-0.25,0.25), asp = 1) 
#abline(h=0, lty="dotted")
#abline(v=0, lty="dotted")
axis(side=1, tck = -0.015, padj=-0.5)
axis(side=2, tck= -0.015, las=1)
title(xlab="RDA 1", ylab="RDA 2") 
box()

## fill in the plot
ordiellipse(ord, net_df_plot$Genotype, kind = "se", draw = "lines", col = "grey")
#ordiellipse(ord, net_df_plot$Genotype, kind = "se", show.groups = "K", col = "grey", draw = "polygon")
ordiellipse(ord, net_df_plot$Genotype, kind = "se", show.groups = "X", col = "grey", draw = "polygon")
ordiellipse(ord, net_df_plot$Genotype, kind = "se", show.groups = "Q", col = "grey", draw = "polygon")
#ordiellipse(ord, net_df_plot$Genotype, kind = "se", show.groups = "W", col = "grey", draw = "polygon")
ordiellipse(ord, net_df_plot$Genotype, kind = "se", show.groups = "S", col = "grey", draw = "polygon")
text(x = centroids[ ,1], centroids[ ,2], labels = row.names(centroids))
dev.off()

## Create subplots of qualitative networks
# source in functions
source('~/Documents/ggnet/bipartite_plot_info.R')
source('~/Documents/ggnet/ggnet_bipartite.R')

library(ggplot2)
library(grid)

## this only retains plants with known gall-parasitoid interactions (77 plants). 
tree_net_vis <- cast(tree_net_add, gall.sp ~ gall_contents_collapse | Genotype, sum) 
tree_net_vis <- tree_net_vis[-21] # remove genotype U because it had no interactions
for(i in 1:length(tree_net_vis)) {
  rownames(tree_net_vis[[i]]) <- as.character(tree_net_vis[[i]]$gall.sp) # change row names to gall specie names
  tree_net_vis[[i]] <- as.data.frame(select(tree_net_vis[[i]], -gall.sp)) # remove gall.sp as a column of data. Turning this into a matrix is important for easier analysis with betalink package
}

tree_net_vis$Q <- tree_net_vis$Q[c("rG", "aSG"), ] # wanted rG to come first because it is easier to see turnover in gall species.

S_info <- bipartite_plot_info(web = tree_net_vis$S, order.type = "normal")
X_info <- bipartite_plot_info(web = tree_net_vis$X, order.type = "normal")
Q_info <- bipartite_plot_info(web = tree_net_vis$Q, order.type = "normal")
arrange(S_info[[2]], vertex.names)
arrange(X_info[[2]], vertex.names)
arrange(Q_info[[2]], vertex.names)

## Shape assignments
# parasitoids: Eulo = 3, Eury = 4, Mesopol = 7, Platy = 8, Tory = 9, Lestodip = 13
# galls: rG = 21, rsLG = 22, vLG = 23, aSG = 24


## S plot
S_plot <- ggnet_bipartite(S_info) + 
  geom_point(data = S_info[[2]], aes(x = x, y = y, shape = vertex.names), size = 12, fill = "grey") +
  scale_shape_manual(values = c(3, 4, 7, 8, 21, 22, 9, 23), guide = "none") +
  theme(plot.margin = structure (c(2, 2, 2, 2), unit = "lines", valid.unit = 3L, class = "unit"))
S_plot
gtS <- ggplot_gtable(ggplot_build(S_plot))
gtS$layout$clip[gtS$layout$name == "panel"] <- "off"
grid.draw(gtS)
  
## X plot

max.seg.wt <- 22
# X: quantitative color
library(RColorBrewer)
X.df <- X_info[[1]] %>%
  filter(Sequence == 500 | Sequence == 1) %>%
  reshape(idvar = "Group", timevar = "Sequence", direction = "wide") %>%
  mutate(weight.trans = Weight.1) # unscaled weights
X.nodes <- X_info[[2]] %>%
  mutate(colrs = c(rab.bud.col, iteo.col, eulo.col, meso.col, platy.col, tory.col)) %>%
  arrange(vertex.names)

metaplot.col <- brewer.pal(11, "Spectral")
cecid.col <- metaplot.col[1]
eulo.col <- metaplot.col[2]
lesto.col <- metaplot.col[3]
meso.col <- metaplot.col[4]
mym.col <- metaplot.col[5]
platy.col <- metaplot.col[7]
rab.bud.col <- metaplot.col[8]
rab.stem.col <- metaplot.col[9]
tory.col <- metaplot.col[10]
iteo.col <- metaplot.col[11]

X.plot <- ggplot(X.df) + 
  geom_segment(data = X.df, aes(x = x.1, xend = x.500, y = y.1, yend = y.500),
               color = "grey",
               size = X.df$weight.trans/max.seg.wt*15,
               alpha = 0.75) +
  new_theme_empty + 
  geom_point(data = filter(X.nodes, y > 1), aes(x = x, y = y, fill = vertex.names), 
             shape = 25, size = 30) +
  geom_point(data = filter(X.nodes, y == 1), aes(x = x, y = y, fill = vertex.names), 
             shape = 21, size = 30) +
  scale_fill_manual(values = X.nodes$colrs, guide = "none")
X.plot
gtX <- ggplot_gtable(ggplot_build(X.plot))
gtX$layout$clip[gtX$layout$name == "panel"] <- "off"
grid.draw(gtX)
  #scale_fill_brewer(palette = "Spectral")

## F plot
# F: quantitative color
F_info <- bipartite_plot_info(web = tree_net_vis$F, order.type = "normal")
F.df <- F_info[[1]] %>%
  filter(Sequence == 500 | Sequence == 1) %>%
  reshape(idvar = "Group", timevar = "Sequence", direction = "wide") %>%
  mutate(weight.trans = Weight.1) # unscaled weights
F.nodes <- F_info[[2]] %>%
  mutate(colrs = c(cecid.col, iteo.col, eulo.col, meso.col, tory.col)) %>%
  arrange(vertex.names)

F.plot <- ggplot(F.df) + 
  geom_segment(data = F.df, aes(x = x.1, xend = x.500, y = y.1, yend = y.500),
               color = "grey",
               size = F.df$weight.trans/max.seg.wt*15,
               alpha = 0.75) +
  new_theme_empty + 
  geom_point(data = filter(F.nodes, y > 1), aes(x = x, y = y, fill = vertex.names), 
             shape = 25, size = 30) +
  geom_point(data = filter(F.nodes, y == 1), aes(x = x, y = y, fill = vertex.names), 
             shape = 21, size = 30) +
  scale_fill_manual(values = F.nodes$colrs, guide = "none")
F.plot
gtF <- ggplot_gtable(ggplot_build(F.plot))
gtF$layout$clip[gtF$layout$name == "panel"] <- "off"
grid.draw(gtF)

## P plot
# P: quantitative color
P_info <- bipartite_plot_info(web = tree_net_vis$P, order.type = "normal")
P.df <- P_info[[1]] %>%
  filter(Sequence == 500 | Sequence == 1) %>%
  reshape(idvar = "Group", timevar = "Sequence", direction = "wide") %>%
  mutate(weight.trans = Weight.1) # unscaled weights
P.nodes <- P_info[[2]] %>%
  mutate(colrs = c(iteo.col, platy.col)) %>%
  arrange(vertex.names)

P.plot <- ggplot(P.df) + 
  geom_segment(data = P.df, aes(x = x.1, xend = x.500, y = y.1, yend = y.500),
               color = "grey",
               size = P.df$weight.trans/max.seg.wt*15,
               alpha = 0.75) +
  new_theme_empty + 
  geom_point(data = filter(P.nodes, y > 1), aes(x = x, y = y, fill = vertex.names), 
             shape = 25, size = 30) +
  geom_point(data = filter(P.nodes, y == 1), aes(x = x, y = y, fill = vertex.names), 
             shape = 21, size = 30) +
  scale_fill_manual(values = P.nodes$colrs, guide = "none")
P.plot
gtP <- ggplot_gtable(ggplot_build(P.plot))
gtP$layout$clip[gtP$layout$name == "panel"] <- "off"
grid.draw(gtP)

# X: qualtiative shapes
X_plot <- ggnet_bipartite(X_info) + 
  geom_point(data = X_info[[2]], aes(x = x, y = y, shape = vertex.names), size = 12, fill = "grey") +
  scale_shape_manual(values = c(3, 7, 8, 21, 9, 23), guide = "none") +
  theme(plot.margin = structure (c(2, 2, 2, 2), unit = "lines", valid.unit = 3L, class = "unit"))
X_plot
gtX <- ggplot_gtable(ggplot_build(X_plot))
gtX$layout$clip[gtX$layout$name == "panel"] <- "off"
grid.draw(gtX)

## Q plot
Q_plot <- ggnet_bipartite(Q_info) +  
  geom_point(data = Q_info[[2]], aes(x = x, y = y, shape = vertex.names), size = 12, fill = "grey") +
  scale_shape_manual(values = c(24, 13, 21, 9), guide = "none") +
  theme(plot.margin = structure (c(2, 2, 2, 2), unit = "lines", valid.unit = 3L, class = "unit"))
Q_plot
gtQ <- ggplot_gtable(ggplot_build(Q_plot))
gtQ$layout$clip[gtQ$layout$name == "panel"] <- "off"
grid.draw(gtQ)

#### Look at dissimilarity in herbivore community and network dissimilarity as an explanation for network structure
library(vegclust)

vLG.height.total <- filter(tree_interaxns_filter_add, gall.sp == "vLG")$g.height
vLG.height.surv <- filter(tree_interaxns_filter_add, gall.sp == "vLG", gall_contents_collapse == "vLG.pupa")$g.height
vLG.height.para <- filter(tree_interaxns_filter_add, gall.sp == "vLG", gall_contents_collapse %in% c("Eulo","Mesopol","Mymarid","Platy", "Tory"))$g.height

## calculate selection differential for vLG size
mean(vLG.height.total, na.rm = TRUE) - mean(vLG.height.para, na.rm = TRUE)

## calculate intensity of selection for vLG size
(mean(vLG.height.total, na.rm = TRUE) - mean(vLG.height.para, na.rm = TRUE))/sd(vLG.height.total, na.rm = TRUE) # intensity of selection.

hist(filter(tree_interaxns_filter_add, gall.sp == "rsLG")$g.height)
boxplot.stats(filter(tree_interaxns_filter_add, gall.sp == "rsLG")$g.height)

hist(filter(tree_interaxns_filter_add, gall.sp == "rG")$g.height)
boxplot.stats(filter(tree_interaxns_filter_add, gall.sp == "rG")$g.height)

hist(filter(tree_interaxns_filter_add, gall.sp == "SG")$g.height)
boxplot.stats(filter(tree_interaxns_filter_add, gall.sp == "SG")$g.height)

hist(filter(tree_interaxns_filter_add, gall.sp == "aSG")$g.height)
boxplot.stats(filter(tree_interaxns_filter_add, gall.sp == "aSG")$g.height)

size.strata <- c(0, 4, 8, 12, 16)#quantile(tree_interaxns_filter_add$g.height, probs = c(0, 0.25, 0.5, 0.75), na.rm = TRUE) # only testing right now...
stratify.galls <- stratifyvegdata(tree_interaxns_filter_add,
                                  sizes = size.strata,
                                  plotColumn = "plant.position",
                                  sizeColumn = "g.height",
                                  speciesColumn = "gall.sp",
                                  counts = TRUE)
names(stratify.galls)
plot.CAP(stratify.galls, sizes = size.strata)

stratify.galls.prof <- vegdiststruct(stratify.galls, type = "profile", method = "bray")
stratify.galls.abund <- vegdiststruct(stratify.galls, type = "abundance", method = "bray")
stratify.galls.vol <- vegdiststruct(stratify.galls, type = "volume", method = "bray")

## align tree info with tree numbering for 'tree_net' for PERMANOVA
pp.df.stratify <- data.frame(plant.position = as.character(names(stratify.galls)))
tree_info_df_dis_stratify <- join(x = pp.df.stratify, 
                         y = select(tree_level_interaxn_all_plants, plant.position, 
                                    shootEst.no18, Gender, Genotype, aSG_Tory:rG_Platy, rG_Tory:rsLG_Lathro, SG_Platy, vLG_Eulo:vLG_Tory), 
                         by = "plant.position")
table(tree_info_df_dis_stratify$Genotype) # Genotype U only has one sample. And N and J have only 2 samples.

stratify.df <- stratify.galls.abund
cap.gall.stratify <- capscale(stratify.df ~ Condition(log(shootEst.no18)) + Genotype, tree_info_df_dis_stratify)
RsquareAdj(cap.gall.stratify)
plot(cap.gall.stratify, display = "cn")

#tree_net_for_gall.profile <- filter(tree_level_interaxn_all_plants, abund_ptoids_gallsurv > 0)

# make sure all of them match up.Interesting, not matching up at all right now. Maybe I need to be a bit more careful...
mantel(stratify.df, vegdist(log(tree_info_df_dis_stratify[ ,ptoids]+1), method = "bray"))

## START A NEW SCRIPT FOR LINKING TRAITS TO COMMUNITY RESPONSES AT THE TREE AND GENOTYPE LEVEL.
net_tree_sum <- tree_level_interaxn_all_plants[ ,ptoids]/tree_level_interaxn_all_plants$shootEst.no18*100

net_sum_geno <- aggregate(net_tree_sum, list(Genotype = tree_level_interaxn_all_plants$Genotype), mean)

net_sum_geno[c(17,19,24), ]

library(psych)
net_corrs <- corr.test(x = sqrt(net_sum_geno[ ,-1]))
round(net_corrs$p, 3)
scatterplotMatrix(net_sum_geno[ ,-1])

net_tree_corrs <- corr.test(x = net_tree_sum)
net_corrs

library(car)
scatterplotMatrix(net_tree_sum)

plot(y = sqrt(net_sum_geno$vLG_Platy), x = sqrt(net_sum_geno$vLG_Mesopol), type = "n")
text(y = sqrt(net_sum_geno$vLG_Platy), x = sqrt(net_sum_geno$vLG_Mesopol), labels = net_sum_geno$Genotype)
## dataset for Greg to play around with during PRIMER course
write.csv(select(tree_level_interaxn_all_plants, Gender, Genotype, plant.position, shootEst.no18, aSG_Tory:rG_Platy, rG_Tory:rsLG_Lathro, SG_Platy, vLG_Eulo:vLG_Tory), "~/Desktop/gall_network_data.csv")


