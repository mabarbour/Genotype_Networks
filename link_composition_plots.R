source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')
source('~/Documents/miscellaneous_R/ggplot_themes.R')
source('~/Documents/ggnet/bipartite_plot_info.R')

require(ggplot2)
require(gridExtra)
require(dplyr)
require(tidyr)

## create theme for figures
theme_links <- theme_bw() +
  theme(axis.title.x = element_text(size = 16, vjust = 0.1), 
        axis.title.y = element_text(size = 16, vjust = 1),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        panel.grid = element_blank(),
        legend.position = "none")
point.size.range <- c(1,8)
label.x.pos <- 1.5
label.y.galls <- 23
label.y.gallsize <- 14
label.y.ptoids <- 9.5
label.y.ptism <- 1
ABCD.allgalls <- data.frame(x = rep(label.x.pos,4), y = rep(label.y.galls,4), 
                            variable = c("Leaf gall", "apical-Stem gall", 
                                         "Bud gall", "mid-Stem gall"), 
                            labels = c("(A)","(B)","(C)","(D)"))
AB.domgalls <- data.frame(x = rep(label.x.pos,2), y = rep(label.y.galls,2), 
                          variable = c("Leaf gall","Bud gall"), 
                          labels = c("(A)","(B)"))
C.gallsize <- data.frame(x = label.x.pos, y = label.y.gallsize, labels = "(C)")
ABC.domptoids <- data.frame(x = rep(label.x.pos,3), y = rep(label.y.ptoids,3), 
                            Parasitoid = c("Platygaster", "Mesopolobus", "Torymus"), 
                            labels = c("(A)","(B)","(C)"))
D.ptism <- data.frame(x = label.x.pos, y = label.y.ptism, labels = "(D)")
A.linkabund <- data.frame(x = 0.5, y = 10, labels = "(A)")
B.ptism <- data.frame(x = 4.25, y = label.y.ptism, labels = "(B)")
#point.size <- 6
line.widths <- 3 # for link plots
## create color-blind friendly palette with grey (taken from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## change Genotype * to C for plotting aesthetics
levels(tree_level_interaxn_all_plants_traits_size$Genotype)[1] <- "C"

## tidy up gall composition data
gall.df <- as.data.frame(tree_level_interaxn_all_plants_traits_size) %>%
  select(Genotype, "Leaf gall" = vLG_abund, 
         "apical-Stem gall" = aSG_abund, 
         "Bud gall" = rG_abund, 
         "mid-Stem gall" = SG_abund) %>%
  gather(Genotype)

gall.summary <- gall.df %>% group_by(Genotype, variable) %>% summarise_each(funs(mean)) %>% arrange(variable)
max(filter(gall.summary, variable == "Leaf gall")$value) # 10
max(filter(gall.summary, variable == "Bud gall")$value) # 8
max(filter(gall.summary, variable == "apical-Stem gall")$value) # 1.4

## order genotypes by median leaf gall abundance for all plots with genotype on X axis. 
vLG.df.sum <- gall.df %>%
  filter(variable == "Leaf gall") %>%
  group_by(Genotype) %>%
  summarise(mean_abundance = mean(value)) %>%
  arrange(mean_abundance)
gall.df$Genotype <- factor(gall.df$Genotype, levels = vLG.df.sum$Genotype)

## create gall composition plot with all four gall species for supplementary info.
gall.community.plot <- ggplot(gall.df, aes(x = Genotype, y = value, fill = variable)) +
  geom_boxplot() + facet_wrap( ~ variable, ncol = 2) + 
  ylab("Gall abundance (#/branch)") +
  scale_fill_manual(values = cbPalette[c(7,8,3,1)]) +
  geom_text(data = ABCD.allgalls, aes(x = x, y = y, label = labels)) +
  theme_links

## create gall composition plot with the leaf and bud galls. The two species that are driving the differences in gall community composition. For Figure 2 in manuscript.
gall.dominants.plot <- ggplot(filter(gall.df, variable %in% c("Leaf gall", "Bud gall")), 
       aes(x = Genotype, y = value, fill = variable)) +
  geom_boxplot() + facet_wrap( ~ variable, ncol = 1) + 
  ylab("Gall abundance (#/branch)") + xlab("") +
  geom_text(data = AB.domgalls, aes(x = x, y = y, label = labels)) +
  scale_fill_manual(values = cbPalette[c(7,3)]) +
  theme_links

## tidy up leaf gall size data
vLG.size.df <- tree_level_interaxn_all_plants_traits_size %>%
  select(Genotype, vLG.height.mean, vLG.gall.count) %>%
  mutate(type = "mean")

lm.vLG.size <- lm(vLG.height.mean ~ Genotype, data = vLG.size.df, weights = vLG.gall.count)
summary(lm.vLG.size)
anova(lm.vLG.size)

## calculate weighted mean for leaf gall size based on the number of galls found on each replicate willow
vLG.size.df.2 <- vLG.size.df %>%
  group_by(Genotype) %>%
  summarise(vLG.height.mean = weighted.mean(vLG.height.mean, vLG.gall.count, 
                                              na.rm = TRUE),
            vLG.gall.count = 10) %>% # point size for weighted mean
  mutate(type = "weighted.mean")

## combine the original leaf gall size data with the weighted mean estimates.
vLG.size.df.3 <- rbind.data.frame(vLG.size.df.2,
                                  vLG.size.df)

## order leaf gall size data by median leaf gall abundance. This enables me to visually see whether there is a correlation between leaf gall abundance and size among genotypes.
vLG.size.df.3$Genotype <- factor(vLG.size.df.3$Genotype, levels = vLG.df.sum$Genotype)

## plot leaf gall size data. Size of points correspond to weights used to calculated the weighted mean of gall size for each genotype.
leaf.gall.size <- ggplot(vLG.size.df.3, 
       aes(x = Genotype, y = vLG.height.mean, size = vLG.gall.count, 
           fill = type, color = type, shape = type)) +
  geom_point() +
  scale_shape_manual(values = c(1, 23)) + 
  scale_color_manual(values = c("#666666", "black")) +
  scale_size(range = point.size.range) +
  scale_fill_manual(values = cbPalette[c(1,7)]) +
  geom_text(data = C.gallsize, aes(x = x, y = y, label = labels), inherit.aes = FALSE) +
  ylab("Leaf gall diameter (mm)") + xlab("") +
  theme_links

## Create a two-panel figure of dominant drivers of link composition and variation in leaf gall parasitism.
grid.arrange(gall.dominants.plot, leaf.gall.size, ncol = 2, 
             sub = textGrob("Willow genotype", vjust = -0.5, gp = gpar(cex = 1.25)))

## tidy up data for visualizing differences in link composition among genotypes.
link.df <- as.data.frame(tree_level_interaxn_all_plants_traits_size) %>%
  select(Genotype, 
         vLG_Platy, vLG_Mesopol, vLG_Tory, vLG_Eulo, vLG_Mymarid, 
         rG_Tory, rG_Eulo, rG_Platy, rG_Mesopol, rG_Lestodip, rG_Mesopol,
         aSG_Tory, SG_Platy) %>%
  gather(Genotype) %>%
  separate(variable, into = c("Gall", "Parasitoid"), sep = "_")

## revalue link composition data from species codes to common/species names
link.df$Gall <- factor(revalue(link.df$Gall, c("vLG" = "Leaf gall", 
                                               "rG" = "Bud gall",
                                               "aSG" = "apical-Stem gall",
                                               "SG" = "mid-Stem gall")))
link.df$Gall <- factor(link.df$Gall, levels = c("Leaf gall", "Bud gall", "apical-Stem gall", "mid-Stem gall"))
link.df$Parasitoid <- factor(revalue(link.df$Parasitoid, c("Platy" = "Platygaster",
                                                           "Mesopol" = "Mesopolobus",
                                                           "Tory" = "Torymus",
                                                           "Eulo" = "Eulophid",
                                                           "Lestodip" = "Lestodiplosis",
                                                           "Mymarid" = "Mymarid")))
link.df$Parasitoid <- factor(link.df$Parasitoid, levels = c("Platygaster",
                                                            "Mesopolobus",
                                                            "Torymus",
                                                            "Eulophid",
                                                            "Lestodiplosis",
                                                            "Mymarid"))
## order genotypes by mean leaf gall abundance
link.df$Genotype <- factor(link.df$Genotype, levels = vLG.df.sum$Genotype)

## link composition plot for all possible interactions. For supplementary materials.
link.composition.plot <- ggplot(link.df, aes(x = Genotype, y = value, fill = Parasitoid)) +
  facet_grid(Parasitoid ~ Gall) + 
  geom_boxplot() +
  scale_fill_manual(values = c(cbPalette[c(6,4,2,5,1)], "black")) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10,2)) +
  ylab("Link abundance (#/branch)") +
  theme_links

## link composition plot for Platygaster, Mesopolobous, and Torymus attacking leaf galls. These are the links that were the dominant drivers of variation in link composition among genotypes.
link.dominants.plot <- ggplot(filter(link.df, Gall == "Leaf gall", 
                           Parasitoid %in% c("Platygaster", "Mesopolobus", "Torymus")), 
                    aes(x = Genotype, y = value, fill = Parasitoid)) +
  facet_grid(Parasitoid ~ Gall) + 
  geom_boxplot() +
  geom_text(data = ABC.domptoids, aes(x = x, y = y, label = labels)) +
  scale_fill_manual(values = cbPalette[c(6,4,2)]) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10,2)) +
  ylab("Link abundance (#/branch)") + xlab("") +
  theme_links

## tidy up the data for visualizing the differences in leaf gall parasitism among willow genotypes.
vLG_ptized.df <- as.data.frame(tree_level_interaxn_all_plants_traits_size) %>%
  mutate(vLG_ptized = vLG_Platy + vLG_Tory + vLG_Mesopol + vLG_Eulo + vLG_Mymarid) %>%
  mutate(vLG_prop.ptized = vLG_ptized/vLG_abund, type = "mean") %>%
  select(Genotype, vLG_prop.ptized, vLG_abund, type)

glm.vLG_ptized <- glm(vLG_prop.ptized ~ Genotype, data = vLG_ptized.df, weights = vLG_abund, family = "quasibinomial")
summary(glm.vLG_ptized)
anova(glm.vLG_ptized, test = "F")

## calculated the weighted mean of leaf gall parasitism on each genotype. I used the number of galls sampled on each replicate willow as the weight, because these will have a more reliable estimate of mean gall size for that willow.
vLG_ptized.df.2 <- vLG_ptized.df %>%
  group_by(Genotype) %>%
  summarise(vLG_prop.ptized = weighted.mean(vLG_prop.ptized, vLG_abund, 
                                            na.rm = TRUE),
            vLG_abund = 15) %>% # point size to use for plotting
  mutate(type = "weighted.mean")

## bind the original and summarized data together, then order the genotypes according to median leaf gall abundance
vLG_ptized.df.3 <- rbind.data.frame(vLG_ptized.df.2, vLG_ptized.df)
vLG_ptized.df.3$Genotype <- factor(vLG_ptized.df.3$Genotype, levels = vLG.df.sum$Genotype)

## plot leaf gall parasitism among willow genotypes. Size of points correspond to weights used to calculated the weighted mean of leaf gall parasitism for each genotype.
vLG_ptized.plot <- ggplot(vLG_ptized.df.3, aes(x = Genotype, y = vLG_prop.ptized, 
                            size = vLG_abund, fill = type, color = type, shape = type)) +
  geom_point() +
  scale_shape_manual(values = c(1, 23)) + 
  scale_color_manual(values = c(cbPalette[1], "black")) +
  scale_size(range = point.size.range) +
  scale_fill_manual(values = cbPalette[c(1,7)]) +
  geom_text(data = D.ptism, aes(x = x, y = y, label = labels), inherit.aes = FALSE) +
  ylab("Proportion of leaf galls parasitized") + xlab("") +
  theme_links


grid.arrange(link.dominants.plot, vLG_ptized.plot, ncol = 2,
             sub = textGrob("Willow genotype", vjust = -0.5, gp = gpar(cex = 1.25)))



## create a focused dataset of how the probability of parasitoid attack of the leaf gall depends on gall size.
attack.df <- as.data.frame(tree_level_interaxn_all_plants_traits_size) %>% 
  filter(vLG.height.mean > 0) %>%
  #mutate(vLG_Platy.prop = vLG_Platy/vLG_abund,
   #      vLG_Mesopol.prop = vLG_Mesopol/vLG_abund,
    #     vLG_Tory.prop = vLG_Tory/vLG_abund) %>%
  select(vLG.height.mean, vLG_abund, vLG_Platy, vLG_Mesopol, vLG_Tory) %>%
  gather(vLG.height.mean, vLG_abund)

platy <- glm(value/vLG_abund ~ vLG.height.mean, data = filter(attack.df, variable == "vLG_Platy"), family = 'binomial', weights = vLG_abund)
summary(platy)
anova(platy, test = "Chi")

mesopol <- glm(value/vLG_abund ~ vLG.height.mean, data = filter(attack.df, variable == "vLG_Mesopol"), family = 'binomial', weights = vLG_abund)
summary(mesopol)
anova(mesopol, test = "Chi")

tory <- glm(value/vLG_abund ~ vLG.height.mean, data = filter(attack.df, variable == "vLG_Tory"), family = 'binomial', weights = vLG_abund)
summary(tory)
anova(tory, test = "Chi")

attack.plots <- ggplot(attack.df, aes(x = vLG.height.mean, y = value/vLG_abund, 
                      color = variable, shape = variable, linetype = variable)) + 
  geom_point(aes(size = vLG_abund), alpha = 0.5) +
  geom_smooth(method = "glm", family = binomial, 
              aes(weight = vLG_abund), se = FALSE, size = line.widths) +
  geom_text(data = B.ptism, aes(x = x, y = y, label = labels), inherit.aes = FALSE) +
  xlab("Leaf gall diameter (mm)") + ylab("Proportion parasitized") +
  scale_color_manual(values = cbPalette[c(6,4,2)]) + 
  theme_links 

library(MASS) # necessary for method = "glm.nb"
vLG_link.plot <- ggplot(attack.df, aes(x = vLG_abund, y = value, 
                                       color = variable, linetype = variable, shape = variable)) + 
  geom_jitter(position = position_jitter(width = 0.25, height = 0.25), 
              alpha = 0.5, size = 3) +
  geom_smooth(method = "glm.nb", formula = y ~ log(x), se = FALSE, size = line.widths) +
  geom_text(data = A.linkabund, aes(x = x, y = y, label = labels), inherit.aes = FALSE) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10,2)) +
  scale_color_manual(values = cbPalette[c(6,4,2)]) + 
  xlab("Leaf gall abundance (#/branch)") + ylab("Link abundance (#/branch)") + theme_links

## combine plots together
## resize graphs and arrange for a multipanel figure. Note that I manually added plot labels in a different program because I had difficulty at getting them all to have the appropriate spacing.
gp1<- ggplot_gtable(ggplot_build(vLG_link.plot))
gp2<- ggplot_gtable(ggplot_build(attack.plots))

maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3])

gp1$widths[2:3] <- maxWidth
gp2$widths[2:3] <- maxWidth

grid.arrange(gp1, gp2, ncol = 1)

### could use the code below for a supplement.
## try individual species plots for bud gall links
rG.df2 <- as.data.frame(tree_level_interaxn_all_plants_traits_size) %>%
  dplyr::select(rG_abund, rG_Tory, rG_Eulo, rG_Platy) %>%
  filter(rG_abund > 0) #%>%
  #gather(rG_abund)

rG_Platy.glm <- glm(rG_Platy ~ log(rG_abund), rG.df2, family = poisson)
summary(rG_Platy.glm) 
#plot(rG_Platy.glm) # not great

rG_Tory.glm <- glm(rG_Tory ~ log(rG_abund), rG.df2, family = poisson)
summary(rG_Tory.glm)
#plot(rG_Tory.glm) # okay

rG_Eulo.glm <- glm(rG_Eulo ~ log(rG_abund), rG.df2, family = poisson)
summary(rG_Eulo.glm)

rG_link.plot.df <- rG.df2 %>%
  gather(rG_abund)

rG_link.plot <- ggplot(rG_link.plot.df, aes(x = rG_abund, y = value, 
                                       color = variable, fill = variable, linetype = variable, shape = variable)) + 
  geom_jitter(position = position_jitter(width = 0.25, height = 0.25), 
              size = 3, color = "black") +
  geom_smooth(method = "glm", family = "poisson", formula = y ~ log(x), se = FALSE, size = line.widths) +
  scale_color_manual(values = cbbPalette[c(4,6,7)]) + 
  scale_fill_manual(values = cbbPalette[c(4,6,7)]) + 
  scale_shape_manual(values = c(21,22,24)) +
  xlab("Bud gall density (#/branch)") + ylab("Link density (#/branch)") + theme_links

## combine plots together
grid.arrange(vLG_link.plot, rG_link.plot, attack.plots, ncol = 3)

#### Ordination of link composition of the entire food web.

## tidy up full link composition data
full.links.df <- as.data.frame(tree_level_interaxn_all_plants_traits_size) %>%
  select(Genotype, 
         vLG_abund, aSG_abund, rG_abund, SG_abund, # willow-gall links
         vLG_Platy, vLG_Tory, vLG_Mesopol, vLG_Eulo, vLG_Mymarid, # vLG-ptoid links
         rG_Platy, rG_Tory, rG_Mesopol, rG_Eulo, rG_Lestodip, # rG-ptoid links
         SG_Platy, aSG_Tory) # other gall-ptoid links
all.links <- names(full.links.df)[-1]
trees.with.links <- which(rowSums(full.links.df[ ,all.links]) > 0)
table(full.links.df$Genotype[trees.with.links]) # J, N, and U have 2 or less replicates with any links.

full.links.df.sub <- filter(full.links.df[trees.with.links, ], 
                            Genotype != "J",
                            Genotype != "N",
                            Genotype != "U")

## perform RDA analyses and extract centroid scores for plotting.
library(vegan)
cap.geno <- capscale(full.links.df.sub[ ,all.links] ~ Genotype, 
                     data = full.links.df.sub,
                     distance = "bray")
summary(cap.geno)
centroids.cap.geno <- scores(cap.geno, choices = c(1,2), display = "cn")
rownames(centroids.cap.geno) <- levels(full.links.df$Genotype)[-c(10,14,21)] # remove Genotypes J, N, and U

## generate plot. Coded out the manner in which I saved the plot

#png(file="~/Documents/Genotype_Networks/figures/full.link.composition.png",width=8, height=8, units = "in", res = 300)
plot.new()
plot.window( xlim=c(-1.75,1.5), ylim = c(-1.5,2.5), asp = 1, cex.axis = 0.8) 
abline(h=0, lty="dotted")
abline(v=0, lty="dotted")
axis(side=1, tck = -0.015, padj=-0.5)
axis(side=2, tck= -0.015, las=1)
title(xlab="CAP 1", ylab="CAP 2", line = 1.75)
box()

# fill in the plot
ordiellipse(cap.geno, groups = full.links.df.sub$Genotype, 
            kind = "se", draw = "polygon",
            col= "gray50", #"gainsboro", 
            border = NA)
#ordiellipse(cap.geno, groups = full.links.df.sub$Genotype, kind = "se", show.groups="X", col = "grey", draw = "polygon") # code for highlighting specific groups
text(x = centroids.cap.geno[ ,1], centroids.cap.geno[ ,2], labels = row.names(centroids.cap.geno))
text(x = -1.5, y = 3, labels = "(B)")

#dev.off() # turn off pdf device.



