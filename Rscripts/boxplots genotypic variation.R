## source in required data
source('~/Documents/Genotype_Networks/Rscripts/network_management_genotype_level.R')

## load required libraries for plotting
library(ggplot2)
library(grid)
library(gridExtra)

#### Generate boxplots

## create data frame for plotting
plot.df <- genotype_level_interaxn_traits_size %>%
  melt(id.vars = "Genotype", measure.vars = c("vLG_abund", "vLG.height.mean_mean.na.rm", "rG_abund", "rG.height.mean_mean.na.rm", "vLG_egg", "vLG_ecto", "rG_ecto")) %>%
  mutate(species = factor(c(rep("Iteomyia", 2*26), rep("Rabdophaga (bud)", 2*26), 
                     rep("Iteomyia - Egg parasitoids", 26), rep("Iteomyia - Larval parasitoids", 26),
                     rep("Rabdophaga (bud) - Larval parasitoids", 26))),
         variable.replace = c(rep("mean_gall.density", 26), rep("mean_height", 26), 
                              rep("mean_gall.density", 26), rep("mean_height", 26),
                              rep("mean_interaxn.density", 26*3)))
levels(plot.df$species) <- gsub("- ", "-\n", levels(plot.df$species))
#levels(plot.df$species) <- rev(levels(plot.df$species))
plot.df$species <- factor(plot.df$species, levels = rev(levels(plot.df$species)))

## create theme for boxplots
theme_boxplot <- theme_classic() + theme(axis.title = element_text(size = 14, vjust = -0.5),
                                         axis.text = element_text(size = 12))

## boxplots
gall.density.plot <- ggplot(data = filter(plot.df, variable.replace == "mean_gall.density"), aes(x = species, y = value)) +
  geom_boxplot() + ylab("Gall density (count per 100 shoots)") + xlab("") + theme_boxplot + coord_flip()

gall.size.plot <- ggplot(data = filter(plot.df, variable.replace == "mean_height"), aes(x = species, y = value)) +
  geom_boxplot() + ylab("Gall diameter (mm)") + xlab("") + theme_boxplot + coord_flip()

interaction.density.plot <- ggplot(data = filter(plot.df, variable.replace == "mean_interaxn.density"), 
       aes(x = species, y = value)) +
  geom_boxplot() + ylab("Interaction density (count per 100 shoots)") + xlab("") + theme_boxplot +
  theme(axis.text.x = element_text(size = 10)) + coord_flip()# , angle = 45, hjust = 1

## barplots of heritabilities
barplot.df <- data.frame(variable = factor(c("Density Iteomyia -\nEgg ptoids", #Density\n
                                      "Density Iteomyia -\nLarval ptoids",#Density\n
                                      "Density Rab. (bud) -\nLarval ptoids",#Density\n
                                      "Density\nIteomyia", "Density\nRab. (bud)",#Density\n
                                      "Diameter\nIteomyia", "Diameter\nRab. (bud)"), #Diameter\n
                                      levels = rev(c("Density Iteomyia -\nEgg ptoids", #Density\n
                                                 "Density Iteomyia -\nLarval ptoids",#Density\n
                                                 "Density Rab. (bud) -\nLarval ptoids",#Density\n
                                                 "Density\nIteomyia", "Density\nRab. (bud)",#Density\n
                                                 "Diameter\nIteomyia", "Diameter\nRab. (bud)"))),
                         H2 = c(0.31, 0.27, 0.11, 0.29, 0.17, 0.15, 0.04))

h2.plots <- ggplot(barplot.df, aes(x = variable, y = H2)) + geom_bar(stat = "identity") + 
  theme_boxplot + xlab("") + scale_y_continuous(expand = c(0,0), limits = c(0, 0.35)) + 
  ylab(expression(paste("Broad-sense Heritability (",italic(H^2),")"))) + 
  theme(axis.text.x = element_text(size = 10)) + coord_flip() # , angle = 45, hjust = 1

## resize graphs and arrange for a multipanel figure. 
gp1<- ggplot_gtable(ggplot_build(gall.density.plot))
gp2<- ggplot_gtable(ggplot_build(gall.size.plot))
gp3<- ggplot_gtable(ggplot_build(interaction.density.plot))
gp4 <- ggplot_gtable(ggplot_build(h2.plots))
maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3], gp3$widths[2:3], gp4$widths[2:3])
gp1$widths[2:3] <- maxWidth
gp2$widths[2:3] <- maxWidth
gp3$widths[2:3] <- maxWidth
gp4$widths[2:3] <- maxWidth

grid.arrange(gp3, gp4, gp1, gp2, ncol = 2)#, 
             #sub = textGrob(expression(paste(italic("S. hookeriana"), " genotype")), 
                     #       gp = gpar(cex = 1.5), vjust=0))

## old below
library(psych)
library(car)
corr.df <- select(genotype_level_interaxn_traits_size, vLG_egg, vLG_ecto, rG_ecto, vLG_abund, rG_abund, vLG.height.mean_mean.na.rm, rG.height.mean_mean.na.rm)

scatterplotMatrix((corr.df))

corr.test(corr.df)

plot(vLG_egg ~ vLG.height.mean_mean.na.rm, filter(genotype_level_interaxn_traits_size, vLG.height.mean_mean.na.rm > 5))
summary(lm(vLG_egg ~ vLG.height.mean_mean.na.rm, filter(genotype_level_interaxn_traits_size, vLG.height.mean_mean.na.rm > 5)))
summary(lm(vLG_ecto ~ vLG.height.mean_mean.na.rm + I(vLG.height.mean_mean.na.rm^2), genotype_level_interaxn_traits_size))
