require(ggplot2)
require(gridExtra)
require(dplyr)

source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')

## create theme for figures
theme_galls <- theme_classic() +
  theme(axis.title.x = element_text(size = 16, vjust = 0.1), 
        axis.title.y = element_text(size = 16, vjust = 1),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1))

## Gall community composition figure

# calculate mean abundance of each gall species for each genotype
gall.comm.df <- tree_level_interaxn_all_plants_traits_size %>%
  select(Genotype,
         "Iteomyia salicisverruca" = vLG_abund,
         "Rabdophaga salicisbrassicoides" = rG_abund,
         "Cecidomyiid sp. A" = aSG_abund,
         "Rabdophaga salicisbattatus" = SG_abund) %>%
  group_by(Genotype) %>%
  summarise_each(funs(mean.narm = mean(., na.rm = TRUE)))
levels(gall.comm.df$Genotype)[1] <- "C" # changed Genotype * to C for plotting aesthetics
gall.comm.df$Genotype <- factor(gall.comm.df$Genotype, levels = LETTERS[1:26])  # reorder for plotting A->Z

gall.comm.melt <- melt(gall.comm.df) # melt for plotting in ggplot

# create figure
comm.plot <- ggplot(gall.comm.melt, aes(x = Genotype, y = value, fill = variable)) + 
  geom_bar(stat = "identity", color = "black") +
  scale_fill_grey(start = 0, end = 1) +
  scale_y_continuous(expand = c(0,0), limits = c(0,20)) +
  ylab("Gall density (#/branch)") +
  xlab("") +
  theme_galls +
  annotate("text", x = 1, y = 19, label = "(A)")

## Variation in Iteomyia gall size figure
gall.size.df <- tree_level_interaxn_all_plants_traits_size 
levels(gall.size.df$Genotype)[1] <- "C" # changed Genotype * to C for plotting aesthetics
gall.size.df$Genotype <- factor(gall.size.df$Genotype, levels = LETTERS[1:26]) # reorder for plotting A->Z

size.plot <- ggplot(gall.size.df, aes(x = Genotype, y = vLG.height.mean)) + 
  geom_boxplot() +
  ylab("Iteomyia gall diameter (mm)") +
  theme_galls +
  annotate("text", x = 1, y = 14, label = "(B)")

## arrange figure
grid.arrange(comm.plot, size.plot) # saved manually as 'gall_community_size_plot.pdf'

## calculate summary statistics for gall community and gall size data

# basic summary data for the gall community
summary(select(gall.comm.df, -Genotype)) 

# calculate average dissimilarity in gall community composition among willow genotypes using average gall abundance data.
require(vegan)
gall.comm.df.noU <- gall.comm.df %>% filter(Genotype != "U")
gall.comm.bray <- vegdist(select(gall.comm.df.noU, -Genotype), method = "bray")
mean(gall.comm.bray) # on average, gall communities were 52% dissimilar from each other.
sd(gall.comm.bray)
range(gall.comm.bray)

# compare to average dissimilarity estimates using dissimilarity index at tree level.
galls <- full.df %>%
  select(Genotype, vLG_abund, rG_abund, aSG_abund, SG_abund) %>%
  mutate(total_abund = vLG_abund + rG_abund + aSG_abund + SG_abund) %>%
  filter(total_abund > 0)

galls.geno.dist <- meandist(vegdist(galls[ ,-1], "bray"), galls$Genotype)
summary(galls.geno.dist) # on average, gall communities were 57% dissimilar from each other
#plot(galls.geno.dist, cluster = "average")

# 
gall.size.summary <- gall.size.df %>%
  group_by(Genotype) %>%
  summarise(Iteomyia_mean_diameter.mm = mean(vLG.height.mean, na.rm = TRUE))
summary(gall.size.summary)
11.025/4.83 # 2-fold variation in Iteomyia size among willow genotypes
