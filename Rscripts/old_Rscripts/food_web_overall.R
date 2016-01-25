############# Overall Gall-Parasitoid Network Analysis

# load libraries
library(reshape2)
library(reshape)
library(dplyr)
library(bipartite)
library(RColorBrewer)

#source('~/Documents/ggnet/bipartite_plot_info.R')
#source('~/Documents/ggnet/ggnet_bipartite.R')

# upload molten gall network data and manage it.
gall_net_melt <- read.csv("~/Documents/Genotype_Networks/data/gall_network_data.csv")
gall_net_melt <- tbl_df(gall_net_melt)

total_gall_parasitoid_network <- gall_net_melt %>%
  filter(gall_contents %in% c("aSG.larv", "Pont.ad", "Pont.prep", "rG.larv", "SG.larv", "vLG.pupa", "Eulo.fem", "Eulo.mal", "Eury.fem", "Eury.mal", "Lathro.fem", "Lathro.mal", "Lestodip", "Mesopol", "Mymarid", "Platy", "Ptero.2", "Tory.fem", "Tory.mal")) 

total_gall_parasitoid_network <- mutate(total_gall_parasitoid_network, gall_contents_collapse = revalue(gall_contents, c("Pont.ad" = "Pont.surv", "Pont.prep" = "Pont.surv", "Eulo.fem" = "Eulo", "Eulo.mal" = "Eulo", "Eury.fem" = "Eury", "Eury.mal" = "Eury", "Lathro.fem" = "Lathro", "Lathro.mal" = "Lathro", "Ptero.2" = "Mesopol", "Tory.fem" = "Tory", "Tory.mal" = "Tory") ))

gall.ptoid.net <- dcast(total_gall_parasitoid_network, gall.sp ~ gall_contents_collapse, sum)
gall.ptoid.net <- gall.ptoid.net %>%
  mutate(survive = aSG.larv + rG.larv + SG.larv + vLG.pupa + Pont.surv) %>%
  select(gall.sp, Eulo, Eury, Lathro, Lestodip, Mesopol, Mymarid, Platy, Tory, survive)
rownames(gall.ptoid.net) <- gall.ptoid.net$gall.sp
gall.ptoid.net$total <- rowSums(select(gall.ptoid.net, -gall.sp))

total.interactions <- sum(rowSums(select(gall.ptoid.net, Eulo:Tory)))
select(gall.ptoid.net, Eulo:Tory)/total.interactions
hist(as.numeric(as.matrix(select(gall.ptoid.net, Eulo:Tory)/total.interactions)))

specieslevel(select(gall.ptoid.net, Eulo:Tory), low.abun = gall.ptoid.net$survive)

# specialization of each parasitoid species
ptoid.abund <- rowSums(t.gall.ptoid.net)
ptoid.percents <- ptoid.abund/sum(ptoid.abund)
sum(ptoid.percents[c("Platy","Mesopol","Tory")]) # 84%

t.gall.ptoid.net <- t(select(gall.ptoid.net, Eulo:Tory))
t.gall.ptoid.net/rowSums(t.gall.ptoid.net)

# percent parasitism on each gall species
parasitism.df <- select(gall.ptoid.net, Eulo:Tory)/gall.ptoid.net$total
parasitism.df$gall.sp <- rownames(parasitism.df)
parasitism.df.melt <- melt(parasitism.df)
parasitism.df.melt$gall.sp <- factor(as.character(parasitism.df.melt$gall.sp), levels = c("rsLG","aSG","rG","SG","vLG"))
parasitism.df.melt$parasitoid <- factor(as.character(parasitism.df.melt$variable), levels = c("Platy","Tory","Mesopol","Eury","Eulo","Lathro","Lestodip","Mymarid"))

ggplot(parasitism.df.melt, aes(x = gall.sp, y = value*100, fill = variable)) + 
  geom_bar(stat = "identity") +
  ylab("Percent Parasitism") +
  xlab("Gall species") + 
  scale_fill_brewer(palette = "Paired") +
  theme_classic()

gall.ptism <- rowSums(select(parasitism.df, -gall.sp))
mean(gall.ptism[c("aSG","rG","rsLG","SG")]) # 14% average parasitism rate on 4 other gall species, whereas vLG received 42% parasitism. 
0.42105263/0.1375239 # 3.1 fold higher parasitism on vLG compared to other gall species.

gall.ptoid.net$total/(sum(gall.ptoid.net$total)) # vLG made up 44% of the galling insect community, followed by rG at 35%, rsLG at 12%, aSG at 4% and SG at 5%.

########### EVERYTHING BELOW THIS IS OLD
# shoot estimates on a genotype level
shootEsts_df <- gall_net_melt %>%
  mutate(plant.position = as.factor(plant.position)) %>%
  group_by(Gender, Genotype, plant.position) %>%
  summarise(shootEst.no18 = mean(shootEst.no18), shootEst.all = mean(shootEst.all), galls.found = mean(Galls.found)) # needed to take the mean because these shoot estimates are duplicated throughout the molten dataframe
shootEsts_df <- mutate(shootEsts_df, plant.position = as.character(plant.position))

shootEsts_genotype_df <- shootEsts_df %>%
  group_by(Genotype) %>%
  summarise(n = n(), shootEst.no18 = sum(shootEst.no18), shootEst.all = sum(shootEst.all))

genotype_gall_network <- dcast(total_gall_parasitoid_network, Genotype ~ gall.sp, sum)