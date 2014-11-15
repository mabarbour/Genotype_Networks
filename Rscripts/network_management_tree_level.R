######## Tree level Network, Community, and Species level analysis

#### upload necessary packages

# analysis
library(mvabund)
library(vegan)
library(visreg)
library(bipartite)

# data manipulation. order of libraries is important. always load dplyr last.
library(reshape2)
library(reshape)
library(dplyr)

# plotting
library(ggplot2)


# upload molten gall network data
gall_net_melt <- read.csv("~/Documents/Genotype_Networks/data/gall_network_data.csv")
gall_net_melt <- tbl_df(gall_net_melt) # 1,511 galls
sum(gall_net_melt$value[1:1511]) # 1740 total samples
sum(filter(gall_net_melt, gall_contents == "nothing")$value) # 530 "nothing" samples
sum(filter(gall_net_melt, gall_contents == "exit.hole")$value) # 136

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


# create the tree-level data

tree_interaxns_filter <- gall_net_melt %>%
  filter(gall_contents %in% c("aSG.larv", "Pont.ad", "Pont.prep", "rG.larv", "SG.larv", "vLG.pupa", "Eulo.fem", "Eulo.mal", "Eury.fem", "Eury.mal", "Lathro.fem", "Lathro.mal", "Lestodip", "Mesopol", "Mymarid", "Platy", "Ptero.2", "Tory.fem", "Tory.mal")) # 764 galls
sum(tree_interaxns_filter$value) # 991 individuals

tree_interaxns_filter_add <- mutate(tree_interaxns_filter, gall_contents_collapse = revalue(gall_contents, c("Pont.ad" = "Pont.surv", "Pont.prep" = "Pont.surv", "Eulo.fem" = "Eulo", "Eulo.mal" = "Eulo", "Eury.fem" = "Eury", "Eury.mal" = "Eury", "Lathro.fem" = "Lathro", "Lathro.mal" = "Lathro", "Ptero.2" = "Mesopol", "Tory.fem" = "Tory", "Tory.mal" = "Tory") ))

tree_level_interaxn_df <- dcast(tree_interaxns_filter_add, Gender + Genotype + plant.position ~ gall.sp + gall_contents_collapse, sum) 

tree_level_interaxn_df_add <- mutate(tree_level_interaxn_df,
                                     plant.position = as.character(plant.position),
                                     aSG_abund = rowSums(select(tree_level_interaxn_df,
                                                                starts_with("aSG"))),
                                     rG_abund = rowSums(select(tree_level_interaxn_df,
                                                               starts_with("rG"))),
                                     vLG_abund = rowSums(select(tree_level_interaxn_df,
                                                                starts_with("vLG"))),
                                     SG_abund = rowSums(select(tree_level_interaxn_df, 
                                                               starts_with("SG"))),
                                     rsLG_abund = rowSums(select(tree_level_interaxn_df,
                                                                 starts_with("rsLG"))),
                                     Platy_abund = rowSums(select(tree_level_interaxn_df,
                                                                  ends_with("Platy"))),
                                     Tory_abund =  rowSums(select(tree_level_interaxn_df,
                                                                  ends_with("Tory"))),
                                     Mesopol_abund =  rowSums(select(tree_level_interaxn_df,
                                                                     ends_with("Mesopol"))),
                                     Eulo_abund = rowSums(select(tree_level_interaxn_df,
                                                                 ends_with("Eulo"))),
                                     Lestodip_abund = rG_Lestodip,
                                     Eury_abund = rsLG_Eury,
                                     Lathro_abund = rsLG_Lathro,
                                     Mymarid_abund = vLG_Mymarid,
                                     gall_total_abund = aSG_abund + rG_abund + vLG_abund + SG_abund + rsLG_abund)

# had to add separately, because rowSums function wasn't acting funky if I didn't
tree_level_interaxn_df_add$link_abund <- rowSums(select(tree_level_interaxn_df, aSG_Tory:rG_Platy, rG_Tory:rsLG_Lathro, SG_Platy, vLG_Eulo:vLG_Tory)) # same as parasitoid abundance, but wrote this way to maintain consistency with how richness is calculated.
tree_level_interaxn_df_add$link_richness <- rowSums(select(tree_level_interaxn_df_add, aSG_Tory:rG_Platy, rG_Tory:rsLG_Lathro, SG_Platy, vLG_Eulo:vLG_Tory) > 0)
tree_level_interaxn_df_add$gall_total_rich <- rowSums(select(tree_level_interaxn_df_add, aSG_abund:rsLG_abund) > 0)
tree_level_interaxn_df_add$ptoid_total_rich <- rowSums(select(tree_level_interaxn_df_add, Platy_abund:Mymarid_abund) > 0)

tree_level_interaxn_df_add <- mutate(tree_level_interaxn_df_add, network_richness = gall_total_rich + ptoid_total_rich, gall_survive_abund = aSG_aSG.larv + rG_rG.larv + rsLG_Pont.surv + SG_SG.larv + vLG_vLG.pupa, vLG_parasitized = vLG_Eulo + vLG_Mesopol + vLG_Mymarid + vLG_Platy + vLG_Tory) # note that gall_total_abund would be the equivalent of "network abundance", since the abundance data include parasitized individuals.
# still need to add trees with no galls collected!!!

tree_level_interaxn_df_joined <- left_join(shootEsts_df, tree_level_interaxn_df_add)
length(table(gall_net_melt$plant.position)) # double checked and the number of unique plant positions is preserved

tree_level_NAs_to_fill <- tree_level_interaxn_df_joined %>%
  ungroup() %>%
  select(-(Gender:galls.found)) %>%
  as.matrix
tree_level_NAs_to_fill[is.na(tree_level_NAs_to_fill)] <- 0 # these zeros are biologically meaningful, since these trees were sampled, but no galls were ever collected from them.
tree_level_interaxn_all_plants <- cbind.data.frame(select(tree_level_interaxn_df_joined, Gender:galls.found), tree_level_NAs_to_fill)
tree_level_interaxn_all_plants <- tbl_df(tree_level_interaxn_all_plants)

tree_level_interaxn_all_plants <- mutate(tree_level_interaxn_all_plants, gall_total_densityNO18 = gall_total_abund/shootEst.no18, gall_total_densityALL = gall_total_abund/shootEst.no18)

write.csv(tree_level_interaxn_all_plants, '~/Documents/Genotype_Networks/data/tree_level_interaxn_all_plants.csv')

### Getting data for estimateS
estimateS_df <- select(tree_level_interaxn_all_plants, Genotype, plant.position, aSG_Tory:rG_Platy, rG_Tory:rsLG_Lathro, SG_Platy, vLG_Eulo:vLG_Tory)

write.csv(t(estimateS_df), "~/Documents/Genotype_Networks/data/estimateS_gall_network_data.csv")
