## This R script organizes the data for gall-parasitoid networks on individual trees of different willow genotypes.
## Author: Matt Barbour
## Email: barbour@zoology.ubc.ca

## upload packages for data manipulation. order of libraries is important. always load dplyr last.
library(reshape2)
library(reshape)
library(plyr)
library(dplyr)

## upload molten gall network data and calculate some basic summary data
gall_net_melt <- read.csv("~/Documents/Genotype_Networks/data/gall_network_data.csv")
gall_net_melt <- tbl_df(gall_net_melt) # 1,495 galls
sum(gall_net_melt$value, na.rm = TRUE) # 1709 total samples
sum(filter(gall_net_melt, gall_contents == "nothing")$value) # 518 "nothing" samples
sum(filter(gall_net_melt, gall_contents == "exit.hole")$value) # 136 exit holes. These indicate ectoparasitism, but I was unable to reliably determine attack by the particular parasitoid.

## create a data frame for the estimated number of shoots surveyed on each willow.
shootEsts_df <- gall_net_melt %>%
  mutate(plant.position = as.factor(plant.position)) %>%
  group_by(Gender, Genotype, plant.position) %>%
  summarise(shootEst.no18 = mean(shootEst.no18), shootEst.all = mean(shootEst.all), galls.found = mean(Galls.found)) # needed to take the mean because these shoot estimates are duplicated throughout the molten dataframe
shootEsts_df <- mutate(shootEsts_df, plant.position = as.character(plant.position))

## create a shoot data frame resolved for each genotype. May not need this later
#shootEsts_genotype_df <- shootEsts_df %>%
 # group_by(Genotype) %>%
  #summarise(n = n(), shootEst.no18 = sum(shootEst.no18), shootEst.all = sum(shootEst.all))

## create the tree-level interaction data
tree_interaxns_filter <- gall_net_melt %>%
  filter(gall_contents %in% c("aSG.larv", "Pont.ad", "Pont.prep", "rG.larv", "SG.larv", "vLG.pupa", "Eulo.fem", "Eulo.mal", "Eury.fem", "Eury.mal", "Lathro.fem", "Lathro.mal", "Lestodip", "Mesopol", "Mymarid", "Platy", "Ptero.2", "Tory.fem", "Tory.mal")) # 747 galls. I'm only retaining gall contents for which we could reliably determine the source of mortality or survival. 
sum(tree_interaxns_filter$value) # 973 individuals collected from these galls.

tree_interaxns_filter_add <- mutate(tree_interaxns_filter, 
                                    gall_contents_collapse = revalue(gall_contents, 
                                                                     c("Pont.ad" = "Pont.surv", 
                                                                       "Pont.prep" = "Pont.surv", 
                                                                       "Eulo.fem" = "Eulo", 
                                                                       "Eulo.mal" = "Eulo", 
                                                                       "Eury.fem" = "Eury", 
                                                                       "Eury.mal" = "Eury", 
                                                                       "Lathro.fem" = "Lathro", 
                                                                       "Lathro.mal" = "Lathro", 
                                                                       "Ptero.2" = "Mesopol", 
                                                                       "Tory.fem" = "Tory", 
                                                                       "Tory.mal" = "Tory") )) # this collapses the different parasitoid sexes into the same species. I also consider seeing Pontania adult and Pontania prepupa as survival.

##### The code below organizes the molten gall network data into a data frame with each plant as a replicate and is therefore the most straightforward for most analyses.
tree_level_interaxn_df <- dcast(tree_interaxns_filter_add, 
                                Gender + Genotype + plant.position ~ gall.sp + gall_contents_collapse,
                                sum) # reshape the data into a data frame for analysis

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
                                     gall_total_abund = aSG_abund + rG_abund + vLG_abund + SG_abund + rsLG_abund) # adding new columns of total abundance to the data frame.

## had to add these separately because rowSums(select()) wasn't working correctly with different sequences of columns to combine.
tree_level_interaxn_df_add$link_abund <- rowSums(select(tree_level_interaxn_df, 
                                                        aSG_Tory:rG_Platy, rG_Tory:rsLG_Lathro,
                                                        SG_Platy, vLG_Eulo:vLG_Tory)) # same as parasitoid abundance, but wrote this way to maintain consistency with how richness is calculated. 
tree_level_interaxn_df_add$link_richness <- rowSums(select(tree_level_interaxn_df_add, 
                                                           aSG_Tory:rG_Platy, rG_Tory:rsLG_Lathro,
                                                           SG_Platy, vLG_Eulo:vLG_Tory) > 0)
tree_level_interaxn_df_add$gall_total_rich <- rowSums(select(tree_level_interaxn_df_add, 
                                                             aSG_abund:rsLG_abund) > 0)
tree_level_interaxn_df_add$ptoid_total_rich <- rowSums(select(tree_level_interaxn_df_add, 
                                                              Platy_abund:Mymarid_abund) > 0)

tree_level_interaxn_df_add <- mutate(tree_level_interaxn_df_add, 
                                     network_richness = gall_total_rich + ptoid_total_rich, 
                                     #gall_survive_abund = aSG_aSG.larv + rG_rG.larv + 
                                       #rsLG_Pont.surv + SG_SG.larv + vLG_vLG.pupa, # not clear why I would need this
                                     vLG_parasitized = vLG_Eulo + vLG_Mesopol + 
                                       vLG_Mymarid + vLG_Platy + vLG_Tory,
                                     rG_parasitized = rG_Eulo + rG_Lestodip + rG_Mesopol + rG_Platy + rG_Tory,
                                     rsLG_parasitized = rsLG_Lathro + rsLG_Eury,
                                     SG_parasitized = SG_Platy,
                                     aSG_parasitized = aSG_Tory,
                                     vLG_ecto = vLG_Eulo + vLG_Tory + vLG_Mesopol,
                                     rG_ecto = rG_Eulo + rG_Tory + rG_Mesopol,
                                     vLG_egg = vLG_Platy + vLG_Mymarid) 

## join the tree level interaction data with the estimated number of shoots sampled per tree.
tree_level_interaxn_df_joined <- left_join(shootEsts_df, tree_level_interaxn_df_add)
length(table(gall_net_melt$plant.position)) # double checked and the number of unique plant positions is preserved

## fill NAs with zeros that are biologically meaningful, since these trees were sampled but no galls were ever collected from them.
tree_level_NAs_to_fill <- tree_level_interaxn_df_joined %>%
  ungroup() %>%
  select(-(Gender:galls.found)) %>%
  as.matrix
tree_level_NAs_to_fill[is.na(tree_level_NAs_to_fill)] <- 0

tree_level_interaxn_all_plants <- cbind.data.frame(select(tree_level_interaxn_df_joined, 
                                                          Gender:galls.found), tree_level_NAs_to_fill)
tree_level_interaxn_all_plants <- tbl_df(tree_level_interaxn_all_plants)

#tree_level_interaxn_all_plants <- mutate(tree_level_interaxn_all_plants, 
 #                                        gall_total_densityNO18 = gall_total_abund/shootEst.no18,
  #                                       gall_total_densityALL = gall_total_abund/shootEst.no18)

## load plant trait data and gall data from 2011. Consider downloading the data straight from dryad once the repository becomes open.
tree_level_traits <- read.csv("~/Documents/Genotype_Networks/data/plant.trait.galls.2011.tree.df.csv")
tree_level_traits <- select(tree_level_traits, -X, Genotype, plant.position, vLG.2011 = vLG, rsLG.2011 = rsLG, Total_Area:flavanonOLES.PC1)
tree_level_traits$plant.position <- factor(tree_level_traits$plant.position)

## manage gall trait data. The code below resolves everything down to the gall level and not the larva level, since there may be multiple larva within a gall.
gall.size.df <- tree_interaxns_filter_add %>%
  group_by(Genotype, plant.position, gall.id.nest, gall.sp) %>% # this makes sure I don't double count galls with multiple larvae, and is also why I take the "average" below, which is just the same number.
  summarise(gall.height = mean(g.height))

#gall.contents.df <- tree_interaxns_filter_add %>%
 # dcast(gall.id.nest + gall.sp ~ gall_contents_collapse, sum) %>%
  #tbl_df()

#gall.size.contents.df <- left_join(gall.size.df, gall.contents.df)

## this data frame resolves everything to gall size at the tree level.
tree_level_gall.size.df <- gall.size.df %>% # gall.size.contents.df
  group_by(Genotype, plant.position, gall.sp) %>% 
  summarise(gall.height.mean = mean(gall.height, na.rm = TRUE), gall.count = n())

## cast the data frame so each row is a unique plant position.
tree_level_gall.size.cast <- dcast(tree_level_gall.size.df, Genotype + plant.position ~ gall.sp, 
                                   value.var = "gall.height.mean") %>% 
  select(Genotype, plant.position, 
         aSG.height.mean = aSG, rG.height.mean = rG, rsLG.height.mean = rsLG, 
         SG.height.mean = SG, vLG.height.mean = vLG) %>%
  tbl_df()

## same as above but retaining data fro gall count per tree, which is relevant in estimating the accuracy of each mean estimate.
tree_level_gall.count.cast <- dcast(tree_level_gall.size.df, Genotype + plant.position ~ gall.sp,
                                    value.var = "gall.count") %>% 
  select(Genotype, plant.position, 
         aSG.gall.count = aSG, rG.gall.count = rG, rsLG.gall.count = rsLG, 
         SG.gall.count = SG, vLG.gall.count = vLG) %>%
  tbl_df()

tree_level_gall.size.count.join <- left_join(tree_level_gall.size.cast, tree_level_gall.count.cast)

## merge interaction and tree-trait data
tree_level_interaxn_all_plants_traits_size <- join_all(list(tree_level_interaxn_all_plants, tree_level_traits, tree_level_gall.size.count.join), by = "plant.position") %>%
  tbl_df()

write.csv(tree_level_interaxn_all_plants_traits_size, '~/Documents/Genotype_Networks/data/tree_level_interaxn_all_plants_traits_size.csv')

## for reference in analyses
interaxns_gallsurv <- tree_level_interaxn_all_plants_traits_size %>% 
  select(aSG_aSG.larv:vLG_vLG.pupa) %>%
  names()

interaxns_gallsurv_noPont <- tree_level_interaxn_all_plants_traits_size %>% 
  select(aSG_aSG.larv:rG_Tory, SG_Platy:vLG_vLG.pupa) %>%
  names()

interaxns <- tree_level_interaxn_all_plants_traits_size %>%
  select(aSG_Tory:rG_Platy, rG_Tory:rsLG_Lathro, SG_Platy, vLG_Eulo:vLG_Tory) %>%
  names()

interaxns_noPont <- tree_level_interaxn_all_plants_traits_size %>%
  select(aSG_Tory:rG_Platy, rG_Tory, SG_Platy, vLG_Eulo:vLG_Tory) %>%
  names()

interaxns_guild <- tree_level_interaxn_all_plants_traits_size %>%
  select(aSG_Tory, rG_ecto, rG_Lestodip, vLG_ecto, vLG_egg, SG_Platy) %>%
  names()

### Getting data for estimateS
#estimateS_df <- select(tree_level_interaxn_all_plants, Genotype, plant.position, aSG_Tory:rG_Platy, rG_Tory:rsLG_Lathro, SG_Platy, vLG_Eulo:vLG_Tory)

#write.csv(t(estimateS_df), "~/Documents/Genotype_Networks/data/estimateS_gall_network_data.csv")
