## This script manages network data at genotype level

## source in required scripts and data
source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')
source('~/Documents/miscellaneous_R/length2_function.R')

## create genotype-level interaction data. I accounted for density by dividing by the number of estimated shoots, because that was the most straightforward way to account for differences in sampling among trees. RIGHT NOW I HAVEN'T CALCULATE RICHNESS RESPONSES AND PARASITISM YET. OR 2011 GALL DATA. I SHOULD NOT CONTROL FOR RICHNESS OR PARASITISM IN TERMS OF SHOOT COUNTS, BUT I SHOULD DIVIDE 2011 GALL DATA BY TOTAL AREA.
tree_level_interaxn_density <- select(tree_level_interaxn_all_plants_traits_size, aSG_aSG.larv:link_abund, vLG_ecto:rG_ecto, vLG_egg)/tree_level_interaxn_all_plants$shootEst.no18*100 # scale to a per 100 shoot basis. 

geno_level_interaxn_density <- aggregate(tree_level_interaxn_density, 
                                         by = list(Genotype = tree_level_interaxn_all_plants_traits_size$Genotype), 
                                         FUN = mean) %>%
  tbl_df()

## create genotype-level trait data.
geno_level_traits <- tree_level_traits %>%
  group_by(Genotype) %>%
  select(Total_Area:flavanonOLES.PC1) %>%
  summarise_each(funs(mean.na.rm = mean(., na.rm = TRUE))) 

## create genotype-level gall size data
geno_level_gall.size.count.join <- select(tree_level_gall.size.count.join, 
                                          Genotype:vLG.height.mean) %>% # counts per tree not needed
  group_by(Genotype) %>% 
  select(-plant.position) %>%
  summarise_each(funs(mean.na.rm = mean(., na.rm = TRUE), 
                      length2.na.rm = length2(., na.rm = TRUE))) # length2 enables me to calculate the number of plant replicates that were used to estimate mean gall size.

## merge interaction and genotype-trait data
genotype_level_interaxn_traits_size <- join_all(list(geno_level_interaxn_density, geno_level_traits, geno_level_gall.size.count.join), by = "Genotype") %>%
  tbl_df()

#tree_level_interaxn_all_plants_traits_size <- read.csv('~/Documents/Genotype_Networks/data/tree_level_interaxn_all_plants_traits_size.csv')


## EVERYTHING BELOW THIS IS OLD

tree_level_interaxn_all_plants <- read.csv('~/Documents/Genotype_Networks/data/tree_level_interaxn_all_plants.csv')

### Genotype level data including larva survival and shoot estimates as separate variables
genotype_level.df <- tree_level_interaxn_all_plants %>%
  group_by(Genotype) %>%
  select(shootEst.no18, shootEst.all, aSG_aSG.larv:gall_survive_abund) %>%
  summarise_each(funs(sum)) %>%
  ungroup()

write.csv(genotype_level.df, "~/Documents/Genotype_Networks/data/genotype_level.df.csv")


### Modularity and nestedness analysis. 
modularity_df <- tree_level_interaxn_all_plants %>%
  group_by(Genotype) %>%
  select(aSG_Tory:rG_Platy, rG_Tory:rsLG_Lathro, SG_Platy, vLG_Eulo:vLG_Tory, shootEst.no18, shootEst.all) %>%
  filter(Genotype != "U") %>% # removing because there were no observed interactions
  summarise_each(funs(sum)) %>%
  ungroup()

rownames(modularity_df) <- as.character(modularity_df$Genotype)
modularity_df_interaxns <- modularity_df %>%
  select(-Genotype) %>%
  apply(MARGIN = 2, FUN = function(x) round(x/modularity_df$shootEst.no18*1000)) # standardize the frequency of all interactions by dividing by the number of shoots sampled, mulitplying it by 1000, then rounding to the nearest whole integer. Therefore, this represents the frequency of each interaction per 1000 shoots sampled. I chose 1000 shoots because it was near the minimum number of shoots sampled per genotype.
modularity_df_interaxns <- as.data.frame(modularity_df_interaxns)
modularity_df_interaxns <- select(modularity_df_interaxns, aSG_Tory:vLG_Tory)

write.csv(modularity_df_interaxns, "~/Documents/Genotype_Networks/data/genotype_gall_parasitoid_network_data.csv")

