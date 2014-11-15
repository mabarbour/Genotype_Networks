library(dplyr)

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

