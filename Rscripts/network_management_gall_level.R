### prepare data that is resolved to the individual gall level
source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')

gall.size.df <- tree_interaxns_filter_add %>%
  group_by(Genotype, plant.position, gall.id.nest, gall.sp) %>% # this makes sure I don't double count galls with multiple larvae, and is also why I take the "average" below, which is just the same number.
  summarise(gall.height = mean(g.height)) %>%
  ungroup()

gall.contents.df <- tree_interaxns_filter_add %>%
  dcast(gall.id.nest ~ gall.sp + gall_contents_collapse, sum) %>%
  tbl_df()

## merge data sets together and make some new columns of data
gall.size.contents.df <- left_join(gall.size.df, gall.contents.df) 
gall.size.contents.df <- mutate(gall.size.contents.df, plant.position = factor(plant.position))
gall.size.contents.df$vLG_total <- rowSums(select(gall.size.contents.df, starts_with("vLG")))
gall.size.contents.df$rG_total <- rowSums(select(gall.size.contents.df, starts_with("rG")))
gall.size.contents.df$rsLG_total <- rowSums(select(gall.size.contents.df, starts_with("rsLG")))
gall.size.contents.df$aSG_total <- rowSums(select(gall.size.contents.df, starts_with("aSG")))
gall.size.contents.df$SG_total <- rowSums(select(gall.size.contents.df, starts_with("SG")))


## create gall density data and plant architecture and trichome density. These were the only plant traits that I expected a priori to affects parasitoids. I would expect the other mechanisms to act indirectly by changing gall densities.
gall_density.df <- with(tree_level_interaxn_all_plants_traits_size, 
                        data.frame(plant.position = plant.position, 
                                   Genotype = Genotype, 
                                   Gender = Gender, 
                                   vLG_density = vLG_abund/shootEst.no18*100,
                                   rG_density = rG_abund/shootEst.no18*100,
                                   aSG_density = aSG_abund/shootEst.no18*100,
                                   rsLG_density = rsLG_abund/shootEst.no18*100,
                                   SG_density = SG_abund/shootEst.no18*100,
                                   Total_Area = Total_Area,
                                   Density = Density,
                                   Height = Height,
                                   D = D_mean_smoothed,
                                   Trichomes = Trichome.No.))

# marge gall size and gall density data for mixed effects modelling
gall_mech.df <- left_join(gall.size.contents.df, gall_density.df, by = "plant.position")


