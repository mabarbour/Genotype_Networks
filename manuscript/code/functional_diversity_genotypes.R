
## load required data
tree_level_interaxn_all_plants_traits_size <- read.csv("manuscript/Dryad_data_copy/empirical_data/tree_level_interaxn_all_plants_traits_size.csv")

## Calculate functional-trait diversity of willow genotypes
library(FD)
genotype.trait.means <- as.data.frame(tree_level_interaxn_all_plants_traits_size) %>%
  select(Genotype, Total_Area:HCH.tremulacin.__A220nm,water_content:N_imputed) %>%
  group_by(Genotype) %>%
  summarise_each(funs(mean.narm = mean(., na.rm = TRUE))) %>%
  ungroup() %>%
  as.data.frame()

FD.traits <- dbFD(genotype.trait.means[ ,-1])
FD.traits$FEve
FD.traits$FDiv
