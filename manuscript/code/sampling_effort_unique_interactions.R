#############################################################
#  Description: This script estimates the total number of unique gall-parasitoid interactions based on my sampling effort (results presented in supplementary material).
#  Code author: Matthew A. Barbour
#  Email: barbour@zoology.ubc.ca
#############################################################

## load required data
tree_level_interaxn_all_plants_traits_size <- read.csv("manuscript/Dryad_data_copy/empirical_data/tree_level_interaxn_all_plants_traits_size.csv")

## estimating unique number of gall-parasitoid interactions
library(vegan)
net.df <- select(tree_level_interaxn_all_plants_traits_size, aSG_Tory:rG_Platy, rG_Tory, SG_Platy, vLG_Eulo:vLG_Tory)

specpool(net.df)

