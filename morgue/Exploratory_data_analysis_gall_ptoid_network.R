setwd("~/Documents/Genotype_Networks")
source('modularity.R')
library(bipartite)

willow_gall_ptoid_network_original <- read.csv("~/Documents/Genotype_Networks/Gall_Interaction_network_metaweb_March27_2013.csv")

willow_gall_ptoid_network <- as.matrix(willow_gall_ptoid_network_original[1:4, 2:8])
row_names_willow_gall_ptoid_network <- list(c("Iteomyia salicisverruca", "Rabdophaga salicisbrassicoides", "Pontania californica", "twG"))
dimnames(willow_gall_ptoid_network)[1] <- row_names_willow_gall_ptoid_network

modularity_bipartite(willow_gall_ptoid_network)
rel_modularity_bipartite(willow_gall_ptoid_network) # modularity is marginally significant with this early dataset.  I bet this will change once I find more parasitoid species.

visweb(willow_gall_ptoid_network, type="diagonal", frame=F, labsize=0.5)