source('~/Documents/vegan/R/nestednodf.R')
source('~/Documents/Genotype_Networks/Rscripts/null_model_analysis_WNODF_QuaBiMo.R')
source('~/Documents/Genotype_Networks/Rscripts/best_QuaBiMo_function.R')
source('~/Documents/Genotype_Willow_Community/datasets_&_Rscripts/beta.div.R')

library(dplyr)
library(bipartite)
library(GGally)
library(ggplot2)

# possibly required packages for plotting


genotype_gall_parasitoid_network <- read.csv('~/Documents/Genotype_Networks/data/genotype_gall_parasitoid_network_data.csv')
rownames(genotype_gall_parasitoid_network) <- genotype_gall_parasitoid_network$X 
genotype_gall_parasitoid_network <- select(genotype_gall_parasitoid_network, -X)


###### Plots of Contributions to network structure

### Genotypes
geno.contrib.df <- data.frame(contrib = rowSums(genotype_gall_parasitoid_network),
                              percent.contrib = rowSums(genotype_gall_parasitoid_network)/sum(rowSums(genotype_gall_parasitoid_network))*100,
                              genotypes = rownames(genotype_gall_parasitoid_network))
geno.contrib.df <- mutate(geno.contrib.df, geno.ord = reorder(genotypes, percent.contrib))
mean.percent.contrib <- with(geno.contrib.df, mean(percent.contrib))

hist(log(geno.contrib.df$contrib)) # log-normal distribution?

ggplot(geno.contrib.df, aes(x = geno.ord, y = percent.contrib)) +
  geom_bar(stat = "identity") +
  xlab("Genotype") +
  ylab("Percent Contribution to Network Structure") +
  coord_flip() +
  theme_classic() +
  geom_hline(yintercept = mean.percent.contrib, linetype = "dashed")

# Gall-parasitoid interactions
gall.ptoid.contrib.df <- data.frame(contrib = colSums(genotype_gall_parasitoid_network),
                              percent.contrib = colSums(genotype_gall_parasitoid_network)/sum(colSums(genotype_gall_parasitoid_network))*100,
                              interactions = colnames(genotype_gall_parasitoid_network))
gall.ptoid.contrib.df <- mutate(gall.ptoid.contrib.df, interaction.ord = reorder(interactions, percent.contrib))
mean.gp.percent.contrib <- with(gall.ptoid.contrib.df, mean(percent.contrib))

hist(log(gall.ptoid.contrib.df$contrib)) # log-normal distribution?

ggplot(gall.ptoid.contrib.df, aes(x = interaction.ord, y = percent.contrib)) +
  geom_bar(stat = "identity") +
  xlab("Interactions") +
  ylab("Percent Contribution to Network Structure") +
  coord_flip() +
  theme_classic() +
  geom_hline(yintercept = mean.gp.percent.contrib, linetype = "dashed")

# betadiversity of genotype-gall-parasitoid network
beta.div(genotype_gall_parasitoid_network, method = c("hellinger")) # vLG_Platy, vLG_Mesopol, and rG_Tory contribute the most to beta-diversity.

geno.order <- arrange(genotypes, Vertex)

adonis(genotype_gall_parasitoid_network ~ factor(geno.order$Module.ID), method = "horn") # looks the same as the results 
anova(betadisper(vegdist(genotype_gall_parasitoid_network, method = "horn"), geno.order$Module.ID))

plot(capscale(genotype_gall_parasitoid_network, method = "horn"))

# Notes for Null Model Analysis:
# r2dtable = Patefield algorithm (preserves marginal totals of network)
# shuffle.web = preserves connectance of observed web
# swap.web = preserves marginal total (r2dtable) and connectance of observed web (shuffle.web)
# note the qualitative connectance, and not weighted connectance is preserved by these null models.
# According to Dormann et al. 2009, using the different null models may give us insight to different processes influencing the structure of the network.

### Visualize how the null models randomize the network and alter network structure

# create 1 rep of the different null models
null_r2dtable <- nullmodel(genotype_gall_parasitoid_network, N = 1, method = "r2dtable")
null_shuffle.web <- nullmodel(genotype_gall_parasitoid_network, N = 1, method = "shuffle.web")
null_swap.web <- nullmodel(genotype_gall_parasitoid_network, N = 1, method = "swap.web")

# visualize how these models may alter the compartmentalization of the network
visweb(genotype_gall_parasitoid_network, type = "diagonal")
visweb(null_r2dtable[[1]], type = "diagonal")
visweb(null_shuffle.web[[1]], type = "diagonal")
visweb(null_swap.web[[1]], type = "diagonal")

# visuazlize how these null models may alter the nestedness of the network
visweb(genotype_gall_parasitoid_network, type = "nested")
visweb(null_r2dtable[[1]], type = "nested")
visweb(null_shuffle.web[[1]], type = "nested")
visweb(null_swap.web[[1]], type = "nested")

### Null Model Analysis

# Nestedness: weighted NODF (Almeida-Neto & Ulrich, 2010, Environmental Modelling Software)
qual_genotype_gall_parasitoid_network <- genotype_gall_parasitoid_network
qual_genotype_gall_parasitoid_network[qual_genotype_gall_parasitoid_network > 0] <- 1
nestednodf(qual_genotype_gall_parasitoid_network, order = T, weighted = F)
WNODF_info <- nestednodf(genotype_gall_parasitoid_network, order = T, weighted = T) # 23.69687
WNODF = WNODF_info$statistic[3]
N_null_webs = 1000

z_WNODF_r2dtable <- null_model_analysis_WNODF_QuaBiMo(web = genotype_gall_parasitoid_network, 
                                                      observed_value = WNODF,
                                                      type = "WNODF", 
                                                      N_null_webs = N_null_webs, 
                                                      null_model = "r2dtable")

z_WNODF_shuffle.web <- null_model_analysis_WNODF_QuaBiMo(web = genotype_gall_parasitoid_network,
                                                         observed_value = WNODF,
                                                         type = "WNODF", 
                                                         N_null_webs = N_null_webs,
                                                         null_model = "shuffle.web")

z_WNODF_swap.web <- null_model_analysis_WNODF_QuaBiMo(web = genotype_gall_parasitoid_network, 
                                                      observed_value = WNODF,
                                                      type = "WNODF",
                                                      N_null_webs = N_null_webs, 
                                                      null_model = "swap.web")

null.models <- c("r2dtable",
                 "shuffle.web",
                 "swap.web")

results_WNODF <- data.frame(null.models = null.models, 
                            WNODF_z_scores = as.vector(c(z_WNODF_r2dtable,
                                                         z_WNODF_shuffle.web,
                                                         z_WNODF_swap.web)))

write.csv(results_WNODF, "~/Documents/Genotype_Networks/data/results_WNODF_genotype_gall_parasitoid_network.csv")

t <- computeModules(t(genotype_gall_parasitoid_network))
tt <- czvalues(t)
plot(tt[[1]], tt[[2]], xlab = "c", ylab = "z")

# Modularity: QuaBiMo algorithm (Dormann & Strauss 2014, Methods in Ecology and Evolution)
# start time 6:50 pm on Friday, July 18, 2014
# Note: when I removed rG_Mesopol, the algorithm converges on the same set of interacting genotypes and gall-parasitoid interactions. When I remove vLG_Mymarid, Genotype V is placed into the A,J,O,K,Z,S,N,H
# Using the default number of steps, but searching for the best partition at least 100 times seemed to result in a stable modularity configuration. Boosting the number of steps to 1e10 does not seem to improve this detection, and appears to actually be worse than the method I'm using.
remove_weak_interactions_df <- select(genotype_gall_parasitoid_network, -vLG_Mymarid)
best_partitions_weak_removed <- best_QuaBiMo(remove_weak_interactions_df, QuaBiMo_reps = 100)
plotModuleWeb(best_partitions_weak_removed$best_module_info[[1]])

t <- best_QuaBiMo(select(genotype_gall_parasitoid_network, starts_with("vLG")), QuaBiMo_reps = 10)
plotModuleWeb(t$best_module_info[[1]])
t$best_observed_value

t2 <- null_model_analysis_WNODF_QuaBiMo(web = select(genotype_gall_parasitoid_network, starts_with("vLG")), observed_value = 0.305034, type = "QuaBiMo", N_null_webs = 100, null_model = "swap.web")

best_partitions <- best_QuaBiMo(genotype_gall_parasitoid_network, QuaBiMo_reps = 100) # 0.327593
plotModuleWeb(best_partitions$best_module_info[[1]])
czvalues(best_partitions$best_module_info[[1]], weighted = TRUE)
QuaBiMo_Q = best_partitions$best_observed_value
N_null_webs = 1000

z_QuaBiMo_r2dtable <- null_model_analysis_WNODF_QuaBiMo(web = genotype_gall_parasitoid_network, 
                                                      observed_value = QuaBiMo_Q,
                                                      type = "QuaBiMo", 
                                                      N_null_webs = N_null_webs,
                                                      null_model = "r2dtable")

z_QuaBiMo_shuffle.web <- null_model_analysis_WNODF_QuaBiMo(web = genotype_gall_parasitoid_network,
                                                         observed_value = QuaBiMo_Q,
                                                         type = "QuaBiMo",
                                                         N_null_webs = N_null_webs,
                                                         null_model = "shuffle.web")

z_QuaBiMo_swap.web <- null_model_analysis_WNODF_QuaBiMo(web = genotype_gall_parasitoid_network, 
                                                        observed_value = QuaBiMo_Q, 
                                                        type = "QuaBiMo",
                                                        N_null_webs = N_null_webs, 
                                                        null_model = "swap.web")
z_QuaBiMo_swap.web # last time: z-score = 2.4138801, mean null values = 0.28276243, SD null values = 0.01857199, 1000 reps of null model.
results_QuaBiMo <- data.frame(null.models = null.models, 
                              QuaBiMo_z_scores = as.vector(c(z_QuaBiMo_r2dtable,
                                                            z_QuaBiMo_shuffle.web,
                                                            z_QuaBiMo_swap.web)))

write.csv(results_QuaBiMo, "~/Documents/Genotype_Networks/data/results_QuaBiMo_genotype_gall_parasitoid_network.csv")


### Matrix visualization of network
visweb(genotype_gall_parasitoid_network, type = "diagonal")
source('~/Documents/miscellaneous_R/ggplot_themes.R')

module.info <- read.csv('~/Documents/Genotype_Networks/data/module_info_genotype_gall_parasitoid_network.csv')
module.info <- arrange(module.info, Trophic, Module.ID, Vertex) # reorder for easy plotting
genotypes <- filter(module.info, Trophic == "Genotype")
gall.parasitoid <- filter(module.info, Trophic == "Gall-Parasitoid")
#gall.parasitoid$Vertex.replace <- c("Rabdophaga Bud-Torymus",
 #                                   "Iteomyia-Eulophid",
  #                                  "Iteomyia-Mymarid",
   #                                 "Iteomyia-Torymus",
    #                                "Rabdophaga Stem-Platygaster",
     #                               "Rabdophaga Bud-Mesopolobus",
      #                              "Iteomyia-Platygaster",
       #                             "Rabdophaga Bud-Eulophid",
        #                            "Rabdophaga Bud-Lestodiplosis",
         #                           "Rabdophaga Bud-Platygaster",
          #                          "Cecidomyiid A-Torymus",
           #                         "Pontania-Eurytoma",
            #                        "Pontania-Lathrostizus",
             #                       "Iteomyia-Mesopolobus")


# reorder based on module ID
web <- genotype_gall_parasitoid_network
#colnames(web) <- c("Cecidomyiid A-Torymus",
 #                  "Rabdophaga Bud-Eulophid",
  #                 "Rabdophaga Bud-Lestodiplosis",
   #                "Rabdophaga Bud-Mesopolobus",
    #               "Rabdophaga Bud-Platygaster",
     #              "Rabdophaga Bud-Torymus",
      #             "Pontania-Eurytoma",
       #            "Pontania-Lathrostizus",
         #          "Rabdophaga Stem-Platygaster",
          #         "Iteomyia-Eulophid",
           #        "Iteomyia-Mesopolobus",
            #       "Iteomyia-Mymarid",
             #      "Iteomyia-Platygaster",
              #     "Iteomyia-Torymus")

web.modular.order <- web[as.character(genotypes$Vertex), # reorder within web
           as.character(gall.parasitoid$Vertex)] # reorder within web
colSums(web.modular.order)
rowSums(web.modular.order)

web.modular_genotypes <- mutate(web.modular.order, Genotypes = factor(rownames(web.modular.order), levels = rev(rownames(web.modular.order)), ordered = TRUE))
levels(web.modular_genotypes$Genotypes)[7] <- "C" # replace * for plotting aesthetics
web.modular_genotypes_melt <- melt(web.modular_genotypes)

# alphabetical order
rownames(web)[1] <- "C" # replace * for plotting aesthetics
web <- web[c(2,3,1,4:25), ]
colSums(web)
rowSums(web)

web.alpha_melt <- mutate(web, Genotypes = factor(rownames(web), levels = rev(rownames(web)), ordered = TRUE))
web.alpha_melt <- melt(web.alpha_melt)


ggplot(web.alpha_melt, aes(y = Genotypes, x = variable, fill = value)) + 
  geom_tile(color = "black") + 
  scale_fill_gradient(low = "white", high = "steelblue", guide = "none") + 
  geom_text(label = web.alpha_melt$value, 
            family = "Verdana") +
  theme_heatmap + 
  theme(axis.ticks = element_blank(),
        axis.text.y = element_text(size = 18,
                                   family = "Verdana")) +
  xlab("") +
  ylab("")
ggsave("~/Documents/Genotype_Networks/heatmap_genotype_gall_parasitoid_alpha_web.png", width = 8, height = 8, units = "in", dpi = 300)



ggplot(web.modular_genotypes_melt, aes(y = Genotypes, x = variable, fill = value)) + 
  geom_tile(color = "black") + 
  scale_fill_gradient(low = "white", high = "steelblue", guide = "none") + 
  geom_text(label = web.modular_genotypes_melt$value, 
            family = "Verdana") +
  theme_heatmap + 
  theme(axis.ticks = element_blank(),
        axis.text.y = element_text(size = 18,
                                   family = "Verdana")) +
  xlab("") +
  ylab("")
ggsave("~/Documents/Genotype_Networks/heatmap_genotype_gall_parasitoid_modular_web.png", width = 8, height = 8, units = "in", dpi = 300)

web.swap_null <- swap.web(N = 1, web = web)
rownames(web.swap_null[[1]]) <- rownames(web)
colnames(web.swap_null[[1]]) <- colnames(web)
web.swap_null <- data.frame(web.swap_null[[1]])

web.swap_melt <- mutate(web.swap_null, Genotypes = factor(rownames(web), levels = rev(rownames(web)), ordered = TRUE))
web.swap_melt <- melt(web.swap_melt)

ggplot(web.swap_melt, aes(y = Genotypes, x = variable, fill = value)) + 
  geom_tile(color = "black") + 
  scale_fill_gradient(low = "white", high = "steelblue", guide = "none") + 
  geom_text(label = web.swap_melt$value, 
            family = "Verdana") +
  theme_heatmap + 
  theme(axis.ticks = element_blank(),
        axis.text.y = element_text(size = 18,
                                   family = "Verdana")) +
  xlab("") +
  ylab("")
ggsave("~/Documents/Genotype_Networks/heatmap_genotype_gall_parasitoid_swap_alpha_web.png", width = 8, height = 8, units = "in", dpi = 300)

# cca ordering of null swap web. Mimics modularity ordering for visualization purposes
web.swap_null_modules <- computeModules(web.swap_null) # Q = 0.243
printoutModuleInformation(web.swap_null_modules)

swap.module.geno.order <- c("B","F","G","I","P","V","W",
                            "D","J","L","N","O","R","S","Y",
                            "A","E","M","Q","T","X","Z",
                            "C","H","K")
swap.module.inter.order <- c("aSG_Tory","rG_Platy","vLG_Mesopol","vLG_Mymarid",
                             "rG_Eulo","rG_Tory","rsLG_Eury","vLG_Eulo","vLG_Tory",
                             "rG_Mesopol","rsLG_Lathro","vLG_Platy",
                             "rG_Lestodip","SG_Platy")
#ca <- cca(web.swap_null)
web.swap_null_modular <- web.swap_null[swap.module.geno.order,
                                       swap.module.inter.order] # reorder null swap web based on modularity algorithm
#<- web.swap_null[order(summary(ca)$sites[, 1], decreasing = TRUE), 
          # order(summary(ca)$species[, 1], decreasing = TRUE)]

web.swap.module_melt <- mutate(web.swap_null_modular, Genotypes = factor(rownames(web.swap_null_modular), levels = rev(rownames(web.swap_null_modular)), ordered = TRUE))
web.swap.module_melt <- melt(web.swap.module_melt)

ggplot(web.swap.module_melt, aes(y = Genotypes, x = variable, fill = value)) + 
  geom_tile(color = "black") + 
  scale_fill_gradient(low = "white", high = "steelblue", guide = "none") + 
  geom_text(label = web.swap.module_melt$value, 
            family = "Verdana") +
  theme_heatmap + 
  theme(axis.ticks = element_blank(),
        axis.text.y = element_text(size = 18,
                                   family = "Verdana")) +
  xlab("") +
  ylab("")
ggsave("~/Documents/Genotype_Networks/heatmap_genotype_gall_parasitoid_swap.module_alpha_web.png", width = 8, height = 8, units = "in", dpi = 300)

# testing out new shuffle web null model. I guess this null model will not preserve the species richness of the web, and therefore not connectance, which is why they fill in the diagonals first, although this may artificially enhance the modularity of the web. This new model though may be relevant for identifying the contribution of the number of links vs. the number of species in contributing to modularity.
# Another useful model would be to be to assess the contribution of interaction distributions would be to calculate QuaBiMo on presence/absence data of the original matrix formulation.
test <- as.matrix(genotype_gall_parasitoid_network)
test.v <- as.vector(test)
test.v.new <- sample(test.v, size = length(test.v))
sum(test.v.new)
test.v.new.matrix <- matrix(test.v.new, ncol = ncol(test), dimnames = dimnames(genotype_gall_parasitoid_network))
test.v.new.matrix.Q <- computeModules(test.v.new.matrix)
plotModuleWeb(test.v.new.matrix.Q)
