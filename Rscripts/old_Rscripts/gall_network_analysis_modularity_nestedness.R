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

# node contribution
node.Q.contribs <- node.Q.contribution(web = genotype_gall_parasitoid_network, n_null_webs = 100, Null_model = "swap.web")
rownames(node.Q.contribs) <- colnames(genotype_gall_parasitoid_network)

write.csv(node.Q.contribs, "~/Documents/Genotype_Networks/data/node.Q.contribs.csv")

tree_level_interaxn_all_plants <- read.csv("~/Documents/Genotype_Networks/data/tree_level_interaxn_all_plants.csv")
tree_level_interaxn_all_plants <- mutate(tree_level_interaxn_all_plants, plant.position = factor(plant.position))

genotype_gall_web <- tree_level_interaxn_all_plants %>%
  group_by(Genotype) %>%
  select(vLG_abund, rG_abund, rsLG_abund, aSG_abund, SG_abund, shootEst.no18) %>%
  summarise_each(funs(sum)) %>%
  mutate(vLG_dens = round(vLG_abund/shootEst.no18*1000),
         rG_dens = round(rG_abund/shootEst.no18*1000),
         rsLG_dens = round(rsLG_abund/shootEst.no18*1000),
         aSG_dens = round(aSG_abund/shootEst.no18*1000),
         SG_dens = round(SG_abund/shootEst.no18*1000)) %>%
  select(Genotype, vLG_dens:SG_dens)
rownames(genotype_gall_web) <- genotype_gall_web[ ,1]

visweb(genotype_gall_web[ ,-1]) # not nested, even after null model
visweb(genotype_gall_web[ ,-1], type = "diagonal")
nestednodf(comm = genotype_gall_web[ ,-1], order = T, weighted = T) # 38.48
g.gall.web.best <- best_QuaBiMo(web = genotype_gall_web[ ,-1], QuaBiMo_reps = 10)
plotModuleWeb(g.gall.web.best$best_module_info[[1]]) # interesting, doesn't strongly correspond to modules identified from tri-trophic interaction.
czvalues(g.gall.web.best$best_module_info[[1]], weighted = TRUE)
g.gall.web.best$best_observed_value # 0.199335
g.gall.web.mod.null <- null_model_analysis_WNODF_QuaBiMo(web = genotype_gall_web[ ,-1], 
                                  observed_value = 0.199335,
                                  type = "QuaBiMo", 
                                  N_null_webs = 100, 
                                  null_model = "swap.web")
null_model_analysis_WNODF_QuaBiMo(web = genotype_gall_web[ ,-1], 
                                  observed_value = 38.48,
                                  type = "WNODF", 
                                  N_null_webs = 1000, 
                                  null_model = "swap.web")




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
QuaBiMo_Q = best_partitions$best_observed_value
N_null_webs = 1000

genotype.ptoid.network <- genotype_gall_parasitoid_network %>%
  mutate(Tory = aSG_Tory + rG_Tory + vLG_Tory,
         Platy = SG_Platy + vLG_Platy + rG_Platy,
         Mesopol = vLG_Mesopol + rG_Mesopol,
         Eulo = vLG_Eulo + rG_Eulo,
         Lathro = rsLG_Lathro,
         Eury = rsLG_Eury,
         Lestodip = rG_Lestodip,
         Mymarid = vLG_Mymarid) %>%
  select(Tory:Mymarid)
best_partitions.ptoids <- best_QuaBiMo(genotype.ptoid.network, QuaBiMo_reps = 100) # 0.327593
plotModuleWeb(best_partitions.ptoids$best_module_info[[1]])

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


### Examine c-z values of genotype-gall-parasitoid networks
source("~/Documents/Genotype_Networks/Rscripts/module.contribution.R")

# roles of herbivore-parasitoid nodes. Definitely a bit of a different perspective with looking at % module contribution.
cz.upper <- module.contribution(best_partitions$best_module_info[[1]], weighted = TRUE, level = "higher")
plot(cz.upper$module.percent.contribution ~ cz.upper$c, xlab = "c", ylab = "z", type = "n")
text(x = cz.upper$c, y = cz.upper$module.percent.contribution, labels = names(cz.upper$c)) 

#cz.upper.geno.ptoids <- czvalues(best_partitions.ptoids$best_module_info[[1]], weighted = TRUE, level = "higher")
#plot(cz.upper.geno.ptoids$z ~ cz.upper.geno.ptoids$c, xlab = "c", ylab = "z", type = "n")
#text(x = cz.upper.geno.ptoids$c, y = cz.upper.geno.ptoids$z, labels = names(cz.upper.geno.ptoids$c))

# qualitative role 
cz.upper.qual <- czvalues(best_partitions$best_module_info[[1]], weighted = FALSE, level = "higher")
plot(cz.upper.qual$z ~ cz.upper.qual$c, xlab = "c", ylab = "z", type = "n")
text(x = cz.upper.qual$c, y = cz.upper.qual$z, labels = names(cz.upper.qual$c)) 

# roles of genotypes. 
cz.lower <- module.contribution(best_partitions$best_module_info[[1]], weighted = TRUE, level = "lower")
plot(cz.lower$module.percent.contribution ~ cz.lower$c, xlab = "c", ylab = "z", type = "n")
text(x = cz.lower$c, y = cz.lower$module.percent.contribution, labels = names(cz.lower$c)) # I, K, and X appear to be module hubs. And possibly connector hubs, because they often have links to many other modules, but not a homogenous distribution (kinless hubs).

#cz.lower.geno.ptoids <- czvalues(best_partitions.ptoids$best_module_info[[1]], weighted = TRUE, level = "lower")
#plot(cz.lower.geno.ptoids$z ~ cz.lower.geno.ptoids$c, xlab = "c", ylab = "z", type = "n")
#text(x = cz.lower.geno.ptoids$c, y = cz.lower.geno.ptoids$z, labels = names(cz.lower.geno.ptoids$c))

### null model analysis of cz-values. May take a really long time depending on the number of null models.
source("~/Documents/Genotype_Networks/Rscripts/null_model_analysis_czvalues.R")
null_model_czvalues.df <- null_model_analysis_czvalues(web = genotype_gall_parasitoid_network, N_null_webs = 100, null_model = "swap.web")

null.high.c.means <- colMeans(null_model_czvalues.df$high.null_c.values, na.rm = TRUE)
null.high.c.sds <- apply(null_model_czvalues.df$high.null_c.values, 2, function(x) sd(x, na.rm = TRUE))
null.high.c.95quant <- apply(null_model_czvalues.df$high.null_c.values, 1, function(x) quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))

null.high.z.means <- colMeans(null_model_czvalues.df$high.null_z.values, na.rm = TRUE)
null.high.z.sds <- apply(null_model_czvalues.df$high.null_z.values, 2, function(x) sd(x, na.rm = TRUE))
null.high.z.95quant <- apply(null_model_czvalues.df$high.null_z.values, 1, function(x) quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))

null.low.c.means <- colMeans(null_model_czvalues.df$low.null_c.values, na.rm = TRUE)
null.low.c.sds <- apply(null_model_czvalues.df$low.null_c.values, 2, function(x) sd(x, na.rm = TRUE))
null.low.c.95quant <- apply(null_model_czvalues.df$low.null_c.values, 1, function(x) quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))

null.low.z.means <- colMeans(null_model_czvalues.df$low.null_z.values, na.rm = TRUE)
null.low.z.sds <- apply(null_model_czvalues.df$low.null_z.values, 2, function(x) sd(x, na.rm = TRUE))
null.low.z.95quant <- apply(null_model_czvalues.df$low.null_z.values, 1, function(x) quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))

null.cz.05.95.quantile.df <- data.frame(null.low.c.05.mean = mean(null.low.c.95quant[1, ]),
           null.low.c.95.mean = mean(null.low.c.95quant[2, ]),
           null.low.z.05.mean = mean(null.low.z.95quant[1, ]),
           null.low.z.95.mean = mean(null.low.z.95quant[2, ]),
           null.high.c.05.mean = mean(null.high.c.95quant[1, ]),
           null.high.c.95.mean = mean(null.high.c.95quant[2, ]),
           null.high.z.05.mean = mean(null.high.z.95quant[1, ]),
           null.high.z.95.mean = mean(null.high.z.95quant[2, ]),
           null.low.c.05.sd = sd(null.low.c.95quant[1, ]),
           null.low.c.95.sd = sd(null.low.c.95quant[2, ]),
           null.low.z.05.sd = sd(null.low.z.95quant[1, ]),
           null.low.z.95.sd = sd(null.low.z.95quant[2, ]),
           null.high.c.05.sd = sd(null.high.c.95quant[1, ]),
           null.high.c.95.sd = sd(null.high.c.95quant[2, ]),
           null.high.z.05.sd = sd(null.high.z.95quant[1, ]),
           null.high.z.95.sd = sd(null.high.z.95quant[2, ]))
           


# create and save dataframe of cz-values
cz.values.df <- data.frame(Nodes = c(names(cz.lower$c), names(cz.upper$c)), 
                           c.values = c(cz.lower$c, cz.upper$c), 
                           null.c.values.mean = c(null.low.c.means, null.high.c.means), 
                           null.c.values.sd = c(null.low.c.sds, null.high.c.sds),

                           z.values = c(cz.lower$module.percent.contribution, cz.upper$module.percent.contribution), 
                           null.z.values.mean = c(null.low.z.means, null.high.z.means), 
                           null.z.values.sd = c(null.low.z.sds, null.high.z.sds))
cz.values.df <- mutate(cz.values.df,
                       Zscore.c.values = (c.values - null.c.values.mean)/null.c.values.sd,
                       Zscore.z.values = (z.values - null.z.values.mean)/null.z.values.sd)

# current null model data based on 100 null models of r2dtable algorithm
write.csv(cz.values.df, "~/Documents/Genotype_Networks/data/cz.values.csv")
write.csv(null.cz.05.95.quantile.df, "~/Documents/Genotype_Networks/data/null.cz.05.95.quantile.df.csv")

# exploring cz.value data and cutoffs for assigning roles in network
cz.values.csv <- read.csv("~/Documents/Genotype_Networks/data/cz.values.csv")
cz.values.csv$trophic <- c(rep("Genotype",25), rep("gall-ptoid", 14))

hist(filter(cz.values.csv, trophic == "Genotype")$c.values)
hist(filter(cz.values.csv, trophic == "Genotype")$z.values, na.rm = T)

hist(filter(cz.values.csv, trophic == "gall-ptoid")$c.values)
hist(filter(cz.values.csv, trophic == "gall-ptoid")$z.values, na.rm = T)

ggplot(cz.values.csv, aes(x = c.values, y = z.values, group = trophic)) + 
  facet_grid( ~ trophic) + 
  geom_point() +
  #geom_text(aes(label = Nodes)) +
  theme_bw()

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
