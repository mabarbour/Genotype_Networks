## This script analyzes the association between plant traits and gall-parasitoid networks.


plot(log(vLG.2011+1) ~ log(Total_Area) + C_N_imputed, tree_level_interaxn_all_plants_traits_size)
test <- with(tree_level_interaxn_all_plants_traits_size, vLG.2011/Total_Area)
visreg(lm(test ~ C_N_imputed, tree_level_interaxn_all_plants_traits_size))
summary(lm(vLG.2011 ~ poly(Total_Area,2) + C_N_imputed, tree_level_interaxn_all_plants_traits_size))

summary(MASS::glm.nb(vLG.2011 ~ offset(log(Total_Area)) + C_N_imputed + Height, tree_level_interaxn_all_plants_traits_size))

test <- MASS::glm.nb(vLG_abund ~ offset(log(shootEst.no18)) + C_N_imputed + Height, tree_level_interaxn_all_plants_traits_size)
summary(test)
library(AER)
dispersiontest(test)

test.df <- filter(tree_level_interaxn_all_plants_traits_size, vLG_abund > 0)
summary(lm(log(vLG_abund) ~ log(shootEst.no18) + C_N_imputed, test.df))
summary(lm(log(vLG_abund) ~ offset(log(shootEst.no18)) + C_N_imputed, test.df))
summary(lm(log(vLG_abund/shootEst.no18) ~ C_N_imputed, test.df))

shoothigh <- which(tree_level_interaxn_all_plants_traits_size$shootEst.no18 > 800)
## upload source code for managing network data at the tree-level
source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')
source('~/Documents/miscellaneous_R/length2_function.R')

## upload required libraries
library(visreg) # regression visualization

## get genotype-level interaction data. I accounted for density by dividing by the number of estimated shoots, because that was the most straightforward way to account for differences in sampling among trees.
tree_level_interaxn_density <- select(tree_level_interaxn_all_plants, aSG_aSG.larv:gall_survive_abund)/tree_level_interaxn_all_plants$shootEst.no18*100 # scale to a per 100 shoot basis

geno_level_interaxn_density <- aggregate(tree_level_interaxn_density, 
                                 by = list(Genotype = tree_level_interaxn_all_plants$Genotype), 
                                 FUN = mean) %>%
  tbl_df()


## load plant trait data. Consider downloading the data straight from dryad once the repository becomes open.
tree_level_traits <- read.csv("~/Documents/Genotype_Networks/data/plant.trait.galls.2011.tree.df.csv")
tree_level_traits <- select(tree_level_traits, -X, Genotype, plant.position, vLG.2011 = vLG, rsLG.2011 = rsLG, Total_Area:flavanonOLES.PC1)
tree_level_traits$plant.position <- factor(tree_level_traits$plant.position)

## get genotype-level trait data.
geno_level_traits <- tree_level_traits %>%
  group_by(Genotype) %>%
  select(Total_Area:flavanonOLES.PC1) %>%
  summarise_each(funs(mean.na.rm = mean(., na.rm = TRUE))) 

## manage gall trait data. The code below resolves everything down to the gall level and not the larva level, since there may be multiple larva within a gall.
gall.size.df <- tree_interaxns_filter_add %>%
  group_by(Genotype, plant.position, gall.id.nest, gall.sp) %>% # this makes sure I don't double count galls with multiple larvae, and is also why I take the "average" below, which is just the same number.
  summarise(gall.height = mean(g.height))

gall.contents.df <- tree_interaxns_filter_add %>%
  dcast(gall.id.nest + gall.sp ~ gall_contents_collapse, sum) %>%
  tbl_df()

gall.size.contents.df <- left_join(gall.size.df, gall.contents.df)

tree_level_gall.size.df <- gall.size.contents.df %>%
  group_by(Genotype, plant.position, gall.sp) %>%
  summarise(gall.height.mean = mean(gall.height, na.rm = TRUE), gall.count = n())

tree_level_gall.size.cast <- dcast(tree_level_gall.size.df, Genotype + plant.position ~ gall.sp, 
      value.var = "gall.height.mean") %>% # may need to alter
  select(Genotype, plant.position, 
         aSG.height.mean = aSG, rG.height.mean = rG, rsLG.height.mean = rsLG, 
         SG.height.mean = SG, vLG.height.mean = vLG) %>%
  tbl_df()

tree_level_gall.count.cast <- dcast(tree_level_gall.size.df, Genotype + plant.position ~ gall.sp,
                                   value.var = "gall.count") %>% # may need to alter
  select(Genotype, plant.position, 
         aSG.gall.count = aSG, rG.gall.count = rG, rsLG.gall.count = rsLG, 
         SG.gall.count = SG, vLG.gall.count = vLG) %>%
  tbl_df()

tree_level_gall.size.count.join <- left_join(tree_level_gall.size.cast, tree_level_gall.count.cast)

geno_level_gall.size.count.join <- select(tree_level_gall.size.count.join, 
                                          Genotype:vLG.height.mean) %>% # counts per tree not needed
  group_by(Genotype) %>% 
  select(-plant.position) %>%
  summarise_each(funs(mean.na.rm = mean(., na.rm = TRUE), 
                      length2.na.rm = length2(., na.rm = TRUE))) # length2 enables me to calculate the number of plant replicates that were used to estimate mean gall size.

## merge interaction and tree-trait data
tree_level_interaxn_all_plants_traits_size <- join_all(list(tree_level_interaxn_all_plants, tree_level_traits, tree_level_gall.size.count.join), by = "plant.position") %>%
  tbl_df()

## merge interaction and genotype-trait data
genotype_level_interaxn_traits_size <- join_all(list(geno_level_interaxn_density, geno_level_traits, geno_level_gall.size.count.join), by = "Genotype") %>%
  tbl_df()

#######################
# there is at least one major tree outlier, and it is the one with the lowest sampled number of shoots. Therefore, I think I'm justified in removing it. Also an outlier for vLG.
test2 <- lm(sqrt(rG_abund) ~ tree_level_interaxn_all_plants_traits_size$Genotype, tree_level_interaxn_density)
plot(test2)
aggregate(tree_level_interaxn_density$rG_abund, by = list(tree_level_interaxn_all_plants_traits_size$Genotype), sd)
library(vegan)
## RDA analysis at tree level
tree_network_traits <- select(tree_level_interaxn_all_plants_traits_size, 
                  aSG_Tory:rG_Platy, rG_Tory:rsLG_Lathro, SG_Platy, vLG_Eulo:vLG_Tory,
                  Total_Area, sal_tannin.PC1:flavanonOLES.PC1, water_content, C_N_imputed, shootEst.no18) %>%
  na.omit()

interactions <- select(tree_network_traits, aSG_Tory:rG_Platy, rG_Tory:rsLG_Lathro, SG_Platy, vLG_Eulo:vLG_Tory) %>%
  names()
rda_network <- rda(log(tree_network_traits[ ,interactions] + 1) ~ Condition(log(shootEst.no18)) + C_N_imputed, tree_network_traits)
RsquareAdj(rda_network)
anova(rda_network)
plot(rda_network, display = c("sp", "bp"))

## correlations at genotype level
library(psych)
geno.corrs <- corr.test(x = select(genotype_level_interaxn_traits_size, 
                                   vLG_Platy, vLG_Mesopol, vLG_Tory, rG_Tory), 
                        y = select(genotype_level_interaxn_traits_size,
                                   Total_Area:specific_leaf_area), adjust = "none")
geno.corrs
plot(sqrt(rG_Tory) ~ sqrt(C_N_ratio), genotype_level_interaxn_traits_size)

## tests
vLG_Platy.lm <- lm(sqrt(vLG_Platy) ~ vLG_abund + sqrt(vLG.height.mean_mean.na.rm), genotype_level_interaxn_traits_size)
summary(vLG_Platy.lm)
visreg(vLG_Platy.lm)
size <- tree_level_interaxn_all_plants_traits_size$vLG.height.mean
vLG_Platy.tree.lm <- lm(vLG_Platy ~ size, tree_level_interaxn_density)
summary(vLG_Platy.tree.lm)
visreg(vLG_Platy.tree.lm, scale = "response")

plot(size ~ tree_level_interaxn_density$vLG_abund)
plot(vLG.height.mean_mean.na.rm ~ vLG_abund, genotype_level_interaxn_traits_size)

vLG_Tory.lm <- lm(vLG_Tory ~ vLG_abund + rG_abund, genotype_level_interaxn_traits_size)
summary(vLG_Tory.lm)
visreg(vLG_Tory.lm)

library(car)
vif(vLG_Tory.lm)

##### load require libraries
library(psych)
library(plyr)
library(dplyr)

##### upload data
genotype_gall_parasitoid_network <- read.csv('~/Documents/Genotype_Networks/data/genotype_gall_parasitoid_network_data.csv')
rownames(genotype_gall_parasitoid_network) <- genotype_gall_parasitoid_network$X 
genotype_gall_parasitoid_network <- select(genotype_gall_parasitoid_network, -X)

glevel.net <- read.csv('~/Documents/Genotype_Networks/data/genotype_level.df.csv')

glevel.traits <- read.csv("~/Documents/Genotype_Networks/data/plant.trait.galls.2011.genotype.df.csv")
glevel.traits <- select(glevel.traits, -X)

glevel.gall.densities <- read.csv("~/Documents/Genotype_Networks/data/gall_density_estimates_genotype_level.csv")
glevel.gall.densities <- select(glevel.gall.densities, -X)

# module ID
modules <- read.csv("~/Documents/Genotype_Networks/data/module_info_genotype_gall_parasitoid_network.csv")
modules <- modules %>%
  filter(Trophic == "Genotype") %>%
  mutate(Genotype = as.character(Vertex), Module.ID = factor(Module.ID)) %>%
  select(Genotype, Module.ID)

# contributions of genotypes to module identity
czvalues <- read.csv("~/Documents/Genotype_Networks/data/cz_values.df.csv")
czvalues.geno <- filter(czvalues, node.set == "lower")
colnames(czvalues.geno)[2] <- "Genotype"

# vLG gall size
vLG_size <- read.csv("~/Documents/Genotype_Networks/data/vLG_gall_size_geno.df.csv")
vLG_size <- select(vLG_size, -X, Genotype, gall.height.vLG.mean = g.height.vLG.mean, large.vLG.count, small.vLG.count)


df <- join_all(list(glevel.net, glevel.traits, modules, czvalues.geno, vLG_size, glevel.gall.densities), by = "Genotype")
df <- mutate(df, vLG_ptoid_attack = vLG_Platy + vLG_Mesopol + vLG_Tory + vLG_Mymarid + vLG_Eulo,
             prop.small.galls = small.vLG.count/(large.vLG.count + small.vLG.count))

##### Dissimilarity of traits and network. At first stab, the herbivore community dissimilarity (just in terms of community composition) is correlated with network dissimilarity (r = 0.19, P = 0.018, perm = 999). However, 
library(vegan)
ggp.net.dis <- vegdist(genotype_gall_parasitoid_network, method = "horn")
trait.dis <- vegdist(decostand(glevel.traits[-21,6:46], method = "standardize"), method = "euclidean")
herb.dis <- vegdist(df[-21,c(110, 115, 118, 120)], method = "horn")
mantel(xdis = herb.dis, ydis = ggp.net.dis, method = "pearson", permutations = 999) # herbivore community and network dissimilarity are correlated.
mantel(xdis = trait.dis, ydis = herb.dis, method = "pearson", permutations = 999) # no correlation between traits and herbivores, likely because these traits driving the dissimilarity among genotypes are not the ones that influence the herbivore community.

# should look at other genetic correlations
cor.test(df$gall.height.vLG.mean[-10], df$sal_tannin.PC1[-10]) # strong negative genetic correlation when genotype J is removed. It also as the smallest gall height by a lot, suggesting that it may be an outlier anyway with some other trait influencing gall size.
plot(df$gall.height.vLG.mean ~ df$sal_tannin.PC1, type = "n")
text(df$gall.height.vLG.mean ~ df$sal_tannin.PC1, label = df$Genotype)

##### Willow trait differences among modules. Will need to check residuals to be careful of variable weights.
summary(lm(log(C_N_ratio) ~ factor(Module.ID), df[-c(21,25),], weights = sqrt(z + abs(min(z, na.rm = TRUE)))))
plot(C_N_ratio ~ factor(Module.ID), df[-c(21,25),])

trait.names <- colnames(glevel.traits)[c(6:46, 50:55)]

trait.lms <- list()
for(i in 1:length(trait.names)){
  trait.lms[[i]] <- summary(lm(df[-c(21,25), trait.names[i]] ~ factor(df[-c(21,25),"Module.ID"]), df[-c(21,25), ], weights = (z + abs(min(z, na.rm = TRUE)))))
}

sig.traits <- trait.names[c(3,4,5,6,8,9,11,15,19,24,25,26,27,30,33,37,39,40,41,42,45)]
marg.sig.traits <- trait.names[c(1,12,14,17,20,28,31,36,47)]

plot(rda(df[-c(21,25), c("Total_Area","Density","Height","D_mean_smoothed")], scale = TRUE))
plot(rda(df[-c(21,25), c(trait.names[c(5:36,39:41)])], scale = TRUE), choices = c(11,12))

library(psych)

prcomp.leaf.quality <- prcomp(df[-c(21,25), c(trait.names[c(5:36,39:41)])], center = TRUE, scale = TRUE)
summary(prcomp.leaf.quality)
screeplot(prcomp.leaf.quality)
biplot(prcomp.leaf.quality, choices = c(7,8))
prcomp.leaf.quality$sdev^2 # first 11 components have eigenvalues > 1
prcomp.leaf.quality$rotation[ ,c(1,4,7,8,11)]

leaf.quality.PCs <- prcomp.leaf.quality$x[ ,1:11]

leaf.quality.PC.lms <- list()
for(i in 1:11){
  leaf.quality.PC.lms[[i]] <- summary(lm(leaf.quality.PCs[ ,i] ~ factor(df[-c(21,25),"Module.ID"]), df[-c(21,25), ], weights = (z + abs(min(z, na.rm = TRUE)))))
}
# PCs: 11,8,7,4(very strong),2(marginal),1

prcomp.architecture <- prcomp(df[-c(21,25), c("Total_Area","Density","Height","D_mean_smoothed")], center = TRUE, scale = TRUE)

architecture.PCs <- prcomp.architecture$x[ ,1:2]

architecture.PC.lms <- list()
for(i in 1:2){
  architecture.PC.lms[[i]] <- summary(lm(architecture.PCs[ ,i] ~ factor(df[-c(21,25),"Module.ID"]), df[-c(21,25), ], weights = (z + abs(min(z, na.rm = TRUE)))))
}
# PCs: 1 & 2


#### Correlations
df.numeric <- select(df, -X, -Genotype, -Module.ID)
df.corrs <- corr.test(df.numeric, use = "pairwise", method = "pearson", adjust = "none") # anything >0.38 is significantly correlated (no p-value adjustments)

plot(df$vLG_density_mean ~ df$C_N_ratio)
summary(lm(sqrt(vLG_density_mean) ~ C_N_ratio, df))

library(car)
scatterplotMatrix(df[ ,c(100,105,107,112,114)])
gall.corrs <- corr.test(df[ ,c(100,105,107,112,114)]) # vLG correlated with rG and rsLG. rsLG correlated with aSG.
library(vegan)
trait.pca <- rda(df[ ,52:88], scale = TRUE)
summary(trait.pca)
plot(trait.pca, choices = c(3,4))

library(ggplot2)
ggplot(df, aes(x = vLG_density_mean.noOutliers, y = vLG_Platy/vLG_abund)) +
  geom_point(size = df$vLG_abund/2, pch = 21, color = "grey") +
  stat_smooth(aes(weight = vLG_abund), method = "glm", family = "binomial") +
  theme_classic() 

ggplot(df, aes(x = prop.small.galls, y = vLG_Platy/vLG_abund)) +
  geom_point(size = df$vLG_abund/2, pch = 21, color = "grey") +
  stat_smooth(aes(weight = vLG_abund), method = "glm", family = "binomial") +
  theme_classic() 

# Best model of percent parasitism by Platy on vLG is vLG density. Gall size does appear to have an independent effect, but not when included with density. This is not due to collinearity either because I tried used the residuals of vLG size not explained by gall density as a predictor variable.
vLG_Platy.glm <- glm(vLG_Platy/vLG_abund ~ vLG_density_mean.noOutliers, data = df, family = "binomial", weights = vLG_abund) # robust to removing some of the datapoint.
summary(vLG_Platy.glm)
plot(vLG_Platy.glm)
vif(vLG_Platy.glm)

vLG.size.resids <- residuals(lm(gall.height.vLG.mean ~ vLG_density_mean.noOutliers, df))
cor(x = df$vLG_density_mean.noOutliers[-c(17,21)], y = df$gall.height.vLG.mean[-c(17,21)])

# vLG_Mesopol
ggplot(df, aes(x = vLG_density_mean.noOutliers, y = vLG_Mesopol/vLG_abund)) +
  geom_point(size = df$vLG_abund/2, pch = 21, color = "grey") +
  stat_smooth(aes(weight = vLG_abund), method = "glm", family = "binomial") +
  theme_classic() # not effect of density

ggplot(df, aes(x = prop.small.galls, y = vLG_Mesopol/vLG_abund)) +
  geom_point(size = df$vLG_abund/2, pch = 21, color = "grey") +
  stat_smooth(aes(weight = vLG_abund), method = "glm", family = "binomial") +
  theme_classic() 


vLG_Mesopol.glm <- glm(vLG_Mesopol/vLG_abund ~ gall.height.vLG.mean, data = df, family = "binomial", weights = vLG_abund) 
summary(vLG_Mesopol.glm)


vLG_attack.glm <- glm(vLG_ptoid_attack/vLG_abund ~ vLG_density_mean.noOutliers + gall.height.vLG.mean, data = df, family = "binomial", weights = vLG_abund) 
summary(vLG_attack.glm)

plot(gall.height.vLG.mean ~ Module.ID, df)
plot(vLG_density_mean.noOutliers ~ Module.ID, df)
plot(rG_density_mean.noOutliers ~ Module.ID, df)
summary(lm(vLG_density_mean.noOutliers ~ Module.ID, df))


###### NEED TO THINK ABOUT THE APPROPRIATE ANALYSES

plot(vLG_Platy/(vLG_Platy + vLG_vLG.pupa) ~ C_N_ratio, df)
t <- glm(cbind(vLG_Platy,vLG_vLG.pupa) ~ C_N_ratio + I(C_N_ratio^2), df, family = "binomial")
summary(t)
plot(t)

visreg(t, scale = "response")

plot(vLG_prop_ptoid_attack ~ vLG_dens, df)
t2 <- glm(vLG_ptoid_attack ~ offset(log(vLG_abund)) + vLG_dens, subset(df, vLG_abund > 0), family = "poisson")
summary(t2)
plot(t2)
visreg(t2, scale = "response")
