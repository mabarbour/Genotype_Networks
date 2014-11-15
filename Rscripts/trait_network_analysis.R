
##### load require libraries
library(psych)
library(plyr)
library(dplyr)

##### upload data
#genotype_gall_parasitoid_network <- read.csv('~/Documents/Genotype_Networks/data/genotype_gall_parasitoid_network_data.csv')
#rownames(genotype_gall_parasitoid_network) <- genotype_gall_parasitoid_network$X 
#genotype_gall_parasitoid_network <- select(genotype_gall_parasitoid_network, -X)
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

# vLG gall size
vLG_size <- read.csv("~/Documents/Genotype_Networks/data/vLG_gall_size_geno.df.csv")
vLG_size <- select(vLG_size, -X, Genotype, gall.height.vLG.mean = g.height.vLG.mean, large.vLG.count, small.vLG.count)


df <- join_all(list(glevel.net, glevel.traits, modules, vLG_size, glevel.gall.densities), by = "Genotype")
df <- mutate(df, vLG_ptoid_attack = vLG_Platy + vLG_Mesopol + vLG_Tory + vLG_Mymarid + vLG_Eulo,
             prop.small.galls = small.vLG.count/(large.vLG.count + small.vLG.count))

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
