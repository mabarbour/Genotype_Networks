#### Discriminant function analysis of modules

library(plyr)
library(dplyr)

genotype_gall_parasitoid_network <- read.csv('~/Documents/Genotype_Networks/data/genotype_gall_parasitoid_network_data.csv')
genotype_gall_parasitoid_network$Genotype <- genotype_gall_parasitoid_network$X
genotype_gall_parasitoid_network <- select(genotype_gall_parasitoid_network, -X)

module.info <- read.csv('~/Documents/Genotype_Networks/data/module_info_genotype_gall_parasitoid_network.csv')
module.info.genotype = module.info %>%
  select(Genotype = Vertex, Module.ID, Trophic) %>%
  filter(Trophic == "Genotype")

gall_density_genotype_summary_df <- read.csv('~/Documents/Genotype_Networks/data/gall_density_genotype_summary_df.csv')
gall_density_genotype_summary_df <- select(gall_density_genotype_summary_df, -X)

gall_size_genotype_summary_df <- read.csv('~/Documents/Genotype_Networks/data/gall_size_genotype_summary_df.csv')
gall_size_genotype_summary_df <- select(gall_size_genotype_summary_df, -X)

dfa_data <- join_all(list(genotype_gall_parasitoid_network, module.info.genotype, gall_density_genotype_summary_df, gall_size_genotype_summary_df))
dfa_data_other <- dfa_data

### Perhaps finding the regressions for some of the key players, rather than trying to see which ones predict the modules the best is the key...

cor.test(dfa_data_other$vLG_volume_mean, dfa_data_other$vLG_abund, na.action = na.omit)
plot(vLG_Platy ~ vLG_volume_mean, dfa_data_other)
vLG_Platy_glm <- glm(rsLG_Eury ~ rsLG_abund, dfa_data_other, family = "poisson")
summary(vLG_Platy_glm)
plot(vLG_Platy_glm)

library(visreg)
visreg(vLG_Platy_glm, scale = "response")

dfa_data <- dfa_data %>%
  select(Genotype, Module.ID, aSG_abund:rsLG_abund, vLG_volume_mean, rG_volume_mean) %>%
  mutate(Module.ID = as.factor(Module.ID))

plot(vLG_abund ~ Module.ID, dfa_data)
plot(vLG_volume_mean ~ Module.ID, dfa_data)
plot(rG_abund ~ Module.ID, dfa_data)
plot(rG_volume_mean ~ Module.ID, dfa_data)

module_dfa <- MASS::lda(Module.ID ~ aSG_abund + rG_abund + vLG_abund + SG_abund + rsLG_abund + vLG_volume_mean + rG_volume_mean, dfa_data, na.action = na.omit, CV = TRUE) # temporarily removing  + rG_volume_mean, because it would omit two replicate genotype
module_dfa

CV_test <- table(dfa_data$Module.ID[-c(5,17,20,21)], module_dfa$class)
diag(prop.table(CV_test,1))
sum(diag(prop.table(CV_test,1)))

module_dfa_values <- predict(module_dfa)

plot(module_dfa_values$x[ ,1], module_dfa_values$x[ ,2])
text(module_dfa_values$x[,1], module_dfa_values$x[,2], dfa_data$Module.ID, cex=0.7, pos=4, col="red")

MASS::ldahist(data = module_dfa_values$x[,1], factor(dfa_data$Module.ID)) # currently not working

library(randomForest)
module_rf <- randomForest(Module.ID ~ ., dfa_data[ ,-1], importance = TRUE, na.action = na.omit)
varImpPlot(module_rf)
