# Genotype level trait-gall and gall-ptoid associations
library(dplyr)
library(psych)
library(ggplot2)
library(visreg)

# load data
plant.traits <- read.csv('~/Documents/Genotype_Networks/data/plant.trait.galls.2011.genotype.df.csv')
plant.traits <- select(plant.traits, -X)

genotype_level_df <- read.csv('~/Documents/Genotype_Networks/data/genotype_level.df.csv')
genotype_level_df <- select(genotype_level_df, -X)

module.info <- read.csv('~/Documents/Genotype_Networks/data/module_info_genotype_gall_parasitoid_network.csv')
module.info <- arrange(module.info, Trophic, Module.ID, Vertex) # reorder for easy plotting
genotypes <- filter(module.info, Trophic == "Genotype")
genotypes <- select(genotypes, Genotype = Vertex, Module.ID)

#galldens <- read.csv("~/Documents/Genotype_Networks/data/gall_density_genotype_summary_df.csv")
#galldens <- select(galldens, -X)

vLG_gallsize <- read.csv("~/Documents/Genotype_Networks/data/vLG_gall_size_geno.df.csv")
vLG_gallsize <- select(vLG_gallsize, -X)
#gallsize <- read.csv("~/Documents/Genotype_Networks/data/gall_size_genotype_summary_df.csv")
#gallsize <- select(gallsize, -X)

#galldens.size <- left_join(galldens, gallsize)
#galldens.size.noU <- filter(galldens.size, Genotype != "U")

geno.df <- plyr::join_all(list(genotype_level_df, plant.traits, genotypes, vLG_gallsize), type = "left")

# vLG_Platy
plot(vLG_Platy ~ vLG_density, geno.df)
plot(vLG_Platy ~ g.height.vLG.mean, geno.df)
plot(vLG_Platy ~ prop.small.vLG, geno.df)

plot(vLG_Mesopol ~ vLG_density, geno.df)
plot(vLG_Mesopol ~ g.height.vLG.mean, geno.df)
plot(vLG_Mesopol ~ prop.small.vLG, geno.df)

small.vLG.density <- with(geno.df, small.vLG.count/shootEst.no18)

plot(vLG_Platy_dens ~ small.vLG.density)
prop.small.vLG <- with(geno.df, small.vLG.count/(large.vLG.count + small.vLG.count))
vLG_Platy_dens <- with(geno.df, vLG_Platy/shootEst.no18)
vLG_Mesopol_dens <- with(geno.df, vLG_Mesopol/shootEst.no18)

theme <- theme_classic() + theme(axis.text = element_text(size = 20))

ggplot(geno.df, aes(x = vLG_density*1000, y = vLG_Platy_dens*1000)) + geom_point(size = 5, shape = 1) + stat_smooth(method = "lm", size = 2) + theme + xlab("Iteomyia individuals per 1000 shoots") + ylab("Iteomyia-Platygaster interactions per 1000 shoots")
ggsave("~/Documents/Genotype_Networks/vLG_Platy_regression.png", height = 9, width = 9, units = "in", dpi = 300)

ggplot(geno.df, aes(x = vLG_density*1000, y = vLG_Mesopol_dens*1000)) + geom_point(size = 5, shape = 1) + stat_smooth(method = "lm", size = 2) + theme + xlab("Iteomyia individuals per 1000 shoots") + ylab("Iteomyia-Mesopolobus interactions per 1000 shoots")
ggsave("~/Documents/Genotype_Networks/vLG_Mesopol_regression.png", height = 9, width = 9, units = "in", dpi = 300)

summary(lm(vLG_Platy_dens ~ small.vLG.density)) # note that small vLG density doesn't exactly correspond with gall density...because it uses the gall as the replicate rather than the individual...
summary(lm(vLG_Platy_dens ~ vLG_density, geno.df))
summary(lm(vLG_Mesopol_dens ~ vLG_density, geno.df))
summary(lm(vLG_Mesopol_dens ~ small.vLG.density))

with(geno.df, cbind(small.vLG.count, large.vLG.count, vLG_abund))
#### Exploratory data analysis
vLG_density <- with(geno.df, vLG_abund/shootEst.no18)
rG_density <- with(geno.df, rG_abund/shootEst.no18)
rG_vLG_dens.cor <- corr.test(y = cbind(vLG_density, rG_density), x = select(geno.df, Total_Area:specific_leaf_area))

plot(rG_density ~ D_mean_smoothed, geno.df)

summary(lm(sqrt(vLG_abund/shootEst.no18) ~ D_mean_smoothed, geno.df))

summary(MASS::glm.nb(rG_abund ~ offset(log(shootEst.no18)) + D_mean_smoothed, geno.df))

# vLG_Platy model
vLG_density <- with(geno.df, vLG_abund/shootEst.no18)
vLG_density_1000 <- vLG_density*1000
vLG_Platy_density <- round(with(geno.df, vLG_Platy/shootEst.no18*1000)) # rounded for use in glm model

plot(vLG_Platy_density ~ vLG_density_1000)
text(y = vLG_Platy_density, x = vLG_density_1000, labels = geno.df$Genotype)

ggplot(data = geno.df, aes(x = vLG_density_1000, y = vLG_Platy_density)) + geom_point() + stat_smooth(method = "glm", family = "quasipoisson")

vLG_Platy_glm <- MASS::glm.nb(vLG_Platy_density ~ vLG_density_1000 + vLG_volume_mean, geno.df) # still significant after remvoing genotype X
summary(vLG_Platy_glm)
plot(vLG_Platy_glm)

visreg(vLG_Platy_glm, scale = "response")

vLG_Platy_glm <- MASS::glm.nb(vLG_Platy_density ~ sqrt(vLG_density_1000) + sqrt(vLG_volume_mean), data = geno.df)
summary(vLG_Platy_glm)
plot(vLG_Platy_glm)
visreg(vLG_Platy_glm, scale = "response")

plot(vLG_Platy_density ~ vLG_volume_mean, geno.df)

test <- visreg(vLG_Platy_glm, "vLG_density_1000", scale = "response")
test$vLG_density_1000$x$x
test$vLG_density_1000$y$fit

# 2011-2012 vLG correlations
plot(sqrt(vLG_) ~ sqrt(vLG_2011_glm.nb_density), geno.df)
cor.test(log(geno.df$vLG_2011_glm.nb_density[-c(17,21)]), log(geno.df$vLG_abund[-c(17,21)]))

# 2012 vLG, rG correlations. Significant
plot(sqrt(vLG_abund) ~ sqrt(rG_abund), geno.df)
plot(geno.df$vLG_dens_modelfit ~ geno.df$rG_dens_modelfit)
cor.test(sqrt(geno.df$vLG_abund), sqrt(geno.df$rG_abund))
cor.test(sqrt(geno.df$vLG_dens_modelfit), sqrt(geno.df$rG_dens_modelfit))

rG_vLG.df <- geno.df %>%
  select(vLG_2011_glm.nb_density, vLG_dens_modelfit, rG_dens_modelfit, vLG_volume_mean, rG_volume_mean, rG_abund, vLG_abund) %>%
  mutate(log.vLG.2011 = log(vLG_2011_glm.nb_density+1), log.vLG.2012.model = log(vLG_dens_modelfit+1), log.vLG.abund = log(vLG_abund+1), log.vLG.vol = log(vLG_volume_mean), log.rG.model = log(rG_dens_modelfit+1), log.rG.abund = log(rG_abund+1),  log.rG.vol = log(rG_volume_mean) )

cors <- corr.test(x = select(plant.traits, Total_Area:flavanonOLES.PC1), y = rG_vLG.df, adjust = "none")

plot(log(vLG_volume_mean+1) ~ luteolin.der2__A320nm, geno.df)
