library(dplyr)
library(visreg)
library(ggplot2)

# upload and select plant trait data
plant.trait.data <- read.csv('~/Documents/Genotype_Willow_Community/datasets_&_Rscripts/arthropods_&_traits/complete.data.csv')
plant.trait.data <- tbl_df(plant.trait.data)

plant.trait.data <- plant.trait.data %>%
  select(Genotype, plant.position = Plant.Position, vLG, rsLG, SG, PG, Total_Area:flavanonOLES.PC1) # includes 2011 gall data (except for rG)

write.csv(plant.trait.data, "~/Documents/Genotype_Networks/data/plant.trait.galls.2011.tree.df.csv")

vLG_2011 <- MASS::glm.nb(vLG ~ offset(log(Total_Area)) + Genotype, plant.trait.data)
summary(vLG_2011)
vLG_2011_fit <- visreg(vLG_2011, "Genotype", scale = "response", cond = list(Total_Area = 2), ylim = c(0,60))
vLG_2011_density <- vLG_2011_fit$Genotype$y$fit

plant.trait.geno.df <- plant.trait.data %>%
  group_by(Genotype) %>%
  select(-Genotype, -plant.position) %>%
  summarise_each(funs(mean.na.rm = mean(., na.rm = TRUE))) %>%
  mutate(vLG_2011_glm.nb_density = vLG_2011_density)

write.csv(plant.trait.geno.df, '~/Documents/Genotype_Networks/data/plant.trait.galls.2011.genotype.df.csv')

#plant.trait.geno.noQ <- filter(plant.trait.geno.df, Genotype != "Q")

ggplot(plant.trait.geno.df, aes(x = log(C_N_ratio), y = log(vLG_2011_density+1))) + geom_point() + stat_smooth(method = "lm")
vLG_2011_CN <- lm(log(vLG_2011_density+1) ~ log(C_N_ratio), plant.trait.geno.df)
summary(vLG_2011_CN)
plot(vLG_2011_CN)

vLG_2011_height <- lm(vLG_2011_density ~ Height, plant.trait.geno.df)
summary(vLG_2011_height)
visreg(vLG_2011_height)
