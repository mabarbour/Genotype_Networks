## network mixed-effect models

## source in required data
source('~/Documents/Genotype_Networks/Rscripts/network_management_gall_level.R')

## source in required libraries
library(lme4)
library(visreg)
library(gamm4)

gall_mech.df

## vLG mechanisms df
vLG_mech.df <- filter(gall_mech.df, vLG_total > 0) %>%
  mutate(sc.gall.height = scale(gall.height, center = TRUE, scale = TRUE),
         sc.vLG_density = scale(vLG_density, center = TRUE, scale = TRUE),
         vLG_egg = vLG_Platy + vLG_Mymarid,
         vLG_ecto = vLG_Eulo + vLG_Mesopol + vLG_Tory)


## vLG_egg
vLG_egg.glmer <- glmer(cbind(vLG_egg, vLG_total - vLG_egg) ~ sc.vLG_density*sc.gall.height + 
                         (1 | Genotype.x/plant.position), data = vLG_mech.df,
                       family = binomial(link = "logit"))
summary(vLG_egg.glmer)
AIC(vLG_egg.glmer)

vLG_egg.glmer.vis <- glmer(cbind(vLG_egg, vLG_total - vLG_egg) ~ vLG_density*gall.height + 
                         (1 | plant.position), data = vLG_mech.df,
                       family = binomial(link = "logit"))
summary(vLG_egg.glmer.vis) # same overall picture as standardized model that includes Genotype as a random effect.
AIC(vLG_egg.glmer.vis)

## vLG_ecto
vLG_ecto.glmer <- glmer(cbind(vLG_ecto, vLG_total - vLG_ecto) ~ sc.vLG_density + sc.gall.height +
                         (1 | Genotype.x/plant.position), data = vLG_mech.df, 
                        family = binomial(link = "logit"))
summary(vLG_ecto.glmer)
AIC(vLG_ecto.glmer)

vLG_ecto.glmer.vis <- glmer(cbind(vLG_ecto, vLG_total - vLG_ecto) ~ vLG_density + gall.height + 
                             (1 | plant.position), data = vLG_mech.df,
                           family = binomial(link = "logit"))
summary(vLG_ecto.glmer.vis) # same overall picture as standardized model that includes Genotype as a random effect.
AIC(vLG_ecto.glmer.vis)

## rG mechanisms df
rG_mech.df <- filter(gall_mech.df, rG_total > 0) %>%
  mutate(sc.gall.height = scale(gall.height, center = TRUE, scale = TRUE),
         sc.rG_density = scale(rG_density, center = TRUE, scale = TRUE),
         rG_ecto = rG_Tory + rG_Mesopol + rG_Eulo)

## rG_ecto
rG_ecto.glmer <- glmer(cbind(rG_ecto, rG_total - rG_ecto) ~ 1 + 
                          (1 | Genotype.x/plant.position), data = rG_mech.df, family = binomial(link = "logit"))
summary(rG_ecto.glmer)
AIC(rG_ecto.glmer)

## plots for manuscript. Consider theta = 45
grey.scale <- colorRampPalette(colors = c("white","grey","black"))

#pdf(file="~/Documents/Genotype_Networks/iteomyia egg parasitoid mixed effect.png",width=11, height=8.5)
visreg2d(vLG_egg.glmer.vis, y = "vLG_density", x = "gall.height", 
         scale = "response", color.palette = grey.scale, 
         zlim = c(0,1), main = "Probability of Iteomyia - Egg parasitoid interaction",
         xlab = "Gall diameter (mm)", ylab = "Gall density (count per 100 shoots)")
#dev.off()
visreg2d(vLG_ecto.glmer.vis, y = "vLG_density", x = "gall.height", 
         scale = "response", color.palette = grey.scale, 
         zlim = c(0,1), main = "Probability of Iteomyia - Larval parasitoid interaction",
         xlab = "Gall diameter (mm)", ylab = "Gall density (count per 100 shoots)")


