## numerical responses of natural enemies

## source in required data
source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')

library(visreg)
library(lme4)
library(RLRsim)
library(vegan)
#library(pbkrtest)

#### Analyses at the tree-level

## make focal data frame for analysis
df <- tree_level_interaxn_all_plants_traits_size %>%
  select(aSG_Tory, rG_Eulo, rG_Lestodip, rG_Mesopol, rG_Platy, rG_Tory, 
         SG_Platy, vLG_Eulo, vLG_Mesopol, vLG_Mymarid, vLG_Platy, vLG_Tory, 
         vLG_egg, vLG_ecto, rG_ecto, # all gall-parasitoid interactions for both pairwise and dissimilarity analyses
         shootEst.no18,
         Genotype, 
         vLG_abund, rG_abund, # only gall species associated with significant pairwise gall-ptoid interactions
         vLG.height.mean) %>% # only gall species that varied in size among genotypes
  mutate(vLG_density = vLG_abund/shootEst.no18,
         rG_density = rG_abund/shootEst.no18)

df.vLG.size.sub <- na.omit(df) # will use for RDA when I include gall size in the analysis, because not all trees had vLG galls on them, hence the reduced dataset.

## RDA of insect food web
rda_full <- rda(df.vLG.size.sub[ ,interaxns_noPont]/df.vLG.size.sub$shootEst.no18 ~ vLG_density*vLG.height.mean + rG_density, data = df.vLG.size.sub) # much higher explanatory power without sqrt transforming the data.
RsquareAdj(rda_full)
summary(rda_full)
plot(rda_full, display = c("sp","bp"))
vif.cca(rda_full)
anova(rda_full, by = "margin") # why aren't the lower level main effects included?

## Linear and logistic regression models of pairwise gall-parasitoid interactions as a function of gall density and gall height. I used a conditional modelling approach to dissect apart the process influence gall-parasitoid densities. Specifically, I used a logistic regression to first model the presence/absence of a gall-parasitoid interaction. I then used linear regression on the presence-only data, to model the effect of gall density and gall height on the log-transformed density of gall-parasitoid interactions. This is a much more straightforward way of modelling data, creates nice residuals, and is easier to interpret than a mixture-model (e.g., zero-inflated poisson).

# vLG_egg. Conditional model
vLG_egg.pos.df <- filter(df, vLG_egg > 0)
vLG_egg.pos.lm <- lm(log(vLG_egg/shootEst.no18) ~ vLG_density*vLG.height.mean, vLG_egg.pos.df)
summary(vLG_egg.pos.lm)
plot(vLG_egg.pos.lm)

visreg(vLG_egg.pos.lm, xvar = "vLG_density", by = "vLG.height.mean", trans = exp)
visreg2d(vLG_egg.pos.lm, x = "vLG.height.mean", y = "vLG_density", trans = exp)

vLG_egg.pres.glm <- glm(vLG_egg > 0 ~ vLG_density, df, family = binomial)
summary(vLG_egg.pres.glm)
#plot(vLG_egg.pres.glm)
visreg(vLG_egg.pres.glm, scale = "response")

#t <- glmer(vLG_egg > 0 ~ 1 + (1), df, family = binomial)
t2 <- glmer(vLG_egg > 0 ~ 1 + (1|Genotype), df, family = binomial)
t4 <- glmmPQL(vLG_ecto > 0 ~ offset(log(shootEst.no18)), ~ 1|Genotype, df, family = binomial)
t3 <- glmmPQL(vLG_egg > 0 ~ offset(log(shootEst.no18)), ~ 1|Genotype, df, family = binomial)
plot(t3)
qqmath(ranef(t3)$'(Intercept)')
summary(t3)
1.53^2/(1.53^2 + 0.83*(pi^2)/3) # heritability of 0.462? pretty high...Make sure I'm estimating this correctly, but maybe this is the way to go?
summary(t2)
AIC(t2)
AIC(t)
library(pbkrtest)
PBmodcomp(t2, t)
profile(t2, which = 1)
library(arm)

# A function to extract simulated estimates of random effect paramaters from 
# lme4 objects using the sim function in arm
# whichel = the character for the name of the grouping effect to extract estimates for 
# nsims = the number of simulations to pass to sim
# x = model object
REsim <- function(x, whichel=NULL, nsims){
  require(plyr)
  mysim <- sim(x, n.sims = nsims)
  if(missing(whichel)){
    dat <- plyr::adply(mysim@ranef[[1]], c(2, 3), plyr::each(c(mean, median, sd)))
    warning("Only returning 1st random effect because whichel not specified")
  } else{
    dat <- plyr::adply(mysim@ranef[[whichel]], c(2, 3), plyr::each(c(mean, median, sd)))
  }
  return(dat)
}

t3 <- REsim(t2, whichel = "Genotype", nsims = 1000)

# Dat = results of REsim
# scale = factor to multiply sd by
# var = character of "mean" or "median"
# sd = character of "sd"
plotREsim <- function(dat, scale, var, sd){
  require(eeptools)
  dat[, sd] <- dat[, sd] * scale
  dat[, "ymax"] <- dat[, var] + dat[, sd] 
  dat[, "ymin"] <- dat[, var] - dat[, sd] 
  dat[order(dat[, var]), "id"] <- c(1:nrow(dat))
  ggplot(dat, aes_string(x = "id", y = var, ymax = "ymax", 
                         ymin = "ymin")) + 
    geom_pointrange() + theme_dpi() + 
    labs(x = "Group", y = "Effect Range", title = "Effect Ranges") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = 0, color = I("red"), size = I(1.1))
}
library(ggplot2)
library(eeptools)
plotREsim(t3, scale = 1.2, 
          var = "mean", sd = "sd")

# non-conditional linear model
summary(lm(sqrt(vLG_egg/shootEst.no18) ~ vLG_density*vLG.height.mean, df)) # qualitatively the same results as the conditional model

library(pscl)
summary(vLG_egg.zero <- zeroinfl(vLG_egg ~ offset(log(shootEst.no18)) +
                                   vLG_density*vLG.height.mean | 
                                   vLG_density, df))
plot(vLG_egg.zero)
AIC(vLG_egg.zero)
visreg(vLG_egg.zero, xvar = "vLG_density", by = "vLG.height.mean")
vLG_egg.lme <- lmer(log(vLG_egg/shootEst.no18) ~ vLG_density*vLG.height.mean  +
                      (1|Genotype), vLG_egg.pos.df)
summary(vLG_egg.lme)
exactRLRT(vLG_egg.lme)

vLG_egg.all.lme <- lmer(sqrt(vLG_egg/shootEst.no18) ~ vLG_density*vLG.height.mean +
                      (1|Genotype), df)
summary(vLG_egg.all.lme)
exactRLRT(vLG_egg.all.lme)

## vLG_ecto 
vLG_ecto.pos.df <- filter(df, vLG_ecto > 0)
vLG_ecto.pos.lm <- lm(log(vLG_ecto/shootEst.no18) ~ vLG_density + vLG.height.mean, vLG_ecto.pos.df)
summary(vLG_ecto.pos.lm)
#plot(vLG_ecto.pos.lm)

visreg(vLG_ecto.pos.lm, scale = "response")
visreg2d(vLG_ecto.pos.lm, x = "vLG.height.mean", y = "vLG_density")

vLG_ecto.pres.glm <- glm(vLG_ecto > 0 ~ vLG_density, df, family = binomial)
summary(vLG_ecto.pres.glm)
#plot(vLG_ecto.pres.glm)
visreg(vLG_ecto.pres.glm, scale = "response")

summary(lm(vLG_ecto/shootEst.no18 ~ vLG_density + vLG.height.mean, df)) # similar results, except vLG size is only marginally significant now

vLG_ecto.lme <- lmer(log(vLG_ecto/shootEst.no18) ~ vLG_density + vLG.height.mean +
                      (1|Genotype), vLG_ecto.pos.df)
summary(vLG_ecto.lme)
exactRLRT(vLG_ecto.lme)

vLG_ecto.all.lme <- lmer(sqrt(vLG_ecto/shootEst.no18) ~ vLG_density + vLG.height.mean +
                          (1|Genotype), df)
summary(vLG_ecto.all.lme)
exactRLRT(vLG_ecto.all.lme)

## rG_ecto 
rG_ecto.pos.df <- filter(df, rG_ecto > 0)
rG_ecto.pos.lm <- lm(log(rG_ecto/shootEst.no18) ~ rG_density, rG_ecto.pos.df)
summary(rG_ecto.pos.lm)
plot(rG_ecto.pos.lm)
visreg(rG_ecto.pos.lm, scale = "response")

rG_ecto.pres.glm <- glm(rG_ecto > 0 ~ rG_density, df, family = binomial)
summary(rG_ecto.pres.glm)
#plot(rG_ecto.pres.glm)
visreg(rG_ecto.pres.glm, scale = "response")

summary(lm(rG_ecto/shootEst.no18 ~ rG_density, df)) # qualitatively the same results as the conditional model.

rG_ecto.lme <- lmer(log(rG_ecto/shootEst.no18) ~ rG_density +
                       (1|Genotype), rG_ecto.pos.df)
summary(rG_ecto.lme)
exactRLRT(rG_ecto.lme)

rG_ecto.all.lme <- lmer((rG_ecto/shootEst.no18) ~ rG_density +
                           (1|Genotype), df)
summary(rG_ecto.all.lme)
exactRLRT(rG_ecto.all.lme)

##### I don't know if I need these genotype-level analyses. They are one step removed from the process that is actually going on in nature. Plus, if I include genotype as a random effect in the models, and it no longer explains a significant proportion of the variance, then I think I've captured the mechanism by which genotype is affecting the gall-parasitoid interaction.
##### Genotype-level analyses
source('~/Documents/Genotype_Networks/Rscripts/network_management_genotype_level.R')

df.geno <- genotype_level_interaxn_traits_size %>%
  select(aSG_Tory, rG_Eulo, rG_Lestodip, rG_Mesopol, rG_Platy, rG_Tory, 
         SG_Platy, vLG_Eulo, vLG_Mesopol, vLG_Mymarid, vLG_Platy, vLG_Tory, 
         vLG_egg, vLG_ecto, rG_ecto, # all gall-parasitoid interactions for both pairwise and dissimilarity analyses
         #shootEst.no18, # not needed because everything has been transformed into densities
         Genotype,
         vLG_abund, rG_abund, # Note that these are gall densities now! I already did this in the "network_management_genotype_level.R" file. Note that only gall species associated with significant pairwise gall-ptoid interactions
         vLG.height.mean = vLG.height.mean_mean.na.rm) # only gall species that varied in size among genotypes

df.geno.vLG.size.sub <- na.omit(df.geno) # will use for RDA when I include gall size in the analysis, because not all trees had vLG galls on them, hence the reduced dataset.

## RDA of insect food web. Genotype-level
rda.geno_full <- rda(df.geno.vLG.size.sub[ ,interaxns_noPont] ~ vLG_abund*vLG.height.mean + rG_abund, data = df.geno.vLG.size.sub) # much higher explanatory power without sqrt transforming the data.
RsquareAdj(rda.geno_full)
summary(rda.geno_full)
plot(rda.geno_full, display = c("sp","bp"))
vif.cca(rda.geno_full)
anova(rda.geno_full, by = "margin") # rG_abund is only marginally significant now. I'm keeping it in the model though.

## linear models

# vLG_egg
vLG_egg.geno.lm <- lm(vLG_egg ~ vLG_abund*vLG.height.mean, df.geno)
summary(vLG_egg.geno.lm)
plot(vLG_egg.geno.lm) # strong effect of #24 (genotype X)
visreg(vLG_egg.geno.lm, xvar = "vLG_abund", by = "vLG.height.mean")

# vLG_ecto
vLG_ecto.geno.lm <- lm(vLG_ecto ~ vLG_abund, df.geno)
summary(vLG_ecto.geno.lm)
plot(vLG_ecto.geno.lm) 
visreg(vLG_ecto.geno.lm, xvar = "vLG_abund")

# rG_ecto
rG_ecto.geno.lm <- lm(rG_ecto ~ rG_abund, df.geno)
summary(rG_ecto.geno.lm)
plot(rG_ecto.geno.lm) 
visreg(rG_ecto.geno.lm, xvar = "rG_abund")
