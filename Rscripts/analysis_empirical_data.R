## This code replicates the analyses of the empirical data presented in the Results and Discussion and Supplementary Material of the manuscript, "Genetic specificity of a plant-insect food web: Implications for linking genetic variation to food-web complexity"
## Code author: Matthew A. Barbour
## Email: barbour@zoology.ubc.ca

## source in required datasets and functions ----
#source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')
source('Rscripts/fxns_mvabund_diagnostics.R')
#source('~/Documents/miscellaneous_R/model_diagnostic_functions.R')
#source('~/Documents/miscellaneous_R/jensen_magnitude_function.R')

## load required libraries ----
library(vegan) # for adonis analysis
#library(knitr) # for making tables
library(ggplot2) # for exploratory plots
theme_set(theme_bw()) # customize ggplot for prettier default graphs
library(gridExtra) # for multipanel plots
library(mvabund) # multivariate analysis
library(HH) # for variance inflation factor analysis (vif) in multiple regression
library(visreg) # visualize regression output
library(dplyr) # for dataset

## load and prep datasets for analysis ----
full.df <- read.csv('data/tree_level_interaxn_all_plants_traits_size.csv') %>% 
  tbl_df()
interaxns_noPont <- full.df %>%
  select(aSG_Tory:rG_Platy, rG_Tory, SG_Platy, vLG_Eulo:vLG_Tory) %>%
  names() # all gall-parasitoid interactions used for analysis. Excluded interactions associated with the sawfly Pontania californica.
#full.df <- tree_level_interaxn_all_plants_traits_size 
#interaxns_noPont # all gall-parasitoid interactions used for analysis

## Gall size analyses ----
# Testing whether the size of each gall species varies among willow genotypes. 
# We were unable to use the mvabund framework for this analysis, because we did not always find all four gall species on the same sampled branch (i.e. missing data for gall size). Therefore, we analyzed separate linear models to test for these effects.
# There were too few data points to reliably test for differences in mid-Stem gall diameter (SG = Rabdophaga salicisbattatus), but we used weighted linear models to test for differences in gall size among willow genotypes. We weighted models by the number of galls collected from each plant replicate, because a higher number would reflect a more accurate measure of mean gall size for that particular willow.

## exploratory boxplots
vLG.height.mean.plot <- ggplot(data = full.df, aes(x = Genotype, y = vLG.height.mean)) + 
  geom_boxplot() + ylab("Iteomyia salicisverruca diameter (mm)")
rG.height.mean.plot <- ggplot(data = full.df, aes(x = Genotype, y = rG.height.mean)) + 
  geom_boxplot() + ylab("Rabdophaga salicisbrassicoides diameter (mm)")
aSG.height.mean.plot <- ggplot(data = full.df, aes(x = Genotype, y = aSG.height.mean)) + 
  geom_boxplot() + ylab("Cecidomyid sp. A diameter (mm)")
SG.height.mean.plot <- ggplot(data = full.df, aes(x = Genotype, y = SG.height.mean)) + 
  geom_boxplot() + ylab("Rabdophaga salicisbattatus diameter (mm)")

grid.arrange(vLG.height.mean.plot, rG.height.mean.plot, aSG.height.mean.plot, SG.height.mean.plot)

## models
vLG.size.lm <- lm(log(vLG.height.mean) ~ Genotype, data = full.df, weights = vLG.gall.count)
rG.size.lm <- lm(log(rG.height.mean) ~ Genotype, data = full.df, weights = rG.gall.count)
aSG.size.lm <- lm(log(aSG.height.mean) ~ Genotype, data = full.df, weights = aSG.gall.count)

## tests
anova(vLG.size.lm)
anova(rG.size.lm)
anova(aSG.size.lm)

## Gall abundance and composition analyses ----
# We used multivariate GLMs to examine whether gall abundances varied among willow genotypes. 

## exploratory boxplots of gall abundance among the different genotypes.
gall.plot.df <- as.data.frame(full.df) %>% # melt doesn't like 'tbl' class
  select(Genotype, 
         Iteomyia_salicisverruca = vLG_abund,
         Cecidomyiid_sp.A = aSG_abund,
         Rabdophaga_salicisbrassicoides = rG_abund, 
         Rabdophaga_salicisbattatus = SG_abund) %>%
  melt()

ggplot(gall.plot.df, aes(x = Genotype, y = value, color = variable)) + 
  geom_boxplot() + facet_wrap( ~ variable, nrow = 2) + theme_bw() + ylab("Gall density (no. per branch)")

## Assessing which error distribution is best to model the data. Negative binomial appear to provide the best fit.
gall.mvabund <- mvabund(full.df[ ,c("vLG_abund","rG_abund","aSG_abund","SG_abund")])
meanvar.plot(gall.mvabund ~ full.df$Genotype)
poisson_curve(from = 0.1, to = 10) # blue
quasipoisson_curve(from = 0.1, to = 10, quasi.scalar = 2) # red
neg.binomial_curve(from = 0.1, to = 10, theta.negbin = 0.7) # black

## specify the multivariate GLM
manyglm.gall <- manyglm(gall.mvabund ~ Genotype,
                        data = full.df,
                        family = "negative.binomial")

## examine residual plots of the model.
plot(manyglm.gall, which = 1:3) # residuals look really good. No heteroscadisticity or non-normality. Note that replotting the residuals gives qualitatively the same picture (it's important to replot them because the residuals involve random number generation, see ?plot.manyglm)

## Tests -- note that p-values may slightly differ from those reported in the manuscript since this uses a permutation.
anova.gall <- anova.manyglm(manyglm.gall, p.uni = "unadjusted")
anova.gall 

## We used analysis of dissimilarity (Bray-Curtis) to test for quantitative differences  in gall community composition.
gall.comm <- full.df %>% select(Genotype, vLG_abund, rG_abund, aSG_abund, SG_abund)
gall.comm.sub <- filter(gall.comm, rowSums(gall.comm[ ,-1]) > 0, 
                        Genotype != "J", Genotype != "N", Genotype != "U")
table(gall.comm.sub$Genotype) # restricted to Genotypes with greater than 2 replicates for composition analyses.

bray.gall.comm.sub <- vegdist(gall.comm.sub[ ,-1], method = "bray")

adonis(bray.gall.comm.sub ~ Genotype, data = gall.comm.sub)
anova(betadisper(bray.gall.comm.sub, group = gall.comm.sub$Genotype)) # no differences in dispersion among willow genotypes.
summary(meandist(bray.gall.comm.sub, grouping = gall.comm.sub$Genotype)) # average dissimilarity of 69%

## correlation in gall abundances
#gall.corrs.df <- full.df %>% select(Genotype, vLG.height.mean, vLG.gall.count, vLG_abund, rG_abund, aSG_abund)
#gall.corrs.pheno <- select(gall.corrs.df, -Genotype, -vLG.gall.count)
#library(psych)
#library(car)
#scatterplotMatrix(gall.corrs.pheno)
#corr.test(gall.corrs.pheno, use = "pairwise", method = "pearson")

#gall.corrs.geno <- gall.corrs.df %>%
 # group_by(Genotype) %>%
  #summarise(leaf_gall_size = weighted.mean(vLG.height.mean, vLG.gall.count, na.rm = TRUE),
   #         leaf_gall_abund = mean(vLG_abund, na.rm = TRUE),
    #        bud_gall_abund = mean(rG_abund, na.rm = TRUE),
     #       apical_stem_gall_abund = mean(aSG_abund, na.rm = TRUE)) %>%
  #select(leaf_gall_size, leaf_gall_abund, bud_gall_abund, apical_stem_gall_abund)
#scatterplotMatrix(gall.corrs.geno)
#corr.test(gall.corrs.geno)

## Plant trait - gall abundance analysis ----
# Determine which willow traits are resulting in variation in gall abundances.
# selects galls and traits for analysis
galls.traits.df <- full.df %>%
  select(aSG_abund:SG_abund,
         vLG.height.mean, vLG.gall.count,
         Trichome.No., Total_Area, Height, Density,
         C_N_imputed, water_content, specific_leaf_area, sal_tannin.PC1:flavanonOLES.PC1) 

# create dataset of galls and plant traits for analysis
gall.density.df <- galls.traits.df %>%
  mutate(log_size = log(Total_Area),
         log_trichomes = log(Trichome.No.+1)) %>%
  select(aSG_abund:SG_abund, Height:log_trichomes) %>%
  na.omit()

# Used sequential regression to reduce multicollinearity among non-principal component traits. For justification, see Graham 2003 (Ecology) and Barbour et al. 2015 (Functional Ecology).
gall.density.df$Height_resid = resid(lm(Height ~ log_size, gall.density.df))
gall.density.df$Density_resid = resid(lm(Density ~ log_size, gall.density.df))
gall.density.df$SLA_resid = resid(lm(specific_leaf_area ~ water_content, gall.density.df))

# subset the final trait-gall datasets
galls_galls.traits <- mvabund(select(gall.density.df, aSG_abund:SG_abund))
traits_galls.traits <- select(gall.density.df, C_N_imputed, water_content, sal_tannin.PC1:SLA_resid)

# specificy the multivariate GLM
mvabund.galls.traits <- manyglm(galls_galls.traits ~ C_N_imputed + water_content + sal_tannin.PC1 + cinn.PC1 + cinn.PC2 + flavonOLES.PC1 + flavonOLES.PC2 + flavanonOLES.PC1 + log_size + log_trichomes + Height_resid + Density_resid + SLA_resid, 
                                data = traits_galls.traits, 
                                family = "negative.binomial")

## plot residuals to assess model fit
plot(mvabund.galls.traits, which = 1:3)

## Starting with the most complex model, we progressively dropped predictors from the model that resulted in the lowest AIC for the model. We continued this in stepwise fashion until we arrived at the null model (no predictors).

# drop variables that reduce AIC the most and remove them in a stepwise fashion from analysis.
drop1(mvabund.galls.traits)
mod2 <- update(mvabund.galls.traits, .~. - sal_tannin.PC1); drop1(mod2)
mod3 <- update(mod2, .~. - flavonOLES.PC1); drop1(mod3)
mod4 <- update(mod3, .~. -water_content); drop1(mod4)
mod5 <- update(mod4, .~. - cinn.PC2); drop1(mod5)
mod6 <- update(mod5, .~. - cinn.PC1); drop1(mod6)
mod7 <- update(mod6, .~. -SLA_resid); drop1(mod7)
mod8 <- update(mod7, .~. -Height_resid); drop1(mod8)
mod9 <- update(mod8, .~. -Density_resid); drop1(mod9)
mod10 <- update(mod9, .~. -flavonOLES.PC2); drop1(mod10)
mod11 <- update(mod10, .~. -log_trichomes); drop1(mod11)
mod12 <- update(mod11, .~. -C_N_imputed); drop1(mod12)
mod13 <- update(mod12, .~. -log_size); drop1(mod13)
mod14.null <- update(mod13, .~. -flavanonOLES.PC1)

# use AIC to compare models
AIC.models.galls.traits <- data.frame(Model = 1:14, # progressive order in which models were fit
                                      Formula = c(paste(formula(mvabund.galls.traits)[3]), 
                                                  paste(formula(mod2)[3]),
                                                  paste(formula(mod3)[3]), 
                                                  paste(formula(mod4)[3]),
                                                  paste(formula(mod5)[3]), 
                                                  paste(formula(mod6)[3]),
                                                  paste(formula(mod7)[3]),
                                                  paste(formula(mod8)[3]), 
                                                  paste(formula(mod9)[3]),
                                                  paste(formula(mod10)[3]), 
                                                  paste(formula(mod11)[3]),
                                                  paste(formula(mod12)[3]),
                                                  paste(formula(mod13)[3]),
                                                  paste(formula(mod14.null)[3])),
                                      AIC = c(sum(AIC(mvabund.galls.traits)), 
                                              sum(AIC(mod2)), 
                                              sum(AIC(mod3)), 
                                              sum(AIC(mod4)),
                                              sum(AIC(mod5)), 
                                              sum(AIC(mod6)),
                                              sum(AIC(mod7)), 
                                              sum(AIC(mod8)), 
                                              sum(AIC(mod9)),
                                              sum(AIC(mod10)), 
                                              sum(AIC(mod11)),
                                              sum(AIC(mod12)), 
                                              sum(AIC(mod13)), 
                                              sum(AIC(mod14.null))))
arrange(AIC.models.galls.traits, AIC)[1:5,] # display all of the models that had lower AICs than the null model (including the null model).

# Likelihood ratio tests suggest that all of the models with a lower AIC than the null model provide a better fit to the data than the null.
anova.manyglm(mod14.null, mod13)
anova.manyglm(mod14.null, mod12)
anova.manyglm(mod14.null, mod11)
anova.manyglm(mod14.null, mod10)

# Likelihood ratio tests suggest that model 11 provides a better fit to the data than the two other closest models
anova.manyglm(mod13, mod11)
anova.manyglm(mod12, mod11)
anova.manyglm(mod10, mod11)
anova.manyglm(mod14.null, mod13, mod12, mod11, mod10)

# Identify gall-trait associations. C_N_imputed has a marginally significant effect on Iteomyia and R. salicisbrassicoides abundance. FlavononOLES.PC1 has a signficiant effect on Cecidomyiid sp. A abundance, and plant size as a significant effect on R. salicisbrassicoides and a marginal effect on R. salicisbattatus abundance.
uni.galls.traits <-anova.manyglm(mod11, p.uni = "unadjusted")

coef.df <- data.frame(coef(mod11))

range.flavanonOLES.PC1 <- with(gall.density.df, max(flavanonOLES.PC1) - min(flavanonOLES.PC1))
exp(coef.df["flavanonOLES.PC1","aSG_abund"]*range.flavanonOLES.PC1) # over the range of flavanonOLES.PC1, Cecidomyiid sp. A abundance increased 15-fold

range.C_N_imputed <- with(gall.density.df, max(C_N_imputed) - min(C_N_imputed))
exp(coef.df["C_N_imputed", "vLG_abund"]*range.C_N_imputed) # over the range of C:N ratio, Iteomyia abundance increased 3-fold (2.9)

exp(coef.df["C_N_imputed", "rG_abund"]*range.C_N_imputed) # over the range of C:N ratio, R. salicisbrassicoides abundance increased 8-fold (7.6)

1 - 1.10^coef.df["log_size","rG_abund"] # for every 10% increase in plant size, R. salicisbrassicoides density decreased 9%

1 - 1.10^coef.df["log_size","SG_abund"] # for every 10% increase in plant size, R. salicisbattatus density decreasted by 37%


## gall trait linear models
#gall.trait.lms <- manylm(galls_galls.traits ~ C_N_imputed + flavanonOLES.PC1 + log_size, 
#                         data = traits_galls.traits)
#summary(gall.trait.lms)
#anova(gall.trait.lms, p.uni = "unadjusted") # qualitatively the same results as GLM.

#best.r.sq(galls_galls.traits ~ traits_galls.traits$C_N_imputed + traits_galls.traits$water_content + traits_galls.traits$sal_tannin.PC1 + traits_galls.traits$cinn.PC1 + traits_galls.traits$cinn.PC2 + traits_galls.traits$flavonOLES.PC1 + traits_galls.traits$flavonOLES.PC2 + traits_galls.traits$flavanonOLES.PC1 + traits_galls.traits$log_size + traits_galls.traits$log_trichomes + traits_galls.traits$Height_resid + traits_galls.traits$Density_resid + traits_galls.traits$SLA_resid) # results are pretty close to GLM, retains flavanonOLES.PC1, C:N, but has cinn.PC1 and log.trichomes before log.size. Note that these variables were included in the GLM AIC, so it isn't too surprising. Not a big concern overall though.

## Visualize trait - abundance associations
# include model with all of the traits, but only examine one of the variables
vLG.trait <- MASS::glm.nb(vLG_abund ~ C_N_imputed + log_size + flavanonOLES.PC1, gall.density.df)
summary(vLG.trait) 
visreg(vLG.trait, xvar = "C_N_imputed", scale = "response")

aSG.trait <- MASS::glm.nb(aSG_abund ~ C_N_imputed + log_size + flavanonOLES.PC1, gall.density.df)
summary(aSG.trait)
visreg(aSG.trait, xvar = "flavanonOLES.PC1", scale = "response") # may or may not be too important

rG.trait <- MASS::glm.nb(rG_abund ~ C_N_imputed + log_size + flavanonOLES.PC1, gall.density.df)
summary(rG.trait)
visreg(rG.trait, xvar = "log_size", scale = "response")
visreg(rG.trait, xvar = "C_N_imputed", scale = "response")

SG.trait <- MASS::glm.nb(SG_abund ~ C_N_imputed + log_size + flavanonOLES.PC1, gall.density.df)
summary(SG.trait) # some warnings about lack of convergence
visreg(SG.trait, xvar = "log_size", scale = "response") # not a great fit. Outlier seems to be driving most of this relationship.

## Plant trait - gall size analyses ----
# only analyzing Iteomyia salicisverruca (vLG), because it was the only gall that varied significantly in size among will genotypes.
# generate dataset for analysis
gall.size.df <- galls.traits.df %>%
  mutate(log_size = log(Total_Area),
         log_trichomes = log(Trichome.No.+1)) %>%
  select(vLG.height.mean, vLG.gall.count, Height:log_trichomes) %>%
  na.omit()
gall.size.df$Height_resid = resid(lm(Height ~ log_size, gall.size.df))
gall.size.df$Density_resid = resid(lm(Density ~ log_size, gall.size.df))
gall.size.df$SLA_resid = resid(lm(specific_leaf_area ~ water_content, gall.size.df))
gall.size.df <- select(gall.size.df, vLG.height.mean, vLG.gall.count, 
                       C_N_imputed, water_content,
                       sal_tannin.PC1:SLA_resid)

# specify linear model
vLG.size.lm <- lm(vLG.height.mean ~ ., select(gall.size.df, -vLG.gall.count), weights = gall.size.df$vLG.gall.count)

# reiterate AIC stepwise procedure
drop1(vLG.size.lm)
size2 <- update(vLG.size.lm, .~. -C_N_imputed); drop1(size2)
size3 <- update(size2, .~. -Density_resid); drop1(size3)
size4 <- update(size3, .~. -flavonOLES.PC2); drop1(size4)
size5 <- update(size4, .~. -flavanonOLES.PC1); drop1(size5)
size6 <- update(size5, .~. -log_trichomes); drop1(size6)
size7 <- update(size6, .~. -SLA_resid); drop1(size7)
size8 <- update(size7, .~. -Height_resid); drop1(size8)
size9 <- update(size8, .~. -water_content); drop1(size9)
size10 <- update(size9, .~. -log_size); drop1(size10)
size11 <- update(size10, .~. -cinn.PC1); drop1(size11)
size12 <- update(size11, .~. -cinn.PC2); drop1(size12)
size13 <- update(size12, .~. -flavonOLES.PC1); drop1(size13)
size14.null <- update(size13, .~. -sal_tannin.PC1)

## use AIC to compare models
AIC.vLG.size.traits <- data.frame(Models = 1:14, # progressive order in which size models were fit
                                  Formula = c(paste(formula(vLG.size.lm)[3]), 
                                              paste(formula(size2)[3]),
                                              paste(formula(size3)[3]), 
                                              paste(formula(size4)[3]),
                                              paste(formula(size5)[3]), 
                                              paste(formula(size6)[3]),
                                              paste(formula(size7)[3]),
                                              paste(formula(size8)[3]), 
                                              paste(formula(size9)[3]),
                                              paste(formula(size10)[3]), 
                                              paste(formula(size11)[3]),
                                              paste(formula(size12)[3]),
                                              paste(formula(size13)[3]),
                                              paste(formula(size14.null)[3])),
                                  AIC = c(sum(AIC(vLG.size.lm)), 
                                          sum(AIC(size2)), 
                                          sum(AIC(size3)), 
                                          sum(AIC(size4)),
                                          sum(AIC(size5)), 
                                          sum(AIC(size6)),
                                          sum(AIC(size7)), 
                                          sum(AIC(size8)), 
                                          sum(AIC(size9)),
                                          sum(AIC(size10)), 
                                          sum(AIC(size11)),
                                          sum(AIC(size12)), 
                                          sum(AIC(size13)), 
                                          sum(AIC(size14.null))))
arrange(AIC.vLG.size.traits, AIC)

# all of these models give the same picture. Sal_tannin.PC1 and flavonOLES.PC1 are the best predictors.
summary(size11)
summary(size10)
summary(size12) # final model
summary(size9)

anova(size14.null, size12)

summary(lm(vLG.height.mean ~ luteolin.der2__A320nm, full.df, weights = full.df$vLG.gall.count))

anova(size14.null, size13, size12, size11, size10, size9, size8, size7)

plot(size12) # residuals looks pretty good
vLG.size.predict <- visreg(size12) # not currently helping calculate variation below.

# Iteomyia size decreased 1.2 fold over the range in salicylate/tannin chemistry
max(vLG.size.predict$sal_tannin.PC1$y$fit)/min(vLG.size.predict$sal_tannin.PC1$y$fit)

# Iteomyia size decreased 1.2 fold over the range in flavonid chemistry
max(vLG.size.predict$flavonOLES.PC1$y$fit)/min(vLG.size.predict$flavonOLES.PC1$y$fit)

## Gall-parasitoid frequency and composition analyses ----
#Here, we evaluate the assumptions of mvabund and see which error distribution is appropriate. Specifically, we first look at a plot of the mean-variance relationship of our response variables. It is easy to see that the negative binomial model (black line) provides the best fit to this data, suggesting that we should specify this as the error distribution in our model
full.mvabund <- mvabund(full.df[ ,interaxns_noPont]) # create mvabund object for easy analysis
full.meanvar <- meanvar.plot(full.mvabund ~ full.df$Genotype, table = TRUE) # easy to see that the data do not have a poisson distribution, therefore, a negative binomial may be a good fit.
poisson_curve(from = 0.1, to = 3.0, color = "blue")
quasipoisson_curve(from = 0.1, to = 3.0, quasi.scalar = 2, color = "red")
neg.binomial_curve(from = 0.1, to = 3.0, theta.negbin = 0.7, color = "black") # negative binomial appears to provide the best fit to this relationship

#no.relationships <- which(names(colSums(full.mvabund)) %in% c("rG_Lestodip", "rG_Mesopol", "SG_Platy", "vLG_Eulo"))
#sum(colSums(full.mvabund)[no.relationships])/sum(colSums(full.mvabund)) # less than 13% of the links in the network had relationships to gall abundance and size (see results below)

# We then fit a model and used residual plots to diagnose the model fit.
manyglm.full <- manyglm(full.mvabund ~ Genotype,
                        data = full.df,
                        family = "negative.binomial")
plot(manyglm.full, which = 1:3) # residuals aren't quite normally distributed, but there doesn't seem to be any heteroscedasticity in the model fit. Note that replotting the residuals gives qualitatively the same picture (it's important to replot them because the residuals involve random number generation, see ?plot.manyglm)

#Given that a negative binomial distribution seems to provide a good fit to the data, we tested whether the composition of gall-parasitoid interaction varied among willow genotypes. To further diagnose which interactions were driving this response, we conducted univariate analyses on each predictor, but adjusted for multiple comparisons. P-values were adjusted for multiple testing using a step-down resampling procedure. This methods provides strong control of family-wise error rates and makes use of resampling to ensure inferences take into account correlation between variables (Westfall & Young 1993).

#From the table, it is clear to say that genotype has a strong effect on the composition of links in the network. Moreover, differences in community composition are driven primarily by 3 interactions: vLG_Platy, vLG_Tory, and vLG_Mesopol.

anova.full <- anova.manyglm(manyglm.full, p.uni = "unadjusted") # Takes about 1 min and 30 sec to run.
anova.full # vLG_Platy, vLG_Tory, and vLG_Mesopol are driving the community response. rG_Tory is marginally significant and vLG_Eulo is close to marginal as well. aSG_Tory is significant as well.

## test with linear models
#link.lins <- manylm(log(full.mvabund+1) ~ Genotype, data = full.df)
#summary(link.lins)
#anova(link.lins, p.uni = "unadjusted") # qualitatively same results as GLM


## composition analysis
link.comm <- cbind.data.frame(Genotype = full.df$Genotype, full.df[ ,interaxns_noPont])

## filtered the data so that all willows have at least one trophic interaction (necessary for quantifying dissimilarity) and that the Genotypes included have 3 or more replicates
link.comm.sub <- filter(link.comm, rowSums(link.comm[ ,-1]) > 0, # must have at least one trophic interaction
                        Genotype %in% c("*","B", "D", "I", "K", "L", 
                                        "Q", "S", "V", "W", "X", "Y", "Z")) 
table(link.comm.sub$Genotype)

bray.link.comm.sub <- vegdist(link.comm.sub[ ,-1], method = "bray")

adonis(bray.link.comm.sub ~ Genotype, data = link.comm.sub)
anova(betadisper(bray.link.comm.sub, group = link.comm.sub$Genotype)) # no differences in dispersion among willow genotypes.
summary(meandist(bray.link.comm.sub, grouping = link.comm.sub$Genotype))

## Gall-parasitoid and gall density/size analyses ----
## Now we examine how variation in gall densities and gall size (for Iteomyia) affects the network.

# First we created a dataset that contained complete observations of the network and predictor variables
full.predictors <- full.df %>%
  select(aSG_Tory:rG_Platy, rG_Tory, SG_Platy, vLG_Eulo:vLG_Tory, # interactions
         aSG_abund:vLG_abund, # gall counts. Exclusing SG_abund because it didn't vary among genotypes
         Genotype, # keeping for later analysis
         vLG.height.mean) %>% # Iteomyia gall size
  # we transformed some of the variables so we could examine their vif 
  mutate(log.vLG_abund = log(vLG_abund),
         log.1.vLG_abund = log(vLG_abund + 1),
         log.1.rG_abund = log(rG_abund + 1),
         log.1.aSG_abund = log(aSG_abund + 1),
         log.vLG.height.mean = log(vLG.height.mean)) %>%
  na.omit() # removes all instances that don't have vLG.height.mean
dim(full.predictors)[1] # 81 data points

net.trait <- mvabund(full.predictors[ ,interaxns_noPont])


#We then looked for evidence of variance inflation among the predictor variables. But found little evidence for it.
#CURRENTLY NOT WORKING
vif(xx = as.data.frame(full.predictors[ ,c("log.1.aSG_abund","log.1.rG_abund","log.vLG_abund","vLG.height.mean")])) # Little evidence of variance inflation among predictor variables.


#We log transformed all predictor variables because it provided a much better fit to the data as determined by AIC. We then used AIC to compare our most complex model to our least complex. Instead of exploring all possible combinations, we started with the most complex model and used AIC to drop predictor variables with the lowest AIC values (i.e. least impact on removal). Using AIC, we identified 3 equivalenet models (difference in AIC < 2.1 amongst models). 
net.mvabund.1.full <- manyglm(net.trait ~ log.vLG_abund*vLG.height.mean + 
                                log.1.rG_abund + log.1.aSG_abund,
                              data = full.predictors,
                              family = "negative.binomial")

net.mvabund.2 <- update(net.mvabund.1.full, .~. - log.vLG_abund:vLG.height.mean)

drop1(net.mvabund.2)
net.mvabund.3 <- update(net.mvabund.2, .~. - vLG.height.mean)

drop1(net.mvabund.3)
net.mvabund.4 <- update(net.mvabund.3, .~. - log.1.aSG_abund)

drop1(net.mvabund.4)
net.mvabund.5 <- update(net.mvabund.4, .~. - log.vLG_abund)

drop1(net.mvabund.5)
net.mvabund.6.null <- update(net.mvabund.5, .~. -log.1.rG_abund)

AIC.models <- data.frame(Model = 1:6, # progressive order in which models were fit
                         Formula = c(paste(formula(net.mvabund.1.full)[3]), 
                                     paste(formula(net.mvabund.2)[3]),
                                     paste(formula(net.mvabund.3)[3]), 
                                     paste(formula(net.mvabund.4)[3]),
                                     paste(formula(net.mvabund.5)[3]), 
                                     paste(formula(net.mvabund.6.null)[3])),
                         AIC = c(sum(AIC(net.mvabund.1.full)), 
                                 sum(AIC(net.mvabund.2)), 
                                 sum(AIC(net.mvabund.3)), 
                                 sum(AIC(net.mvabund.4)),
                                 sum(AIC(net.mvabund.5)), 
                                 sum(AIC(net.mvabund.6.null))))
arrange(AIC.models, AIC)

anova.manyglm(net.mvabund.1.full, net.mvabund.2, net.mvabund.3, net.mvabund.4, net.mvabund.5, net.mvabund.6.null) # takes about 3 min to run.


# After deciding on the model with all main effects and no interactions, we examined the residuals and everything looked pretty good.
plot(net.mvabund.2, which = 1:3) # residuals look pretty good.


# We then examined the which processes underlied variation in the gall-parasitoid interactions.
anova.net <- anova.manyglm(net.mvabund.2, p.uni = "unadjusted") # takes about 2.5 min to run.
anova.net # not that since the P-values are determined by a resampling procedure, they may differ slightly between runs. Therefore, we retain all P-values < 0.10 (during at least one run) for the coefficient summary below.

t(coef(net.mvabund.2))[,-1][c(5,4,6,2,3,1,7),c(2,1,3,4)]

anova.net.null <- anova.manyglm(net.mvabund.2, net.mvabund.6.null) # takes about 30 sec to run
anova.net.null

sig.coef.df.2 <- mutate(melt(coef(net.mvabund.2)),
                        predictor_response = paste(X1, X2, sep = "_")) %>%
  # select(predictor = X1, response = X2, predictor_response, value) %>%
  subset(predictor_response %in% c("log.1.aSG_abund_aSG_Tory", 
                                   "log.1.rG_abund_rG_Eulo",
                                   "log.1.rG_abund_rG_Platy",
                                   "log.1.rG_abund_rG_Lestodip", # marginal
                                   "log.vLG_abund_SG_Platy", # marginal
                                   "log.vLG_abund_vLG_Mymarid",
                                   "vLG.height.mean_vLG_Eulo", # marginal
                                   "vLG.height.mean_rG_Tory", "log.1.rG_abund_rG_Tory",
                                   "vLG_abund_vLG_Mesopol", "vLG.height.mean_vLG_Mesopol", 
                                   "vLG.height.mean_vLG_Platy", "log.vLG_abund_vLG_Platy",
                                   "vLG.height.mean_vLG_Tory", "log.vLG_abund_vLG_Tory")) %>% # vLG.height.mean is marginal for vLG_Tory
  select(predictor = X1, response = X2, coefficient = value)
arrange(sig.coef.df.2, predictor, response)
1 - exp(-0.2157298) # 19% decrease in vLG_Platy interaction with every one unit increase in gall size.
1 - exp(-0.2697011) # 24% decrease in vLG_Mesopolobous interaction with every one unit increase in gall size.


# linear models
#net.trait.lms <- manylm(net.trait ~ log.vLG_abund + vLG.height.mean + 
 #                         log.1.rG_abund + log.1.aSG_abund,
  #                      data = full.predictors)
#summary(net.trait.lms) # non significant effect of gall height. It was the least important variable though in the GLM as well with mostly marginally significant effects.


# The results from the link composition analysis suggest that increases in gall abundance lead to increasing food web complexity for those individual nodes, whereas differences in leaf gall size lead to fundamental differences in link composition.

## Gall parasitism rates analysis ----
# Does the proportion of galls parasitized vary among willow genotypes?
vLG.ptized.glm <- glm(vLG_parasitized/vLG_abund ~ Genotype, data = full.df, 
                      weights = vLG_abund, family = "binomial")
vLG.Platy.glm <- glm(vLG_Platy/vLG_abund ~ Genotype, data = full.df, 
                     weights = vLG_abund, family = "binomial")
vLG.Mesopol.glm <- glm(vLG_Mesopol/vLG_abund ~ Genotype, data = full.df, 
                       weights = vLG_abund, family = "binomial")
vLG.Tory.glm <- glm(vLG_Tory/vLG_abund ~ Genotype, data = full.df, 
                    weights = vLG_abund, family = "binomial")
vLG.Eulo.glm <- glm(vLG_Eulo/vLG_abund ~ Genotype, data = full.df, 
                    weights = vLG_abund, family = "binomial")
vLG.Mymar.glm <- glm(vLG_Mymarid/vLG_abund ~ Genotype, data = full.df, 
                     weights = vLG_abund, family = "binomial")

rG.ptized.glm <- glm(rG_parasitized/rG_abund ~ Genotype, data = full.df, 
                     weights = rG_abund, family = "binomial")

aSG.ptized.glm <- glm(aSG_parasitized/aSG_abund ~ Genotype, data = full.df,
                      weights = aSG_abund, family = "binomial")


# Summary of results from parasitism models
anova(vLG.ptized.glm, test = "LR")
anova(vLG.Platy.glm, test = "LR")
anova(vLG.Mesopol.glm, test = "LR")
anova(vLG.Tory.glm, test = "LR")
anova(vLG.Eulo.glm, test = "LR")
anova(vLG.Mymar.glm, test = "LR")
anova(rG.ptized.glm, test = "LR")
anova(aSG.ptized.glm, test = "LR")

# Further explore the factors determining gall parasitism rates among willow genotypes. These were the best models as determined by AIC.
vLG_total.ptism <- glm(vLG_parasitized/vLG_abund ~ vLG.height.mean, 
                       data = full.df, weights = vLG_abund, family = "binomial")
summary(vLG_total.ptism)

vLG_Platy.ptism <- glm(vLG_Platy/vLG_abund ~ vLG.height.mean*vLG_abund, 
                       data = full.df, weights = vLG_abund, family = "binomial")
summary(vLG_Platy.ptism)

vLG_Mesopol.ptism <- glm(vLG_Mesopol/vLG_abund ~ vLG.height.mean*vLG_abund, 
                         data = full.df, weights = vLG_abund, family = "binomial")
summary(vLG_Mesopol.ptism)

vLG_Tory.ptism <- glm(vLG_Tory/vLG_abund ~ vLG.height.mean + vLG_abund, 
                      data = full.df, weights = vLG_abund, family = "binomial")
summary(vLG_Tory.ptism)

## summary of models. Used Anova function because they were multiple independent variables and we wanted to examine their marginal effects.
car::Anova(vLG_total.ptism)
car::Anova(vLG_Platy.ptism)
car::Anova(vLG_Mesopol.ptism)
car::Anova(vLG_Tory.ptism)

# The odds of a leaf gall being parasitized decreased by 25% with every one mm increase in gall size.
1-exp(coef(vLG_total.ptism)[2]) # 25% reduction in the odds of vLG being parasitized with every one unit increase in gall size.


