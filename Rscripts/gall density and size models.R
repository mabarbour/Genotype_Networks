# script borrowed heavily from older network_mechanism_models' data sheet

## source required script
source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')
source('~/Documents/Genotype_Willow_Community/datasets_&_Rscripts/functions_ms_willow_community.R') # used for testing for random effects of plant genotype on gall densities

## load required libraries
library(RLRsim) # for testing random effect. Using exactRLRT function gives more accurate P-values.
library(lme4) # for random effect models
library(ggplot2)

## function for evaluating normality of random effects assumption.
ggQQ_ranef <- function(ranef) # argument: vector of random effects
{
  y <- quantile(ranef, c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  p <- ggplot(data.frame(ranef = ranef), aes(sample=ranef)) +
    stat_qq(alpha = 0.5) +
    geom_abline(slope = slope, intercept = int, color="blue")
  
  return(p)
} 

## focused data frame to use for analyses in this script
abund.df <- tree_level_interaxn_all_plants_traits_size %>%
  select(Genotype, shootEst.no18, aSG_aSG.larv:gall_total_abund, vLG_egg, vLG_ecto, vLG.height.mean, rG.height.mean) 

geno_3plus <- c("*","B","D","I","K","L","Q","S","V","W","X","Y","Z") # Genotypes with 3 plus replicates after removing trees with zero gall-parasitoid interactions
#abund.df <- filter(abund.df, Genotype %in% geno_3plus)

## vLG_egg = Iteomyia egg parasitoid interactions (predominantly Platygaster)
abund.df.noX <- filter(abund.df, Genotype != "X")
vLG_egg.lmer <- lmer(sqrt(vLG_egg/shootEst.no18) ~ (1 | Genotype), abund.df)
summary(vLG_egg.lmer)
exactRLRT(vLG_egg.lmer)
genotype_random_model(x = abund.df, y = with(abund.df, sqrt(vLG_egg/shootEst.no18))) # H2 = 0.31. H2 drops to 0.21 when Genotype X is removed, and the trait distribution becomes much more normal.
vLG_egg.lmer.ranef <- ranef(vLG_egg.lmer)$Genotype$'(Intercept)'
ggQQ_ranef(vLG_egg.lmer.ranef) # look okay once Genotype X is removed.

## vLG_Ecto 
#vLG_ecto <- with(abund.df, vLG_Mesopol + vLG_Tory + vLG_Eulo)
vLG_ecto.lmer <- lmer(sqrt(vLG_ecto/shootEst.no18) ~ (1 | Genotype), abund.df)
summary(vLG_ecto.lmer)
exactRLRT(vLG_ecto.lmer)
genotype_random_model(x = abund.df, y = with(abund.df, sqrt(vLG_ecto/shootEst.no18))) # H2 = 0.2695621
vLG_ecto.lmer.ranef <- ranef(vLG_ecto.lmer)$Genotype$'(Intercept)'
ggQQ_ranef(vLG_ecto.lmer.ranef)

## vLG_Mesopol = Iteomyia-Mesopolobus interactions
vLG_Mesopol.lmer <- lmer(sqrt(vLG_Mesopol/shootEst.no18) ~ (1 | Genotype), abund.df)
summary(vLG_Mesopol.lmer)
exactRLRT(vLG_Mesopol.lmer)
genotype_random_model(x = abund.df, y = with(abund.df, sqrt(vLG_Mesopol/shootEst.no18))) # H2 = 0.119
vLG_Mesopol.lmer.ranef <- ranef(vLG_Mesopol.lmer)$Genotype$'(Intercept)'
ggQQ_ranef(vLG_Mesopol.lmer.ranef) # look okay.


## vLG_Tory = Iteomyia-Torymus interactions
vLG_Tory.lmer <- lmer(sqrt(vLG_Tory/shootEst.no18) ~ (1 | Genotype), abund.df)
summary(vLG_Tory.lmer)
exactRLRT(vLG_Tory.lmer)
genotype_random_model(x = abund.df, y = with(abund.df, sqrt(vLG_Tory/shootEst.no18))) # H2 = 0.291572
dotplot(ranef(vLG_Tory.lmer))
vLG_Tory.lmer.ranef <- ranef(vLG_Tory.lmer)$Genotype$'(Intercept)'
ggQQ_ranef(vLG_Tory.lmer.ranef) # doesn't look okay.

## rG_ecto
rG_ecto <- with(abund.df, rG_Mesopol + rG_Tory + rG_Eulo) 
#plot(rG_ecto ~ Genotype, abund.df)
rG_ecto.lmer <- lmer(sqrt(rG_ecto/shootEst.no18) ~ (1 | Genotype), abund.df)
summary(rG_ecto.lmer)
exactRLRT(rG_ecto.lmer)
genotype_random_model(x = abund.df, y = with(abund.df, sqrt(rG_ecto/shootEst.no18))) # H2 = 0.11
rG_ecto.lmer.ranef <- ranef(rG_ecto.lmer)$Genotype$'(Intercept)'
ggQQ_ranef(rG_ecto.lmer.ranef)

## rG_Tory = Rabdophaga salicisbrassicoides-Torymus interactions
abund.df.noKD <- filter(abund.df, Genotype != "K", Genotype != "D") # no estimate of variance due to Genotype once K and D are removed from the dataset.
rG_Tory.lmer <- lmer(sqrt(rG_Tory/shootEst.no18) ~ (1 | Genotype), abund.df)
summary(rG_Tory.lmer)
exactRLRT(rG_Tory.lmer)
genotype_random_model(x = abund.df, y = with(abund.df, sqrt(rG_Tory/shootEst.no18))) # H2 = 0.0928
dotplot(ranef(rG_Tory.lmer))
rG_Tory.lmer.ranef <- ranef(rG_Tory.lmer)$Genotype$'(Intercept)'
ggQQ_ranef(rG_Tory.lmer.ranef)

## Torymus abundance
Tory.lmer <- lmer(sqrt(Tory_abund/shootEst.no18) ~ (1 | Genotype), abund.df)
summary(Tory.lmer)
exactRLRT(Tory.lmer)
genotype_random_model(x = abund.df, y = with(abund.df, sqrt(Tory_abund/shootEst.no18))) 
dotplot(ranef(Tory.lmer))
Tory.lmer.ranef <- ranef(Tory.lmer)$Genotype$'(Intercept)'
ggQQ_ranef(Tory.lmer.ranef)

## for testing other interactions and parasitoid abundance
exactRLRT(lmer(sqrt(Eulo_abund/shootEst.no18) ~ (1 | Genotype), abund.df))

## vLG = Iteomyia salicisverruca.
vLG_lmer <- lmer(sqrt(vLG_abund/shootEst.no18*100) ~ (1 | Genotype),
                 data = abund.df)
summary(vLG_lmer)
exactRLRT(vLG_lmer)
genotype_random_model(x = abund.df, y = with(abund.df, sqrt(vLG_abund/shootEst.no18*100))) # H2 = 0.287

vLG.density.genotypes <- data.frame(Genotypes = rownames(coef(vLG_lmer)$Genotype),
                                 gall.sp = rep("Iteomyia salicisverruca", 26),
                                 gall.density = coef(vLG_lmer)$Genotype$'(Intercept)'^2)

## rG = Rabdophaga salicisbrassicoides
rG_lmer <- lmer(sqrt(rG_abund/shootEst.no18*100) ~ (1 | Genotype),
                 data = abund.df)
summary(rG_lmer)
exactRLRT(rG_lmer)
genotype_random_model(x = abund.df, y = with(abund.df, sqrt(rG_abund/shootEst.no18*100))) # H2 = 0.17

## Plot. NOTE THAT THESE BACK-TRANSFORMATIONS I'M USING MAY NOT BE ENTIRELY CORRECT...SQUARING THE MEAN ISN'T THE SAME AS SQUARING THE INDIVIDUAL DATA POINTS AND THEN TAKING THE MEAN.
rG.density.genotypes <- data.frame(Genotypes = rownames(coef(rG_lmer)$Genotype),
                                    gall.sp = rep("Rabdophaga salicisbrassicoides", 26),
                                    gall.density = coef(rG_lmer)$Genotype$'(Intercept)'^2) 

gall.density.genotypes <- rbind.data.frame(vLG.density.genotypes, rG.density.genotypes)

ggplot(gall.density.genotypes, aes(y = gall.density, x = gall.sp, fill = gall.sp)) + 
  geom_boxplot()

## vLG size. To test this, I took advantage of the full gall data set and used a nested random effect model.
vLG.size.df <- filter(gall.size.df, gall.sp == "vLG")
table(vLG.size.df$Genotype) # large heterogeneity in sample sizes
table(vLG.size.df$plant.position) # large heterogeneity in sample sizes

vLG.size.lmer <- lmer(gall.height ~ 1 + (1 | Genotype) + (1 |plant.position:Genotype),
                      data = vLG.size.df)
vLG.size.lmer.gen <- update(vLG.size.lmer, .~. - (1 | plant.position:Genotype)) # only Genotype as random effect
vLG.size.lmer.pp <- update(vLG.size.lmer, .~. - (1 | Genotype)) # only plant.position as random effect

vLG.size.lmer.sum <- summary(vLG.size.lmer)
summary(vLG.size.lmer.gen)
summary(vLG.size.lmer.pp)

gen.var <- vLG.size.lmer.sum$varcor$Genotype[1] # Variance due to Genotype
pp.var <- vLG.size.lmer.sum$varcor$'plant.position:Genotype'[1] # variance due to plant position
res.var <- 4.1133 # residual variance - don't know how to extract this from summary object.

H2.vLG.size <- gen.var/(gen.var + pp.var + res.var) # 0.147
exactRLRT(m = vLG.size.lmer.gen, mA = vLG.size.lmer, m0 = vLG.size.lmer.pp) # this model tests whether 

vLG.size.genotypes <- data.frame(Genotypes = rownames(coef(vLG.size.lmer)$Genotype),
                                 gall.sp = rep("Iteomyia salicisverruca", 24),
                                 gall.size = coef(vLG.size.lmer)$Genotype$'(Intercept)')

## rG size. To test this, I took advantage of the full gall data set and used a nested random effect model.
rG.size.df <- filter(gall.size.df, gall.sp == "rG")
table(rG.size.df$Genotype) # large heterogeneity in sample sizes
table(rG.size.df$plant.position) # large heterogeneity in sample sizes

rG.size.lmer <- lmer(gall.height ~ 1 + (1 | Genotype) + (1 |plant.position:Genotype),
                      data = rG.size.df)
rG.size.lmer.gen <- update(rG.size.lmer, .~. - (1 | plant.position:Genotype)) # only Genotype as random effect
rG.size.lmer.pp <- update(rG.size.lmer, .~. - (1 | Genotype)) # only plant.position as random effect

rG.size.lmer.sum <- summary(rG.size.lmer)
summary(rG.size.lmer.gen)
summary(rG.size.lmer.pp)

gen.var <- rG.size.lmer.sum$varcor$Genotype[1] # Variance due to Genotype
pp.var <- rG.size.lmer.sum$varcor$'plant.position:Genotype'[1] # variance due to plant position
res.var <- 1.80518 # residual variance - don't know how to extract this from summary object.

H2.rG.size <- gen.var/(gen.var + pp.var + res.var) # 0.0418
exactRLRT(m = rG.size.lmer.gen, mA = rG.size.lmer, m0 = rG.size.lmer.pp) # this model tests whether 

rG.size.genotypes <- data.frame(Genotypes = rownames(coef(rG.size.lmer)$Genotype),
                                 gall.sp = rep("Rabdophaga salicisbrassicoides", 22),
                                 gall.size = coef(rG.size.lmer)$Genotype$'(Intercept)')

## Plot of gall sizes
gall.size.genotypes <- rbind.data.frame(vLG.size.genotypes, rG.size.genotypes)
#gall.size.genotypes <- mutate(gall.size.genotypes, )

ggplot(gall.size.genotypes, aes(y = gall.size, x = gall.sp, fill = gall.sp)) + 
  geom_boxplot()# +
  #geom_text(aes(x = as.numeric(gall.sp) + 0.3), position = position_jitter(width = 0.1))


####### Random effect models of variation in gall densities among willow genotypes. Uses tree-level data. Consider using simulate, confint and boot in lme4 to get confidence intervals and everything for the parameters I'm interested in. Tried using rptR functions and used getME to extract X and Zt components, but it didn't seem to re-estimate the confidence intervals...

##### vLG = Iteomyia salicisverruca

### explore data

# histograms
ggplot(tree_level_interaxn_all_plants) + 
  geom_density(aes(x = vLG_abund), color = "red") + # count data
  geom_density(aes(x = rpois(n = length(tree_level_interaxn_all_plants$vLG_abund), lambda = mean(tree_level_interaxn_all_plants$vLG_abund))), color = "blue") + # simulated poisson data
  geom_density(aes(x = rnegbin(n = length(tree_level_interaxn_all_plants$vLG_abund), mu = mean(tree_level_interaxn_all_plants$vLG_abund), theta = 1.3))) # simulated negative binomial data. Used theta from glm model.

# box plots
vLG_base <- ggplot(tree_level_interaxn_all_plants, aes(x = Genotype)) 
vLG_base + geom_boxplot(aes(y = vLG_abund)) 
vLG_base + geom_boxplot(aes(y = shootEst.no18)) # definitely should account for variation in estimated number of shoots sampled.
vLG_base + geom_boxplot(aes(y = vLG_abund/shootEst.no18*100)) # at least a couple of large outliers

# identify outliers
with(tree_level_interaxn_all_plants, which(vLG_abund/shootEst.no18*100 > 10)) #  data points #17 and #121
tree_level_interaxn_all_plants[c(17,121), "plant.position"] # plant positions 59 and 373.

# genotype means
vLG_means <- tree_level_interaxn_all_plants %>%
  group_by(Genotype) %>%
  summarise(vLG_density_mean = mean(vLG_abund/shootEst.no18*100))

ggplot(vLG_means) +
  geom_density(aes(x = vLG_density_mean)) + # not normally distributed 
  geom_density(aes(x = sqrt(vLG_density_mean)), color = "red") # a bit better


# create alternative data frame. One where Genotypes with zero vLG are removed and the other with outliers removed.
vLG_noQU <- tree_level_interaxn_all_plants %>%
  filter(Genotype != "Q") %>% 
  filter(Genotype != "U")
vLG_noOutliers <- tree_level_interaxn_all_plants[-c(17,121),]
vLG_noOutliersNoQU <- tree_level_interaxn_all_plants[-c(17,121),]  %>%
  filter(Genotype != "Q") %>% 
  filter(Genotype != "U")

# final dataframe chosen for analysis
vLG.df <- vLG_noOutliers

vLG.df %>%
  group_by(Genotype) %>%
  summarise(mean.vLG.dens = mean(vLG_abund/shootEst.no18*100, na.rm = TRUE)) # 0 - 4.18 galls per 100 shoots.
4.18045463/0.05952599 # 70.2-fold

vLG_glm <- glm.nb(vLG_abund ~ offset(log(shootEst.no18)),
                  data = vLG.df)
summary(vLG_glm)

vLG_glmm <- glmmPQL(vLG_abund ~ offset(log(shootEst.no18)), random = ~ 1 | Genotype, data = vLG.df, family = negative.binomial(theta = 0.4482, link = log)) # Double check that theta value corresponds to data frame used in glm

ggplot(data = data.frame(resid = residuals(vLG_glmm, type = "pearson"), fit = fitted(vLG_glmm)), aes(x = fit, y = resid)) + 
  geom_point() +
  geom_hline(y = 0) +
  stat_smooth()
vLG_glmm_ranef <- ranef(vLG_glmm)$"(Intercept)"
ggQQ_ranef(vLG_glmm_ranef) # looks worse than lmer model below...

vLG_lmer <- lmer(sqrt(vLG_abund/shootEst.no18) ~ (1 | Genotype),
                 data = vLG.df)
summary(vLG_lmer)




ggplot(data = data.frame(resid = residuals(vLG_lmer, type = "pearson"), fit = fitted(vLG_lmer)), aes(x = fit, y = resid)) + 
  geom_point() +
  geom_hline(y = 0) +
  stat_smooth() # do these have to be homoscedastic if the I'm only testing a random effect?
vLG_ranef <- ranef(vLG_lmer)$Genotype$"(Intercept)"
ggQQ_ranef(vLG_ranef) # looks pretty good. Looks a bit worse on the low end of the tail, likely because there are a higher proportion of means with low gall abundance.

genotype_random_model(x = vLG.df, y = with(vLG.df, sqrt(vLG_abund/shootEst.no18))) # H2 = 0.36. 
exactRLRT(vLG_lmer) 


#### rG = Rabdophaga salicisbrassicoides

# box plots
rG_base <- ggplot(tree_level_interaxn_all_plants, aes(x = Genotype)) 
rG_base + geom_boxplot(aes(y = rG_abund)) 
rG_base + geom_boxplot(aes(y = rG_abund/shootEst.no18*100)) # one huge outlier for genotype S
ggplot(tree_level_interaxn_all_plants[-121, ], aes(x = Genotype))  + geom_boxplot(aes(y = rG_abund/shootEst.no18*100))

which(with(tree_level_interaxn_all_plants, rG_abund/shootEst.no18*100) > 40) # point #121. Same as one of the ones for vLG

# data frame for rG analyses
rG.noOutliers <- tree_level_interaxn_all_plants[-121, ]
rG.noOutliers.noZeroGenos <- rG.noOutliers %>%
  filter(Genotype != "T") %>%
  filter(Genotype != "U") %>%
  filter(Genotype != "V") %>%
  filter(Genotype != "E")
rG.df <- rG.noOutliers

rG.noOutliers %>%
  group_by(Genotype) %>%
  summarise(mean.rG.dens = mean(rG_abund/shootEst.no18*100, na.rm = TRUE)) 
3.29171346/0.05305478 # 62-fold

# random effect model
rG_lmer <- lmer(sqrt(rG_abund/shootEst.no18) ~ (1 | Genotype),
                data = rG.df)
summary(rG_lmer)

ggplot(data = data.frame(resid = residuals(rG_lmer, type = "pearson"), fit = fitted(rG_lmer)), aes(x = fit, y = resid)) + 
  geom_point() +
  geom_hline(y = 0) +
  stat_smooth() # do these have to be homoscedastic if the I'm only testing a random effect?
rG_ranef <- ranef(rG_lmer)$Genotype$"(Intercept)"
ggQQ_ranef(rG_ranef) # looks okay, except for on high end...

genotype_random_model(x = rG.df, y = with(rG.df, sqrt(rG_abund/shootEst.no18))) # H2 = 0.17
exactRLRT(rG_lmer) # P-values more accurate


