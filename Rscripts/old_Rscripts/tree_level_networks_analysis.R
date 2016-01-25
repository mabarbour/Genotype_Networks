######## Tree level Network, Community, and Species level analysis

#### upload necessary packages

source('~/Documents/Genotype_Willow_Community/datasets_&_Rscripts/functions_ms_willow_community.R') # used for testing for random effects of plant genotype on gall densities

# function for evaluating normality of random effects assumption.
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


# analysis
#library(mvabund)
#library(vegan)
library(visreg)
#library(bipartite)
#library(pscl)
library(car)
library(lme4)
library(MASS)

# data manipulation. order of libraries is important. always load dplyr last.
library(reshape2)
library(reshape)
library(plyr)
library(dplyr)

# plotting
library(ggplot2)

##### upload data

# interaction data
tree_level_interaxn_all_plants <- read.csv("~/Documents/Genotype_Networks/data/tree_level_interaxn_all_plants.csv")

table(tree_level_interaxn_all_plants$Genotype) # sample sizes
sum(table(tree_level_interaxn_all_plants$Genotype))

# plant trait data
tree_level_traits <- read.csv("~/Documents/Genotype_Networks/data/plant.trait.galls.2011.tree.df.csv")
tree_level_traits <- select(tree_level_traits, -X, -Genotype)

# join trait and interaction data
galls_traits <- left_join(tree_level_interaxn_all_plants, tree_level_traits, by = "plant.position")
galls_traits <- select(galls_traits, -Nitrogen, -Carbon, - C_N_ratio)
galls_traits <- filter(galls_traits, salicortin__A270nm >= 0 & Total_Area > 0 & Trichome.No. >= 0 & specific_leaf_area > 0) # remove NAs

# transform some of the plant traits for use in regression.
galls_traits <- transform(galls_traits, 
                                            log_size = log(Total_Area), 
                                            density_resid = residuals(lm(Density ~ log(Total_Area), galls_traits)),
height_resid = residuals(lm(Height ~ log(Total_Area), galls_traits)), 
log_trich = log(Trichome.No. + 1), 
sla_resid = residuals(lm(specific_leaf_area ~ water_content, galls_traits)))


####### Random effect models of variation in gall densities among willow genotypes. Consider using simulate, confint and boot in lme4 to get confidence intervals and everything for the parameters I'm interested in.

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

ggplot(data = data.frame(resid = residuals(vLG_lmer, type = "pearson"), fit = fitted(vLG_lmer)), aes(x = fit, y = resid)) + 
  geom_point() +
  geom_hline(y = 0) +
  stat_smooth() # do these have to be homoscedastic if the I'm only testing a random effect?
vLG_ranef <- ranef(vLG_lmer)$Genotype$"(Intercept)"
ggQQ_ranef(vLG_ranef) # looks pretty good. Looks a bit worse on the low end of the tail, likely because there are a higher proportion of means with low gall abundance.

genotype_random_model(x = vLG.df, y = with(vLG.df, sqrt(vLG_abund/shootEst.no18))) # H2 = 0.36. 


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


#### rsLG = Pontania californica

# box plots
rsLG_base <- ggplot(tree_level_interaxn_all_plants, aes(x = Genotype)) 
rsLG_base + geom_boxplot(aes(y = rsLG_abund)) 
rsLG_base + geom_boxplot(aes(y = rsLG_abund/shootEst.no18*100)) # no huge outliers

# data frame for rsLG analyses
rsLG.df <- tree_level_interaxn_all_plants

# random effect model
rsLG_lmer <- lmer(sqrt(rsLG_abund/shootEst.no18) ~ (1 | Genotype),
                data = rsLG.df)
summary(rsLG_lmer)

ggplot(data = data.frame(resid = residuals(rsLG_lmer, type = "pearson"), fit = fitted(rsLG_lmer)), aes(x = fit, y = resid)) + 
  geom_point() +
  geom_hline(y = 0) +
  stat_smooth() # do these have to be homoscedastic if the I'm only testing a random effect?
rsLG_ranef <- ranef(rsLG_lmer)$Genotype$"(Intercept)"
ggQQ_ranef(rsLG_ranef) # looks okay, except for on low end

genotype_random_model(x = rsLG.df, y = with(rsLG.df, sqrt(rsLG_abund/shootEst.no18))) # H2 = 0.20



#### aSG = Cecidomyiid sp. A (undescribed stem gall)

# box plots
aSG_base <- ggplot(tree_level_interaxn_all_plants, aes(x = Genotype)) 
aSG_base + geom_boxplot(aes(y = aSG_abund)) 
aSG_base + geom_boxplot(aes(y = aSG_abund/shootEst.no18*100)) # no huge outliers


# data frame for aSG analyses
aSG.df <- tree_level_interaxn_all_plants

# random effect model
aSG_lmer <- lmer(sqrt(aSG_abund/shootEst.no18) ~ (1 | Genotype),
                  data = aSG.df)
summary(aSG_lmer)

ggplot(data = data.frame(resid = residuals(aSG_lmer, type = "pearson"), fit = fitted(aSG_lmer)), aes(x = fit, y = resid)) + 
  geom_point() +
  geom_hline(y = 0) +
  stat_smooth() # do these have to be homoscedastic if the I'm only testing a random effect?
aSG_ranef <- ranef(aSG_lmer)$Genotype$"(Intercept)"
ggQQ_ranef(aSG_ranef) # deviation from normality on high end

genotype_random_model(x = aSG.df, y = with(aSG.df, sqrt(aSG_abund/shootEst.no18))) # H2 = 0.12

#### SG = Cecidomyiid sp. A (undescribed stem gall)

# box plots
SG_base <- ggplot(tree_level_interaxn_all_plants, aes(x = Genotype)) 
SG_base + geom_boxplot(aes(y = SG_abund)) 
SG_base + geom_boxplot(aes(y = SG_abund/shootEst.no18*100)) # note that genotype Y is associated with a module for SG_Platy, but that this was likely driven by a single outlying interaction.


# data frame for SG analyses
SG.df <- tree_level_interaxn_all_plants

# random effect model
SG_lmer <- lmer(sqrt(SG_abund/shootEst.no18) ~ (1 | Genotype),
                 data = SG.df)
summary(SG_lmer)

ggplot(data = data.frame(resid = residuals(SG_lmer, type = "pearson"), fit = fitted(SG_lmer)), aes(x = fit, y = resid)) + 
  geom_point() +
  geom_hline(y = 0) +
  stat_smooth() # do these have to be homoscedastic if the I'm only testing a random effect?
SG_ranef <- ranef(SG_lmer)$Genotype$"(Intercept)"
ggQQ_ranef(SG_ranef) # deviation from normality on high end

genotype_random_model(x = SG.df, y = with(SG.df, sqrt(SG_abund/shootEst.no18))) # H2 = 0.06


######### Mixed effect models of which plant traits determine variation in gall densities.

### Explore correlations
library(psych)
library(car)

gall_density_trait.df <- galls_traits %>%
  mutate(sqrt_vLG_density = sqrt(vLG_abund/shootEst.no18),
         sqrt_rG_density = sqrt(rG_abund/shootEst.no18),
         sqrt_rsLG_density = sqrt(rsLG_abund/shootEst.no18),
         sqrt_aSG_density = sqrt(aSG_abund/shootEst.no18)) %>%
  select(sqrt_vLG_density, sqrt_rG_density, sqrt_rsLG_density, sqrt_aSG_density, Total_Area:sla_resid) # didn't include SG, because it didn't display significant variation among genotypes.

gall_density_trait.df.noOutliers <- gall_density_trait.df[-c(15,93), ] # note that data points 15 & 93 were the original outliers for vLG_density that I removed from the previous analysis, and that point #93 was also an outlier for rG. So I decided to remove them both from this dataset

gall.trait.df <- gall_density_trait.df.noOutliers

# vLG correlations
scatterplotMatrix(x = select(gall.trait.df, sqrt_vLG_density, water_content, C_N_imputed:sla_resid))

vLG_density.lm <- lm(sqrt(vLG_abund/shootEst.no18) ~ C_N_imputed + I(C_N_imputed^2), galls_traits[-c(15,19),])
summary(vLG_density.lm)
plot(vLG_density.lm)
plot(residuals(vLG_density.lm) ~ Genotype, galls_traits[-c(15,19),])

vLG_density.glm <- MASS::glm.nb(vLG_abund ~ offset(log(shootEst.no18)) + C_N_imputed, galls_traits[-c(15,19),])
summary(vLG_density.glm)
plot(vLG_density.glm)
plot(residuals(vLG_density.glm) ~ Genotype, galls_traits[-c(15,19),])

pp <- factor(galls_traits$plant.position)
vLG_density.lmer <- lmer(sqrt(vLG_abund/shootEst.no18) ~  C_N_imputed + I(C_N_imputed^2) + (1 | Genotype), data = galls_traits[-c(15,93), ])
summary(vLG_density.lmer)
plot(vLG_density.lmer)

vLG_density.lmer.null <- lmer(sqrt(vLG_abund/shootEst.no18) ~  1 + (1 | Genotype), data = galls_traits[-c(15,93), ])

anova(vLG_density.lmer.null, vLG_density.lmer)

gall_trait.corrs <- corr.test(y = select(gall_density_trait.df, vLG_density:aSG_density), x = select(gall_density_trait.df, Total_Area:sla_resid))
gall_trait.corrs$r

test <- galls_traits


#vLG_density_nullnoQU <- MASS::glm.nb(vLG_abund ~ offset(log(shootEst.no18)), vLG_noQU)
#vLG_density_null <- MASS::glm.nb(vLG_abund ~ offset(log(shootEst.no18)), tree_level_interaxn_all_plants)
#vLG_density_null_noOutliers <- MASS::glm.nb(vLG_abund ~ offset(log(shootEst.no18)), vLG_noOutliers)

#vLG_density_glm_noQU <- MASS::glm.nb(vLG_abund ~ offset(log(shootEst.no18)) + Genotype, vLG_noQU)
vLG_density_glm <- glm.nb(vLG_abund ~ offset(log(shootEst.no18)) + Genotype, tree_level_interaxn_all_plants)
vLG_density_glm_noOutliers <- MASS::glm.nb(vLG_abund ~ offset(log(shootEst.no18)) + Genotype, vLG_noOutliers)
vLG_density_fewzeros <- MASS::glm.nb(vLG_abund ~ offset(log(shootEst.no18)) + Genotype, data = subset(vLG_noOutliers, Genotype != "Q" & Genotype != "U"  & Genotype !=  "A" & Genotype !=  "G" & Genotype !=  "N" & Genotype !=  "P")) # residuals still don't look great.

summary(vLG_density_glm_noQU)
summary(vLG_density_glm)
summary(vLG_density_glm_noOutliers)
summary(vLG_density_fewzeros)

anova(vLG_density_nullnoQU, vLG_density_glm_noQU)
anova(vLG_density_null, vLG_density_glm)
anova(vLG_density_null_noOutliers, vLG_density_glm_noOutliers)

plot(vLG_density_glm_noQU)
plot(vLG_density_glm)
plot(vLG_density_glm_noOutliers)
plot(vLG_density_fewzeros)

vLG_density_visreg_noQU <- visreg(vLG_density_glm_noQU, "Genotype", scale = "response", cond = list(shootEst.no18 = 200), ylim = c(0,40)) # 2 likely outliers
vLG_density_visreg <- visreg(vLG_density_glm, "Genotype", scale = "response", cond = list(shootEst.no18 = 200), ylim = c(0,40)) # 2 likely outliers
vLG_density_visreg_noOutliers <- visreg(vLG_density_glm_noOutliers, "Genotype", scale = "response", cond = list(shootEst.no18 = 200), ylim = c(0,25))

# evaluate the utility of zeroinflated model. I chose a zero-inflated vs. hurdle model, so I could model the "false zeros" in the data as a function of the number of shoots I sampled.
vLG_zeroinfl <- zeroinfl(vLG_abund ~ offset(log(shootEst.no18)) + Genotype | log(shootEst.no18), vLG_noOutliers, dist = "negbin") # this models zeros simply as a function of the number of shoots sampled. In other words, this is trying to detect the "false zeros" in the data. Negative binomial gave a substantially better fit that the poisson model. 
summary(vLG_zeroinfl)  # note that my modelling of the zeros is not significant, which may further justify me not having to use a zero inflated model
vLG_zeroinfl_visreg <- visreg(vLG_zeroinfl)
plot(residuals(vLG_zeroinfl, type = "pearson") ~ fitted(vLG_zeroinfl))

plot(residuals(vLG_density_glm_noQU, type = "pearson") ~ fitted(vLG_density_glm_noQU), col = vLG_noQU$Genotype)
text(residuals(vLG_density_glm_noQU, type = "pearson") ~ fitted(vLG_density_glm_noQU),label = vLG_noQU$Genotype)

AIC(vLG_zeroinfl)
AIC(vLG_density_glm_noOutliers)

vuong(vLG_density_glm_noOutliers, vLG_zeroinfl) # null models are indisinguishable, although AIC favors vLG_zeroinfl and likelihood statistic favors glm.nb


genotype_random_model(x = vLG_noOutliers, y = sqrt(vLG_noOutliers$vLG_abund/vLG_noOutliers$shootEst.no18))
plot(sqrt(vLG_noOutliers$vLG_abund/vLG_noOutliers$shootEst.no18) ~ Genotype, vLG_noOutliers)

vLG_fewzeros.df <- subset(vLG_noOutliers, Genotype != "Q" & Genotype != "U"  & Genotype !=  "A" & Genotype !=  "G" & Genotype !=  "N" & Genotype !=  "P")
vLG_density_lm <- lm(sqrt(vLG_abund/shootEst.no18) ~ Genotype, vLG_fewzeros.df) # wow, adj. R2 is remarkably close to Broad sense heritability estimate. Residuals look great aside from the heterogeneity in variance. Even after I remove genotypes with pretty low zero values, the model is still significant and the residuals look pretty good
plot(residuals(vLG_density_lm) ~ fitted(vLG_density_lm))
plot(residuals(vLG_density_lm) ~  vLG_fewzeros.df$Genotype)

# create a data frame with estimates of densities from models as well as from the original data at the tree level
vLG_raw_means <- tree_level_interaxn_all_plants %>%
  group_by(Genotype) %>%
  summarise(vLG_density_mean = mean(vLG_abund/shootEst.no18*200))
vLG_raw_means.noOutliers <- vLG_noOutliers %>%
  group_by(Genotype) %>%
  summarise(vLG_density_mean.noOutliers = mean(vLG_abund/shootEst.no18*200))

vLG_model_densities <- data.frame(Genotype = vLG_density_visreg$Genotype$x$xx,
                                   vLG_glm.nb.ALL = vLG_density_visreg$Genotype$y$fit,
                                   vLG_glm.nb.noOutliers = vLG_density_visreg_noOutliers$Genotype$y$fit,
                                  vLG_zeroinfl.noOutliers = vLG_zeroinfl_visreg$Genotype$y$fit) # not necessary to include noQU estimates because it doesn't alter the mean estimates of the other genotypes.
vLG_density_estimates <- join_all(list(vLG_model_densities, vLG_raw_means, vLG_raw_means.noOutliers), by = "Genotype")

plot(vLG_glm.nb.noOutliers ~ vLG_glm.nb.ALL, vLG_density_estimates) # visually shows how the mean estimates are altered by the outlying data points.
plot(vLG_glm.nb.noOutliers ~ vLG_density_mean.noOutliers, vLG_density_estimates)
plot(vLG_glm.nb.ALL ~ vLG_density_mean, vLG_density_estimates) # outliers appear to bias the glm model estimates even more (especially, Genotype S)

# rG
hist(tree_level_interaxn_all_plants$rG_abund)
hist(rpois(n = length(tree_level_interaxn_all_plants$rG_abund), lambda = mean(tree_level_interaxn_all_plants$rG_abund)))
hist(rnegbin(n = length(tree_level_interaxn_all_plants$rG_abund), mu = mean(tree_level_interaxn_all_plants$rG_abund), theta = 0.55)) # used theta from glm model

plot(rG_abund/shootEst.no18*200 ~ Genotype, tree_level_interaxn_all_plants[-121,]) # one major outlier
with(tree_level_interaxn_all_plants, which(rG_abund/shootEst.no18*200 > 20)) # one of the same outliers for vLG. Also, one of our lowest shoot estimates for the number of shoots sampled. Data point #121


rG_density_null <-  MASS::glm.nb(rG_abund ~ offset(log(shootEst.no18)), tree_level_interaxn_all_plants) 
rG_density_null_noOutlier <-  MASS::glm.nb(rG_abund ~ offset(log(shootEst.no18)), tree_level_interaxn_all_plants[-121, ]) 

rG_density_glm <- MASS::glm.nb(rG_abund ~ offset(log(shootEst.no18)) + Genotype, tree_level_interaxn_all_plants) 
rG_density_glm_noOutlier <- MASS::glm.nb(rG_abund ~ offset(log(shootEst.no18)) + Genotype, tree_level_interaxn_all_plants[-121, ]) 

summary(rG_density_glm)
summary(rG_density_glm_noOutlier)

anova(rG_density_null, rG_density_glm)
anova(rG_density_null_noOutlier, rG_density_glm_noOutlier)

plot(rG_density_glm)
plot(rG_density_glm_noOutlier)

rG_density_visreg <- visreg(rG_density_glm, "Genotype", scale = "response", cond = list(shootEst.no18 = 200), ylim = c(0,120)) # one major outlier
rG_density_visreg_noOutliers <- visreg(rG_density_glm_noOutlier, "Genotype", scale = "response", cond = list(shootEst.no18 = 200), ylim = c(0,25)) 

genotype_random_model(x = tree_level_interaxn_all_plants[-121, ], y = sqrt(tree_level_interaxn_all_plants$rG_abund[-121]/tree_level_interaxn_all_plants$shootEst.no18[-121]))

rG_density_lm <- lm(sqrt(rG_abund/shootEst.no18) ~ Genotype, tree_level_interaxn_all_plants[-121, ])
summary(rG_density_lm)
plot(rG_density_lm)

rG_removeZeroGenotypes = tree_level_interaxn_all_plants %>%
  filter(Genotype != "T") %>%
  filter(Genotype != "U") %>%
  filter(Genotype != "V") %>%
  filter(Genotype != "E") 

# evaluate the utility of zeroinflated model. I chose a zero-inflated vs. hurdle model, so I could model the "false zeros" in the data as a function of the number of shoots I sampled.
rG_zeroinfl <- zeroinfl(rG_abund ~ offset(log(shootEst.no18)) + Genotype | log(shootEst.no18), tree_level_interaxn_all_plants[-121, ], dist = "negbin") # this models zeros simply as a function of the number of shoots sampled. In other words, this is trying to detect the "false zeros" in the data. Negative binomial gave a substantially better fit that the poisson model. 
summary(rG_zeroinfl) # note that my modelling of the binomial part is not significant, which may further justify me not having to use a zero inflated model
rG_zeroinfl_visreg <- visreg(rG_zeroinfl)

AIC(rG_zeroinfl) # only marginally better than glm.nb
AIC(rG_density_glm_noOutlier)

vuong(rG_density_glm_noOutlier, rG_zeroinfl) # null models are indisinguishable, although both AIC and vuong show slight favor for rG_zeroinfl.

rGnegbin <- MASS::glm.nb(rG_abund ~ offset(log(shootEst.no18)) + Genotype, rG_removeZeroGenotypes)
summary(rGnegbin)
anova(MASS::glm.nb(rG_abund ~ offset(log(shootEst.no18)), rG_removeZeroGenotypes), rGnegbin) # still significant differences among genotypes in mean abundances even after removing the genotypes with zeros.

rG_raw_means <- tree_level_interaxn_all_plants %>%
  group_by(Genotype) %>%
  summarise(rG_density_mean = mean(rG_abund/shootEst.no18*200))

rG_raw_means.noOutliers <- tree_level_interaxn_all_plants[-121, ] %>%
  group_by(Genotype) %>%
  summarise(rG_density_mean.noOutliers = mean(rG_abund/shootEst.no18*200))

rG_model_densities <- data.frame(Genotype = rG_density_visreg$Genotype$x$xx,
                                  rG_glm.nb.ALL = rG_density_visreg$Genotype$y$fit,
                                  rG_glm.nb.noOutliers = rG_density_visreg_noOutliers$Genotype$y$fit,
                                 rG_zeroinfl.noOutliers = rG_zeroinfl_visreg$Genotype$y$fit) 
rG_density_estimates <- join_all(list(rG_model_densities, rG_raw_means, rG_raw_means.noOutliers), by = "Genotype")

scatterplotMatrix(rG_density_estimates[ ,-1]) # outlying data point really appears to bias the Genotype's mean estimate.


# rsLG
hist(tree_level_interaxn_all_plants$rsLG_abund) # possible zero inflation
plot(rsLG_abund/shootEst.no18*200 ~ Genotype, tree_level_interaxn_all_plants) # no really clear outliers.

rsLG_density_null <-  MASS::glm.nb(rsLG_abund ~ offset(log(shootEst.no18)), tree_level_interaxn_all_plants)
rsLG_density_glm <- MASS::glm.nb(rsLG_abund ~ offset(log(shootEst.no18)) + Genotype, tree_level_interaxn_all_plants)
summary(rsLG_density_glm)
anova(rsLG_density_null, rsLG_density_glm) 
plot(rsLG_density_glm)

rsLG_density_visreg <- visreg(rsLG_density_glm, "Genotype", scale = "response", cond = list(shootEst.no18 = 200), ylim = c(0,10))

genotype_random_model(x = tree_level_interaxn_all_plants, y = sqrt(tree_level_interaxn_all_plants$rsLG_abund/tree_level_interaxn_all_plants$shootEst.no18))

rsLG_removeZeroGenotypes = tree_level_interaxn_all_plants %>%
  filter(Genotype != "R") %>%
  filter(Genotype != "P") %>%
  filter(Genotype != "N") %>%
  filter(Genotype != "E") %>%
  filter(Genotype != "T") %>%
  filter(Genotype != "Y") %>%
  filter(Genotype != "Z")

# evaluate the utility of zeroinflated model. I chose a zero-inflated vs. hurdle model, so I could model the "false zeros" in the data as a function of the number of shoots I sampled.
rsLG_zeroinfl <- zeroinfl(rsLG_abund ~ offset(log(shootEst.no18)) + Genotype | log(shootEst.no18), tree_level_interaxn_all_plants, dist = "negbin") # this models zeros simply as a function of the number of shoots sampled. In other words, this is trying to detect the "false zeros" in the data. Negative binomial gave a substantially better fit that the poisson model. 
summary(rsLG_zeroinfl) # note that my modelling of the binomial part is not significant, which may further justify me not having to use a zero inflated model
rsLG_zeroinfl_visreg <- visreg(rsLG_zeroinfl)

AIC(rsLG_zeroinfl) 
AIC(rsLG_density_glm) # favored over zeroinfl

vuong(rsLG_density_glm, rsLG_zeroinfl) # null models are indisinguishable

rsLGnegbin <- MASS::glm.nb(rsLG_abund ~ offset(log(shootEst.no18)) + Genotype, rsLG_removeZeroGenotypes)
rsLGhurdle <- hurdle(rsLG_abund ~ offset(log(shootEst.no18)) + Genotype, rsLG_removeZeroGenotypes, dist = "negbin") # need to remove genotypes with only zeros, because otherwise the solution is singular
rsLGzeroinfl <- zeroinfl(rsLG_abund ~ offset(log(shootEst.no18)) + Genotype, rsLG_removeZeroGenotypes, dist = "negbin") # error - system is exactly singular

AIC(rsLG_density_glm)
AIC(rsLGnegbin) # outperforms hurdle and zeroinflated model.
AIC(rsLGhurdle)
AIC(rsLGzeroinfl) 

#summary(rsLGnegbin)
#anova(MASS::glm.nb(rsLG_abund ~ offset(log(shootEst.no18)), rsLG_removeZeroGenotypes), rsLGnegbin) # no longer significant if zero genotypes are removed.

rsLG_raw_means <- tree_level_interaxn_all_plants %>%
  group_by(Genotype) %>%
  summarise(rsLG_density_mean = mean(rsLG_abund/shootEst.no18*200))

rsLG_model_densities <- data.frame(Genotype = rsLG_density_visreg$Genotype$x$xx,
                                 rsLG_glm.nb.ALL = rsLG_density_visreg$Genotype$y$fit) 
rsLG_density_estimates <- left_join(rsLG_model_densities, rsLG_raw_means, by = "Genotype")


# aSG
table(tree_level_interaxn_all_plants$aSG_abund)
hist(tree_level_interaxn_all_plants$aSG_abund) # binomial model may be more appropriate given the infrequent number of counts greater than 1
hist(rpois(n = length(tree_level_interaxn_all_plants$aSG_abund), lambda = mean(tree_level_interaxn_all_plants$aSG_abund)))
hist(rnegbin(n = length(tree_level_interaxn_all_plants$aSG_abund), mu = mean(tree_level_interaxn_all_plants$rG_abund), theta = 0.55)) # used theta from glm model

plot(aSG_abund/shootEst.no18*200 ~ Genotype, tree_level_interaxn_all_plants)
with(tree_level_interaxn_all_plants, which(aSG_abund/shootEst.no18*200 > 7))

aSG.geno.nonzeros <- subset(tree_level_interaxn_all_plants, 
                            Genotype != "A" &
                            Genotype != "M" &
                            Genotype != "P" &
                            Genotype != "T" &
                            Genotype != "U" &
                            Genotype != "W" &
                            Genotype != "Z")

aSG_density_glm <- glm(aSG_abund ~ offset(log(shootEst.no18)) + Genotype, tree_level_interaxn_all_plants, family = "poisson") 
summary(aSG_density_glm) # underdispersion in poisson model, suggesting it may not be appropriate
plot(aSG_density_glm)
anova(aSG_density_glm, test = "Chi") # still significant after removing "zero value" genotypes.

aSG_bin_glm <- glm(aSG_abund > 0 ~ offset(log(shootEst.no18)) + Genotype, aSG.geno.nonzeros, family = "binomial") 
summary(aSG_bin_glm) # no overdispersion in poisson model
plot(aSG_bin_glm)
anova(aSG_bin_glm, test = "Chi") # still significant with binomial model. Binomial model no longer significant after removing "zero value" genotypes. 

aSG_density_visreg <- visreg(aSG_density_glm, "Genotype", scale = "response", cond = list(shootEst.no18 = 200), ylim = c(0,15))

genotype_random_model(x = tree_level_interaxn_all_plants, y = sqrt(tree_level_interaxn_all_plants$aSG_abund/tree_level_interaxn_all_plants$shootEst.no18))

aSG_raw_means <- tree_level_interaxn_all_plants %>%
  group_by(Genotype) %>%
  summarise(aSG_density_mean = mean(aSG_abund/shootEst.no18*200))

summary(lm(sqrt(aSG_abund/shootEst.no18) ~ Genotype, aSG.geno.nonzeros)) # still a bit of heteroscedasticity after removing genotypes with zero values. 

aSG_model_densities <- data.frame(Genotype = aSG_density_visreg$Genotype$x$xx,
                                   aSG_glm.nb.ALL = aSG_density_visreg$Genotype$y$fit) 
aSG_density_estimates <- left_join(aSG_model_densities, aSG_raw_means, by = "Genotype")


# SG
plot(SG_abund/shootEst.no18*200 ~ Genotype, tree_level_interaxn_all_plants)

SG_density_glm <- glm(SG_abund > 0 ~ offset(log(shootEst.no18)) + Genotype, tree_level_interaxn_all_plants, family = "binomial") 
summary(SG_density_glm) # little overdispersion in poisson model. Same results with quasipoisson. iteration for theta doesn't converge for glm.nb so its likely not necessary.
plot(SG_density_glm)
anova(SG_density_glm, test = "Chi") # I don't trust either the quasipoisson or poisson models... Not significant if you use a binomial model

SG_density_visreg <- visreg(SG_density_glm, "Genotype", scale = "response", cond = list(shootEst.no18 = 200), ylim = c(0,100)) # I think some of these density estimates are highly unrealistic...

genotype_random_model(x = tree_level_interaxn_all_plants, y = sqrt(tree_level_interaxn_all_plants$SG_abund/tree_level_interaxn_all_plants$shootEst.no18))

# didn't obtain density estimates for SG because they were too low and not significant

##### compile data set for gall density estimates at the genotype level through a variety of methods

gall_densities <- join_all(list(vLG_density_estimates, rG_density_estimates, rsLG_density_estimates, aSG_density_estimates), by = "Genotype")
write.csv(gall_densities, "~/Documents/Genotype_Networks/data/gall_density_estimates_genotype_level.csv")

### Figures for vLG and rG abund
theme <- theme_bw() + theme(text = element_text(family = "Verdana", size = 24), legend.key=element_blank(), axis.title.x = element_text(vjust = -0.5), axis.title.y = element_text(vjust = 0.5), axis.line=element_line(colour="black"), panel.grid=element_blank(), panel.border=element_blank(), legend.title=element_blank())

# Total herbivore abundance
gall.geno.sum <- tree_level_interaxn_all_plants %>% # removing plant position #121 for rG, because of the ridiculously high value
  group_by(Genotype) %>%
  summarise(vLG.dens.mean = mean(vLG_abund/shootEst.no18*200), vLG.dens.se = sd(vLG_abund/shootEst.no18*200)/sqrt(length(vLG_abund/shootEst.no18*200)))

gall.geno.sum.rG <- tree_level_interaxn_all_plants[-121, ] %>% # removing plant position #121 for rG, because of the ridiculously high value
  group_by(Genotype) %>%
  summarise(rG.dens.mean = mean(rG_abund/shootEst.no18*200), rG.dens.se = sd(rG_abund/shootEst.no18*200)/sqrt(length(rG_abund/shootEst.no18*200)), N = n()) %>%
  arrange(rG.dens.mean)

gall.geno.sum.rG <- mutate(gall.geno.sum.rG, Geno.ord = factor(Genotype, levels = as.character(Genotype), ordered = TRUE))
levels(gall.geno.sum.rG$Geno.ord)[24] <- "C"

gall.geno.sum <- mutate(gall.geno.sum, Geno.ord = factor(Genotype, levels = as.character(Genotype), ordered = TRUE))
levels(gall.geno.sum$Geno.ord)[10] <- "C" # for aesthetics

ggplot(gall.geno.sum, aes(x = Geno.ord, y = vLG.dens.mean)) + 
  geom_errorbar(aes(ymin = vLG.dens.mean - vLG.dens.se, ymax = vLG.dens.mean + vLG.dens.se), width = 0.2) +
  geom_point(shape = 21, fill = "steelblue", size = 10) + 
  theme + 
  xlab ("Genotype") + 
  ylab ("No. individuals per 200 shoots")
ggsave("~/Documents/Genotype_Networks/vLG_density_variation.png", width = 11.5, height = 8, units = "in", dpi = 300)

rG_plot <- ggplot(gall.geno.sum.rG, aes(x = Geno.ord, y = rG.dens.mean)) + 
  geom_errorbar(aes(ymin = rG.dens.mean - rG.dens.se, ymax = rG.dens.mean + rG.dens.se), width = 0.2) +
  geom_point(shape = 21, fill = "steelblue", size = 10) + 
  theme + 
  xlab ("Genotype") + 
  scale_y_continuous(breaks = seq(0,10,2), limits = c(0,10)) + 
  ylab ("No. individuals per 200 shoots")
ggsave("~/Documents/Genotype_Networks/rG_density_variation.png", width = 11.5, height = 8, units = "in", dpi = 300)

library(gridExtra)
gp1<- ggplot_gtable(ggplot_build(vLG_plot))
gp2<- ggplot_gtable(ggplot_build(rG_plot))

maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
gp1$widths[2:3] <- maxWidth
gp2$widths[2:3] <- maxWidth

grid.arrange(gp1, gp2, ncol = 1)

#ggsave(vLG_rG_density_variation, "~/Documents/Genotype_Networks/vLG_rG_density_variation.png", width = 11.5, height = 8, units = "in", dpi = 300)



### Explore tree-level data
pdf("~/Documents/Genotype_Networks/boxplot_gall_parasitoid_responses_to_genotype.pdf")
for(i in 7:length(tree_level_interaxn_all_plants)){
  plot(tree_level_interaxn_all_plants[ ,i] ~ tree_level_interaxn_all_plants$Genotype, xlab = "Genotype", ylab = paste(names(tree_level_interaxn_all_plants)[i]))
}
dev.off()


### Does gall community composition vary among willow genotypes?
gall.comm <- select(tree_level_interaxn_all_plants, Genotype, aSG_abund:rsLG_abund)
gall.comm <- filter(gall.comm, rowSums(gall.comm[ ,-1]) > 0)
table(gall.comm$Genotype) # J, N, U have sample sizes less than 3

scatterplotMatrix(log(gall.comm[,-1]+1))

gall.comm.genotype <- gall.comm %>%
  group_by(Genotype) %>%
  summarise_each(funs(sum)) %>%
  filter(Genotype %in% c("A","V","X")) %>%
  mutate(total_abund = aSG_abund + rG_abund + vLG_abund + SG_abund + rsLG_abund) 
gall.com.genotype.prop <- round(gall.comm.genotype[,2:6]/gall.comm.genotype$total_abund,2)


gall.comm.sub <- filter(gall.comm, Genotype %in% c("*","A","B","D","E","F","G","H","I","K","L","M","O","P","Q","R","S","T","V","W","X","Y","Z"))

gall.comm.horn <- vegdist(gall.comm[ ,-1], method = "horn")
gall.comm.horn.sub <- vegdist(gall.comm.sub[ ,-1], method = "horn") # results robust to removal of genotypes with low sample sizes

adonis(gall.comm[ ,-1] ~ Genotype, gall.comm, method = "horn")
adonis(gall.comm.sub[ ,-1] ~ Genotype, gall.comm.sub, method = "horn")

gall.comm.cap <- capscale(gall.comm[ ,-1] ~ Genotype, gall.comm, distance = "horn")
plot(gall.comm.cap, display = c("cn","sp"))

gall.comm.meandist <- meandist(gall.comm.horn, grouping = gall.comm$Genotype)
low <- lower.tri(gall.comm.meandist)
gall.comm.meandist.lower <- gall.comm.meandist[which(low == TRUE)]
mean(gall.comm.meandist.lower) # 0.60 mean dissimilarity among genotypes
mean(diag(gall.comm.meandist)[-21]) # 0.53 mean dissimilarity within genotypes.
hist(diag(gall.comm.meandist)[-21])
hist(gall.comm.meandist.lower)
quantile(gall.comm.meandist.lower)

geno.N <- as.data.frame(table(gall.comm$Genotype))
plot(diag(gall.comm.meandist) ~ geno.N$Freq) # shows that estimation of dissimilarity within genotypes decreases with increasing sample size.

### Dissimilarity in parasitoid community composition
ptoid.comm <- select(tree_level_interaxn_all_plants, Genotype, Platy_abund:Mymarid_abund)
ptoid.comm.no.zero <- filter(ptoid.comm, rowSums(ptoid.comm[ ,-1]) > 0)
table(ptoid.comm.no.zero$Genotype)

ptoid.comm.sub <- filter(ptoid.comm.no.zero, Genotype %in% c("*","B","D","I","K","L","Q","S","V","W","X","Y","Z")) # retained genotypes with 3+ replicates. Note that this may restrict the analysis to genotypes that are more similar to each other...

adonis(ptoid.comm[ ,-1] ~ Genotype, ptoid.comm, method = "euclidean")
adonis(ptoid.comm.no.zero[ ,-1] ~ Genotype, ptoid.comm.no.zero, method = "horn")
adonis(ptoid.comm.sub[ ,-1] ~ Genotype, ptoid.comm.sub, method = "horn") # marginally significant if I remove genotypes with less than 3 replicates

plot(rda(ptoid.comm[,-1] ~ Genotype, ptoid.comm), display = c("sp","cn"))
ptoid.betadisper <- betadisper(vegdist(ptoid.comm[,-1], method = "euclidean"), ptoid.comm$Genotype, bias.adjust = TRUE)
anova(ptoid.betadisper)

scatterplotMatrix(log(ptoid.comm[,-1]+1))

### Which traits best predict gall community composition
gall.comm.trait <- select(galls_traits, aSG_abund:rsLG_abund, log_size:sla_resid, C_N_imputed:flavanonOLES.PC1, water_content)
galls <- select(gall.comm.trait, aSG_abund:rsLG_abund)
gall.comm.trait <- filter(gall.comm.trait, rowSums(galls) > 0)
traits <- select(gall.comm.trait, log_size:sla_resid, C_N_imputed:flavanonOLES.PC1, water_content)

library(car)
vif(gall.trait.cap)
gall.trait.null <- capscale(gall.comm.trait[,1:5] ~ 1, traits, distance = "horn")
gall.trait.cap <- capscale(gall.comm.trait[,1:5] ~ ., traits, distance = "horn")
summary(gall.trait.cap)
anova(gall.trait.cap)
plot(gall.trait.cap)

ordiR2step(gall.trait.null, scope = formula(gall.trait.cap))

### rG_Tory
rG_dom <- tree_level_interaxn_all_plants$rG_abund/tree_level_interaxn_all_plants$gall_total_abund
plot(tree_level_interaxn_all_plants$rG_Tory ~ rG_dom)
rG_dens <- tree_level_interaxn_all_plants$rG_abund/tree_level_interaxn_all_plants$shootEst.no18

rG_bin <- glm(cbind(rG_abund, gall_total_abund - rG_abund) ~ Genotype, tree_level_interaxn_all_plants[-c(53,52),], family = "binomial")
summary(rG_bin)
plot(rG_bin)

library(ggplot2)
ggplot(tree_level_interaxn_all_plants, aes(x = rG_dens, y = rG_Tory)) + stat_smooth(method = "glm", family = "poisson")
vLG_dom <- tree_level_interaxn_all_plants$vLG_abund/tree_level_interaxn_all_plants$gall_total_abund
plot(vLG_dom ~ Genotype, tree_level_interaxn_all_plants)
rG_Tory_glm <- glm(rG_Tory ~ rG_dom, tree_level_interaxn_all_plants, family = "poisson")
summary(rG_Tory_glm)
plot(rG_Tory_glm)


### How does the density of different gall species vary among genotypes?

# vLG
plot(vLG_abund/shootEst.no18*200 ~ Genotype, tree_level_interaxn_all_plants)

vLG_density_null <- MASS::glm.nb(vLG_abund ~ offset(log(shootEst.no18)), tree_level_interaxn_all_plants)
vLG_density_glm <- MASS::glm.nb(vLG_abund ~ offset(log(shootEst.no18)) + Genotype, tree_level_interaxn_all_plants)
summary(vLG_density_glm)
anova(vLG_density_null, vLG_density_glm)
# plot(vLG_density_glm)
vLG_density_visreg <- visreg(vLG_density_glm, "Genotype", scale = "response", cond = list(shootEst.no18 = 200), ylim = c(0,40))


# rG
plot(rG_abund/shootEst.no18*200 ~ Genotype, tree_level_interaxn_all_plants)

rG_density_null <-  MASS::glm.nb(rG_abund ~ offset(log(shootEst.no18)), tree_level_interaxn_all_plants[-121,])
rG_density_glm <- MASS::glm.nb(rG_abund ~ offset(log(shootEst.no18)) + Genotype, tree_level_interaxn_all_plants[-121,])
summary(rG_density_glm)
plot(rG_density_glm)
anova(rG_density_null, rG_density_glm)

rG_density_visreg <- visreg(rG_density_glm, "Genotype", scale = "response", cond = list(shootEst.no18 = 200), ylim = c(0,120))

# rsLG
plot(rsLG_abund/shootEst.no18*200 ~ Genotype, tree_level_interaxn_all_plants)

rsLG_density_null <-  MASS::glm.nb(rsLG_abund ~ offset(log(shootEst.no18)), tree_level_interaxn_all_plants)
rsLG_density_glm <- MASS::glm.nb(rsLG_abund ~ offset(log(shootEst.no18)) + Genotype, tree_level_interaxn_all_plants)
summary(rsLG_density_glm)
plot(rsLG_density_glm)
anova(rsLG_density_null, rsLG_density_glm) 

rsLG_density_visreg <- visreg(rsLG_density_glm, "Genotype", scale = "response", cond = list(shootEst.no18 = 200), ylim = c(0,10))

# aSG

plot(aSG_abund/shootEst.no18*200 ~ Genotype, tree_level_interaxn_all_plants)

aSG_density_glm <- glm(aSG_abund ~ offset(log(shootEst.no18)) + Genotype, tree_level_interaxn_all_plants, family = "poisson") 
summary(aSG_density_glm) # no overdispersion in poisson model
plot(aSG_density_glm)
anova(aSG_density_glm, test = "Chi") 

aSG_density_visreg <- visreg(aSG_density_glm, "Genotype", scale = "response", cond = list(shootEst.no18 = 200), ylim = c(0,15))

# SG
plot(SG_abund/shootEst.no18*200 ~ Genotype, tree_level_interaxn_all_plants)

SG_density_glm <- glm(SG_abund ~ offset(log(shootEst.no18)) + Genotype, tree_level_interaxn_all_plants, family = "poisson") 
summary(SG_density_glm) # little overdispersion in poisson model. Same results with quasipoisson. iteration for theta doesn't converge for glm.nb so its likely not necessary.
plot(SG_density_glm)
anova(SG_density_glm, test = "Chi") 

SG_density_visreg <- visreg(SG_density_glm, "Genotype", scale = "response", cond = list(shootEst.no18 = 200), ylim = c(0,100)) # I think some of these density estimates are highly unrealistic...

# using mvabund package to test variation in gall abundances simultaneously. Takes less than a minute to run the anova model.
gall_community <- select(tree_level_interaxn_all_plants, aSG_abund:rsLG_abund)
gall_community_mvabund <- mvabund(gall_community)

gall_community_manyglm <- manyglm(gall_community_mvabund ~ offset(log(tree_level_interaxn_all_plants$shootEst.no18)) + tree_level_interaxn_all_plants$Genotype, family = "negative.binomial")
gall_community_manyglm_anova <- anova(gall_community_manyglm, p.uni = "unadjusted") # significant multivariate community response, driven primarily by vLG, although this response is only marginally significant using adjusted P-values.

# gall density summary information
gall_density_genotype_summary_df <- tree_level_interaxn_all_plants %>%
  group_by(Genotype) %>%
  select(aSG_abund:rsLG_abund, shootEst.no18) %>%
  summarise_each(funs(density200shoots = mean(./shootEst.no18*200))) %>% # note that rG's high estimate for S is driven by a single datapoint.
  mutate(vLG_dens_modelfit = vLG_density_visreg$Genotype$y$fit,
         rG_dens_modelfit = rG_density_visreg$Genotype$y$fit,
         rsLG_dens_modelfit = rsLG_density_visreg$Genotype$y$fit,
         aSG_dens_modelfit = aSG_density_visreg$Genotype$y$fit,
         SG_dens_modelfit = SG_density_visreg$Genotype$y$fit)

write.csv(gall_density_genotype_summary_df, "~/Documents/Genotype_Networks/data/gall_density_genotype_summary_df.csv")

### Does total gall density vary among genotypes?
gall_total_density_glm.nb_null <- MASS::glm.nb(gall_total_abund ~ offset(log(shootEst.no18)), tree_level_interaxn_all_plants[-c(121), ])
gall_total_density_glm.nb <- MASS::glm.nb(gall_total_abund ~ offset(log(shootEst.no18)) + Genotype, tree_level_interaxn_all_plants[-c(121),])

visreg(gall_total_density_glm.nb, "Genotype", scale = "response", cond = list(shootEst.no18 = 200))

summary(gall_total_density_glm.nb)
anova(gall_total_density_glm.nb_null, gall_total_density_glm.nb)
plot(gall_total_density_glm.nb)

### Does plant genotype affect link density?
plot(link_abund ~ Genotype, tree_level_interaxn_all_plants)
link_density_null <- MASS::glm.nb(link_abund ~ offset(log(shootEst.no18)), tree_level_interaxn_all_plants)
link_density_glm <- MASS::glm.nb(link_abund ~ offset(log(shootEst.no18)) + Genotype, tree_level_interaxn_all_plants)
summary(link_density_glm)
anova(link_density_null, link_density_glm)
plot(link_density_glm)

visreg(link_density_glm, "Genotype", scale = "response", ylim = c(0,15), cond = list(shootEst.no18 = 200))

### Does plant genotype affect link richness?
plot(link_richness ~ Genotype, tree_level_interaxn_all_plants)
linkrich.sumdf <- tree_level_interaxn_all_plants %>%
  group_by(Genotype) %>%
  summarise(mean.linkrich = mean(link_richness)) %>%
  arrange(mean.linkrich)

linkrich.df <- tree_level_interaxn_all_plants %>%
  select(Genotype, link_richness) %>%
  mutate(Genotype.ord = factor(Genotype, levels = linkrich.sumdf$Genotype, ordered = TRUE))

ggplot(linkrich.df, aes(x = Genotype.ord, y = link_richness)) + geom_boxplot()

link_rich_glm <- glm(link_richness ~ Genotype, tree_level_interaxn_all_plants, family = "poisson") # considering not using offset for richness as I should use rarefied information anyway.
summary(link_rich_glm)
anova(link_rich_glm, test = "Chi")
plot(link_rich_glm)

visreg(link_rich_glm, "Genotype", scale = "response", ylim = c(0,10)) # weird that the data points aren't exactly on each value...Oh wait, this is likely because it is richness density.


# Does the proportion of parasitized galls vary among genotypes?
test <- tree_level_interaxn_all_plants %>%
  group_by(Genotype) %>%
  summarise(link_abund_sum = sum(link_abund), gall_survive_abund_sum = sum(gall_survive_abund), gall_total_density_mean = mean(gall_total_densityNO18), vLG_density_mean = mean(vLG_abund/shootEst.no18))

test.glm <- glm(cbind(link_abund_sum, gall_survive_abund_sum) ~ vLG_density_mean, test[-c(19,24),], family = "binomial")
summary(test.glm)
plot(test.glm)
visreg(test.glm, scale = "response")


prop_parasitized_glm <- glm(cbind(link_abund, gall_survive_abund) ~ Genotype, filter(tree_level_interaxn_all_plants, Genotype != "U"), family = "binomial") # removed Genotype "U" because there was never any parasitoid attack...
summary(prop_parasitized_glm)
visreg(prop_parasitized_glm, "Genotype", scale = "response") # note that the points on the line of each "mean" are not real!!!!! Should figure out how to remove them if I want to plot this data.

anova(prop_parasitized_glm, test = "Chi")
plot(prop_parasitized_glm)

# Does plant genotype affect vLG density?
plot(vLG_abund/shootEst.no18 ~ Genotype, tree_level_interaxn_all_plants)
plot(vLG_abund ~ Genotype, tree_level_interaxn_all_plants)
vLG_density_null <- MASS::glm.nb(vLG_abund ~ offset(log(shootEst.no18)), tree_level_interaxn_all_plants)
vLG_density_glm <- MASS::glm.nb(vLG_abund ~ offset(log(shootEst.no18)) + Genotype, tree_level_interaxn_all_plants)
summary(vLG_density_glm)
anova(vLG_density_null, vLG_density_glm)
plot(vLG_density_glm)
test <- visreg(vLG_density_glm, "Genotype", scale = "response", ylim = c(0,40), cond = list(shootEst.no18 = 200))

round(test$Genotype$y$fit,1)

vLG_dataset <- tree_level_interaxn_all_plants %>%
  group_by(Genotype) %>%
  summarise(vLG_ptoid_count = sum(vLG_parasitized), vLG_surv_count = sum(vLG_vLG.pupa)) %>%
  mutate(vLG_genotype_mean_density = test$Genotype$y$fit)

vLG_geno.ptized = glm(cbind(vLG_ptoid_count, vLG_surv_count) ~ vLG_genotype_mean_density, vLG_dataset, family = "binomial")
summary(vLG_geno.ptized)
visreg(vLG_geno.ptized, scale = "response")
anova(vLG_geno.ptized, test = "Chi") # cool, looks like parasitoid attack rate increases with increasing gall density. Suggesting density dependence.
plot(vLG_geno.ptized)

