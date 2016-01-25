#####################################################################
# Random and Mixed effect models of plant genotype correlations with gall density and gall size, as well as how these factors shape gall survival.
#####################################################################

##### Upload required packages and functions for analysis

### packages

# visualization
library(visreg)
library(car)
library(ggplot2)

# analysis
library(psych)
library(RLRsim)
library(MASS)
library(mgcv)
library(gamm4)
library(nlme)
library(lme4)
library(pbkrtest)
library(vegan)

# data manipulation
library(reshape)
library(reshape2) 
library(plyr)
library(dplyr)

### functions

source('~/Documents/Genotype_Willow_Community/datasets_&_Rscripts/functions_ms_willow_community.R') # used for testing for random effects of plant genotype on gall densities

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

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

# function for calculating adjusted repeatability for nested random effect model. 
genotype_nested_model <- function(data, y){
  require(nlme)
  
  y <<- y
  # calculate broad-sense heritability of genotype
  random_geno_lme <- lme(y ~ 1, 
                         data, 
                         random = ~ 1 | Genotype/plant.position, 
                         method = "REML", 
                         na.action = na.omit) # use lme function because it is much easier to extract variance of random effects
  
  gen_var <- as.numeric(VarCorr(random_geno_lme)[2])
  tot_var <-  gen_var +  as.numeric(VarCorr(random_geno_lme)[5]) # not taking account into the variation explained by plant position. Therefore, this is an estimate of the adjusted repeatability or broad-sense heritability
  her <- gen_var / tot_var
  names(her) <- "Broad sense Heritability"
  print(her)
  
  # examining significance of random effect
  no_geno <- lm(y ~ 1, data) # no random effect to assess the effect of plant genotype. I tested this against a null model instead of a model with just plant.position as a random effect, because the variation due to plant position is primarily a result of plant genotype.
  rand_anova <- anova(random_geno_lme, no_geno)
  rand.Lratio <- rand_anova$L.Ratio[2]
  names(rand.Lratio) <- "Rand eff chi.sq"
  rand.pval <- rand_anova$"p-value"[2]
  names(rand.pval) <- "Rand eff p-val"
  print(rand_anova)
  
  # Return all objects in a list
  output = c(rand.Lratio, rand.pval, her)
}


###### Data management

### Upload data that is resolved to the tree level

# interaction data
tree_level_interaxn_all_plants <- read.csv("~/Documents/Genotype_Networks/data/tree_level_interaxn_all_plants.csv")
tree_level_interaxn_all_plants <- mutate(tree_level_interaxn_all_plants, plant.position = factor(plant.position))

# test analysis. don't know how promising the tree-level approach is going to be
sub.tree <- select(tree_level_interaxn_all_plants, plant.position, vLG_Tory, vLG_Platy, vLG_Mesopol)#aSG_aSG.larv:vLG_vLG.pupa, -aSG_aSG.larv, -rsLG_Pont.surv, -SG_SG.larv, - vLG_vLG.pupa, -rG_rG.larv)
row.names(sub.tree) <- sub.tree$plant.position
visweb(sub.tree[ ,-1], "diagonal")
library(bipartite)
sub.tree.mod <- computeModules(sub.tree[ ,-1])
plotModuleWeb(sub.tree.mod)
printoutModuleInformation(sub.tree.mod)

# plant trait data
tree_level_traits <- read.csv("~/Documents/Genotype_Networks/data/plant.trait.galls.2011.tree.df.csv")
tree_level_traits <- select(tree_level_traits, -X, -Genotype, plant.position, vLG.2011 = vLG, rsLG.2011 = rsLG, Total_Area:flavanonOLES.PC1)
tree_level_traits$plant.position <- factor(tree_level_traits$plant.position)

geno_level_traits <- read.csv("~/Documents/Genotype_Networks/data/plant.trait.galls.2011.genotype.df.csv")
geno_level_traits <- select(geno_level_traits, -X)

# merge interaction and tree-trait data
tree_level_interaxn_all_plants <- left_join(tree_level_interaxn_all_plants, tree_level_traits, by = "plant.position")

module.info <- read.csv('~/Documents/Genotype_Networks/data/module_info_genotype_gall_parasitoid_network.csv')


### prepare data that is resolved to the individual gall level

# upload molten gall network data
gall_net_melt <- read.csv("~/Documents/Genotype_Networks/data/gall_network_data.csv")
gall_net_melt <- tbl_df(gall_net_melt)

# filter data to focus on surviving and parasitized larva
gall_size_filter <- gall_net_melt %>%
  filter(gall_contents %in% c("aSG.larv", "Pont.ad", "Pont.prep", "rG.larv", "SG.larv", "vLG.pupa", "Eulo.fem", "Eulo.mal", "Eury.fem", "Eury.mal", "Lathro.fem", "Lathro.mal", "Lestodip", "Mesopol", "Mymarid", "Platy", "Ptero.2", "Tory.fem", "Tory.mal"))  

gall_size_niche_data <- gall_size_filter %>%
  group_by(plant.position, gall.id.nest) %>%
  summarise(g.height = mean(g.height)) %>% # focused on gall height, because exploring other measures of gall size (e.g. gall volume, reduced the predictive power, if anything)
  ungroup()

parasitoid_size_niche_data <- gall_size_filter %>%
  dcast(gall.id.nest ~ gall.sp + gall_contents, sum)

# merge data sets together and make some new columns of data
gall_parasitoid_size_niche_data <- left_join(gall_size_niche_data, parasitoid_size_niche_data)
gall_parasitoid_size_niche_data <- mutate(gall_parasitoid_size_niche_data, plant.position = factor(plant.position))
gall_parasitoid_size_niche_data$vLG_total <- rowSums(select(gall_parasitoid_size_niche_data, starts_with("vLG")))
gall_parasitoid_size_niche_data$rG_total <- rowSums(select(gall_parasitoid_size_niche_data, starts_with("rG")))
gall_parasitoid_size_niche_data$rsLG_total <- rowSums(select(gall_parasitoid_size_niche_data, starts_with("rsLG")))
gall_parasitoid_size_niche_data$aSG_total <- rowSums(select(gall_parasitoid_size_niche_data, starts_with("aSG")))
gall_parasitoid_size_niche_data$SG_total <- rowSums(select(gall_parasitoid_size_niche_data, starts_with("SG")))

# merge gall size data with tree-level data. This data set will be used for the mixed-effect models on gall survival as well as the association between plant genotype and gall size.
gall_size.df <- left_join(gall_parasitoid_size_niche_data, select(tree_level_interaxn_all_plants, Genotype, Gender, plant.position), by = "plant.position")
gall_size.df <- mutate(gall_size.df, plant.position = factor(plant.position))

# gall density data and plant architecture and trichome density. These were the only plant traits that I expected a priori to affects parasitoids. I would expect the other mechanisms to act indirectly by changing gall densities.
gall_density.df <- with(tree_level_interaxn_all_plants, 
                        data.frame(plant.position = plant.position, 
                                   Genotype = Genotype, 
                                   Gender = Gender, 
                                   vLG_density = vLG_abund/shootEst.no18,
                                   rG_density = rG_abund/shootEst.no18,
                                   aSG_density = aSG_abund/shootEst.no18,
                                   rsLG_density = rsLG_abund/shootEst.no18,
                                   SG_density = SG_abund/shootEst.no18,
                                   Total_Area = Total_Area,
                                   Density = Density,
                                   Height = Height,
                                   D = D_mean_smoothed,
                                   Trichomes = Trichome.No.))

# marge gall size and gall density data for mixed effects modelling
gall_mech.df <- left_join(gall_size.df, gall_density.df)


###### Random effect model of variation in gall community composition among willow genotypes. Uses tree-level data. I decided to use euclidean distance for 2 reasons. (1) It is directly interpretable in terms of "variation explained". And (2) using it allowed me to preserve all of the sample is my analyses, and (3) My analysis was focused on a small number of species that I thoroughly sampled on each tree, and that did not experience use beta-diversity.

gall.dens.tree.level <- with(tree_level_interaxn_all_plants, data.frame(vLG_dens = vLG_abund/shootEst.no18,
                                                                        rG_dens = rG_abund/shootEst.no18,
                                                                        rsLG_dens = rsLG_abund/shootEst.no18,
                                                                        SG_dens = SG_abund/shootEst.no18,
                                                                        aSG_dens = aSG_abund/shootEst.no18))

gall.adonis.euc <- adonis(sqrt(gall.dens.tree.level[-c(17,121), ]) ~ Genotype, tree_level_interaxn_all_plants[-c(17,121), ], method = "euclidean") # removed outliers for vLG-density. Note that log(x+0.01) transformation was slightly better, but I decided to stay consistent with the sqrt since that is what I used for all of the independent density models
gall.adonis.euc

gall.adonis.ss <- mean(table(tree_level_interaxn_all_plants$Genotype[-c(17,121)])) # 5.5 average sample size for genotypes
gall.gen.ms <- 0.035226
gall.res.ms <- 0.013430 # also equal to residual variance
gall.gen.var <- (gall.gen.ms - gall.res.ms)/gall.adonis.ss

gall.h2 <- gall.gen.var/(gall.gen.var + gall.res.ms)
gall.h2

anova(betadisper(vegdist(sqrt(gall.dens.tree.level[-c(17,121), ]), method = "euclidean"), tree_level_interaxn_all_plants$Genotype[-c(17,121)], bias.adjust = TRUE)) # not significant

#gall.meandist <- meandist(vegdist(sqrt(gall.dens.tree.level[-c(17,121), ]), method = "euclidean"), tree_level_interaxn_all_plants$Genotype)
#mean(diag(gall.meandist))
#mean(gall.meandist[which(lower.tri(gall.meandist, diag = FALSE))])

###### Random effect model of variation in gall-parasitoid network composition among willow genotypes. Uses tree-level data.

gall.ptoid.dens.tree.level <- with(tree_level_interaxn_all_plants, data.frame(aSG_Tory_dens = aSG_Tory/shootEst.no18,
                                                                        rG_Eulo_dens = rG_Eulo/shootEst.no18,
                                                                        rG_Lestodip_dens = rG_Lestodip/shootEst.no18,
                                                                        rG_Mesopol_dens = rG_Mesopol/shootEst.no18,
                                                                        rG_Platy_dens = rG_Platy/shootEst.no18,
                                                                        rG_Tory_dens = rG_Tory/shootEst.no18,
                                                                        rsLG_Lathro_dens = rsLG_Lathro/shootEst.no18,
                                                                        rsLG_Eury_dens = rsLG_Eury/shootEst.no18,
                                                                        SG_Platy_dens = SG_Platy/shootEst.no18,
                                                                        vLG_Eulo_dens = vLG_Eulo/shootEst.no18,
                                                                        vLG_Mesopol_dens = vLG_Mesopol/shootEst.no18,
                                                                        vLG_Mymarid_dens = vLG_Mymarid/shootEst.no18,
                                                                        vLG_Platy_dens = vLG_Platy/shootEst.no18,
                                                                        vLG_Tory_dens = vLG_Tory/shootEst.no18))

#gall.ptoid.tree.level <- select(tree_level_interaxn_all_plants, aSG_Tory, rG_Eulo, rG_Lestodip, rG_Mesopol, rG_Platy, rG_Tory, rsLG_Eury, rsLG_Lathro, SG_Platy, vLG_Eulo, vLG_Mesopol, vLG_Mymarid, vLG_Platy, vLG_Tory)
#sub <- which(rowSums(gall.ptoid.tree.level) > 0)

#ss <- table(tree_level_interaxn_all_plants$Genotype[sub])[-21] # several genotypes with only one sample. Note however, that adonis only retains groups with more than 1 sample.
#weights.1rep <- ss/sum(ss)
#mean(ss)
#weighted.mean(ss, weights.1rep)

#ss.2reps <- ss[which(ss > 1)] # sample size of genotypes with more than 1 estimate
#weights.2reps <- ss.2reps/sum(ss.2reps)
#mean(ss.2reps)
#weighted.mean(ss.2reps, weights.2reps)


gall.ptoid.adonis.euc <- adonis(sqrt(gall.ptoid.dens.tree.level[-c(17,121), ]) ~ Genotype, tree_level_interaxn_all_plants[-c(17,121), ], method = "euclidean") # same procedure as for gall community composition
gall.ptoid.adonis.euc

gall.ptoid.adonis.ss <- gall.adonis.ss # same as gall community composition samples

#gall.ptoid.adonis.euc.log <- adonis(log(gall.ptoid.tree.level+1) ~ shootEst.no18 + Genotype, tree_level_interaxn_all_plants, method = "euclidean")
#gall.ptoid.adonis.euc.log

#gall.ptoid.adonis.horn <- adonis(gall.ptoid.tree.level[sub, ] ~  Genotype, tree_level_interaxn_all_plants[sub, ], method = "horn")
#gall.ptoid.adonis.horn

gen.ms <- 0.0110087
res.ms <- 0.0050816
res.var <- res.ms
gen.var <- (gen.ms - res.ms)/gall.ptoid.adonis.ss

gall.ptoid.h2 <- gen.var/(gen.var + res.var)
gall.ptoid.h2


# result below suggest that networks among genotypes were 28% more dissimilar than network variation within genotypes.
#gall.ptoid.meandist.horn <- meandist(vegdist(gall.ptoid.tree.level[sub, ], "horn"), tree_level_interaxn_all_plants$Genotype[sub])
#mean(diag(gall.ptoid.meandist.horn), na.rm = TRUE) # 0.60
#sd(diag(gall.ptoid.meandist.horn), na.rm = TRUE) # 0.34

# note that the among network comparisons are made up of a much larger sample size.
#mean(gall.ptoid.meandist.horn[which(lower.tri(gall.ptoid.meandist.horn, diag = FALSE) == TRUE)]) # 0.77
#sd(gall.ptoid.meandist.horn[which(lower.tri(gall.ptoid.meandist.horn, diag = FALSE) == TRUE)]) # 0.20

gall.ptoid.beta.euc <- betadisper(vegdist(sqrt(gall.ptoid.dens.tree.level[-c(17,121), ]), "euclidean"), tree_level_interaxn_all_plants$Genotype[-c(17,121)], bias.adjust = TRUE)
anova(gall.ptoid.beta.euc) # a lot of beta-diversity among genotypes

#gall.ptoid.beta.euc.log <- betadisper(vegdist(log(gall.ptoid.tree.level+1), "euclidean"), tree_level_interaxn_all_plants$Genotype, bias.adjust = TRUE)
#plot(gall.ptoid.beta.euc.log)
#anova(gall.ptoid.beta.euc.log)

#gall.ptoid.beta.horn <- betadisper(vegdist(gall.ptoid.tree.level[sub, ], "horn"), tree_level_interaxn_all_plants$Genotype[sub], bias.adjust = TRUE)
#anova(gall.ptoid.beta.horn) # no longer dignificant dispersion among groups

###### Random effect models of variation in abundance (density) and richness of gall-parasitoid interactions among willow genotypes. Uses tree-level data

plot(link_richness ~ link_abund, tree_level_interaxn_all_plants)
plot(sqrt(link_richness) ~ sqrt(link_abund), tree_level_interaxn_all_plants)

# abundance of links

plot(link_abund/shootEst.no18 ~ Genotype, tree_level_interaxn_all_plants) # one major outlier
with(tree_level_interaxn_all_plants, which(link_abund/shootEst.no18 > 0.20)) # data point #121

with(tree_level_interaxn_all_plants[-121, ], hist(sqrt(link_abund/shootEst.no18))) # very skewed distribution

link_abund_lmer <- lmer(sqrt(link_abund/shootEst.no18) ~ (1 | Genotype), tree_level_interaxn_all_plants[-121, ])
summary(link_abund_lmer)
plot(link_abund_lmer)

exactRLRT(link_abund_lmer)

genotype_random_model(x = tree_level_interaxn_all_plants[-121, ], y = with(tree_level_interaxn_all_plants[-121, ], sqrt(link_abund/shootEst.no18))) # H2 = 0.39


# richness of links
plot(link_richness ~ Genotype, tree_level_interaxn_all_plants)
plot(link_richness/shootEst.no18 ~ Genotype, tree_level_interaxn_all_plants) # 121 is a major outlier again

with(tree_level_interaxn_all_plants, hist(sqrt(link_richness/shootEst.no18)))

# controls for number of shoots sampled
rich_abund_lmer <- lmer(sqrt(link_richness/shootEst.no18) ~ (1 | Genotype), tree_level_interaxn_all_plants[-121, ])
summary(rich_abund_lmer)
plot(rich_abund_lmer)

exactRLRT(rich_abund_lmer)

genotype_random_model(x = tree_level_interaxn_all_plants[-121, ], y = with(tree_level_interaxn_all_plants[-121, ], sqrt(link_richness/shootEst.no18))) # 0.39

# doesn't control for number of shoots sampled
rich_abund_lmer_nocontrol <- lmer(link_richness ~ (1 | Genotype), tree_level_interaxn_all_plants)
summary(rich_abund_lmer_nocontrol)
plot(rich_abund_lmer_nocontrol) # residuals look better without transformation

exactRLRT(rich_abund_lmer_nocontrol)
genotype_random_model(x = tree_level_interaxn_all_plants, y = tree_level_interaxn_all_plants$link_richness) # 0.40

#### rarefied number of links

gall.ptoid.tree.level <- with(tree_level_interaxn_all_plants, cbind.data.frame(Genotype,
                                                                               aSG_Tory,
                                                                               rG_Eulo,
                                                                               rG_Lestodip,
                                                                               rG_Mesopol,
                                                                               rG_Platy,
                                                                               rG_Tory,
                                                                               rsLG_Eury,
                                                                               rsLG_Lathro,
                                                                               SG_Platy,
                                                                               vLG_Eulo,
                                                                               vLG_Mesopol,
                                                                               vLG_Mymarid,
                                                                               vLG_Platy,
                                                                               vLG_Tory))
sub.gall.ptoid <- which(rowSums(gall.ptoid.tree.level[,-1]) > 1) # note that this drastically reduced my sample size.
rare.gall.ptoid <- rarefy(x = gall.ptoid.tree.level[ ,-1], sample = 2) # this allows natural 0 and 1 to stay in dataset.

rrich_gall_ptoid_lmer <- lmer(rare.gall.ptoid ~ (1|Genotype), gall.ptoid.tree.level)
summary(rrich_gall_ptoid_lmer)
plot(rrich_gall_ptoid_lmer) # look funky

exactRLRT(rrich_gall_ptoid_lmer) # not significant. Highly significant though if 0 and 1 are permitted in dataset
genotype_random_model(x = gall.ptoid.tree.level, y = rare.gall.ptoid)

# rarefied richness for subset of genotypes with sufficient samples
very.rare.gall.ptoid <- rarefy(x = gall.ptoid.tree.level[sub.gall.ptoid ,-1], sample = 2)
table(gall.ptoid.tree.level[sub.gall.ptoid, "Genotype"])
geno <- subset(gall.ptoid.tree.level[sub.gall.ptoid, ], Genotype %in% c("B","I","K","S","V","X","Z"))
very.rare.gall.ptoid.geno.sub <- rarefy(x = geno[ ,-1], sample = 2)
summary(lmer(very.rare.gall.ptoid.geno.sub ~ (1|Genotype), geno)) # no variance explained by genotype...

# Chao 1 estimated richness
chao1.rich <- estimateR(x = gall.ptoid.tree.level[ ,-1])
chao1.rich.only <- chao1.rich[2, ]

chao1_gall_ptoid_lmer <- lmer(chao1.rich.only ~ (1|Genotype), gall.ptoid.tree.level)
summary(chao1_gall_ptoid_lmer)
plot(chao1_gall_ptoid_lmer) # residuals look better

exactRLRT(chao1_gall_ptoid_lmer)
genotype_random_model(x = gall.ptoid.tree.level, y = chao1.rich.only) # H2 = 0.29, which seems much more reasonable.

plot(chao1.rich.only ~ chao1.rich[1,])

# richness while controlling for abundance in regression
summary(lm(link_richness ~ link_abund, tree_level_interaxn_all_plants))
summary(lm(link_richness ~ link_abund + I(link_abund^2), tree_level_interaxn_all_plants)) # better fit
plot(lm(link_richness ~ link_abund + I(link_abund^2), tree_level_interaxn_all_plants)) # residuals still look quite funky.

rich.abund.lmer <- lmer(link_richness ~ link_abund + I(link_abund^2) + (1|Genotype), tree_level_interaxn_all_plants)
summary(rich.abund.lmer)
exactRLRT(rich.abund.lmer) # still significant with link abund as linear predictor. No longer significant with link abund^2 as well


### Calculate rarefied richness at genotype-level
genotype.rarerich.df <- gall.ptoid.tree.level %>%
  group_by(Genotype) %>%
  summarise_each(funs(sum))
rownames(genotype.rarerich.df) <- genotype.rarerich.df$Genotype
genotype.rarerich.df$Richobs <- rowSums(genotype.rarerich.df[ ,-1] > 0)

write.csv(genotype.rarerich.df, "~/Documents/Genotype_Networks/data/link richness for estimateS.csv")

geno.sub <- which(rowSums(genotype.rarerich.df[ ,-1]) > 1)

geno.rarerich <- data.frame(Genotype = genotype.rarerich.df$Genotype[geno.sub], t(rarefy(genotype.rarerich.df[geno.sub,-1], 2, se = TRUE)))
geno.rarerich$CI <- geno.rarerich$se*1.96
geno.rarerich$Sprob <- geno.rarerich$S - 1


t(estimateR(genotype.rarerich.df[geno.sub,-1]))

t <- genotype.rarerich.df[ ,-1]
geno.label <- genotype.rarerich.df$Genotype[geno.sub]
rarecurve(x = genotype.rarerich.df[geno.sub, -1])

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


#### rsLG = Pontania californica

# box plots
rsLG_base <- ggplot(tree_level_interaxn_all_plants, aes(x = Genotype)) 
rsLG_base + geom_boxplot(aes(y = rsLG_abund)) 
rsLG_base + geom_boxplot(aes(y = rsLG_abund/shootEst.no18*100)) # no huge outliers

# data frame for rsLG analyses
rsLG.df <- tree_level_interaxn_all_plants

rsLG.df %>%
  group_by(Genotype) %>%
  summarise(mean.rsLG.dens = mean(rsLG_abund/shootEst.no18*100, na.rm = TRUE)) # up to 1.43 galls per 100 shoots.
1.42996530/0.04253578 # 33.6-fold

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
exactRLRT(rsLG_lmer)



#### aSG = Cecidomyiid sp. A (undescribed stem gall)

# box plots
aSG_base <- ggplot(tree_level_interaxn_all_plants, aes(x = Genotype)) 
aSG_base + geom_boxplot(aes(y = aSG_abund)) 
aSG_base + geom_boxplot(aes(y = aSG_abund/shootEst.no18*100)) # no huge outliers


# data frame for aSG analyses
aSG.df <- tree_level_interaxn_all_plants

aSG.df %>%
  group_by(Genotype) %>%
  summarise(mean.aSG.dens = mean(aSG_abund/shootEst.no18*100, na.rm = TRUE)) # up to 1.43 galls per 100 shoots.
0.73410214/0.03219616 # 22.8

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
exactRLRT(aSG_lmer) # better P-value estimate

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
exactRLRT(SG_lmer) # better p-value estimate


###### Phenotypic correlations among gall densities and tree phenotypes

# data frame of gall de+nsities
galls.traits <- select(tree_level_interaxn_all_plants, Genotype, plant.position, shootEst.no18, vLG_abund, rG_abund, rsLG_abund, aSG_abund, SG_abund, vLG.2011, rsLG.2011, Total_Area:flavanonOLES.PC1, Platy_abund:Mymarid_abund)
galls.traits <- mutate(galls.traits,
                vLG_density = vLG_abund/shootEst.no18,
                rG_density = rG_abund/shootEst.no18,
                rsLG_density = rsLG_abund/shootEst.no18,
                aSG_density = aSG_abund/shootEst.no18,
                SG_density = SG_abund/shootEst.no18,
                Platy_density = Platy_abund/shootEst.no18,
                Mesopol_density = Mesopol_abund/shootEst.no18,
                Tory_density = Tory_abund/shootEst.no18,
                Eulo_density = Eulo_abund/shootEst.no18,
                Eury_density = Eury_abund/shootEst.no18,
                Lathro_density = Lathro_abund/shootEst.no18,
                Mymarid_density = Mymarid_abund/shootEst.no18,
                Lestodip_density = Lestodip_abund/shootEst.no18)
galls.tree <- select(galls.traits, vLG_density, rG_density, rsLG_density, aSG_density, SG_density)
ptoids.tree <- select(galls.traits, Platy_density:Lestodip_density)
traits.tree <- select(galls.traits, Total_Area:flavanonOLES.PC1)

# correlations
galls.traits.corrs <- corr.test(y = sqrt(galls.tree), x = traits.tree, adjust = "none") # correlations didn't appear to change after removing influential vLG datapoints (17,121)
galls.traits.corrs$r

galls.ptoids.corrs <- corr.test(y = sqrt(ptoids.tree[-c(17,121),]), x = sqrt(galls.tree[-c(17,121),]), adjust = "none") # note that correlations are not qualitatively influenced by removing the outlying datapoints for vLG and rG_densities. So I have decided to just remove them.
galls.ptoids.corrs$r
round(galls.ptoids.corrs$r, 2)
round(galls.ptoids.corrs$p,3) # the density of all of the parasitoids is significantly positively correlated with each of their dominant host-species.
corr.test(y = sqrt(ptoids.tree[-c(17,121),]), x = sqrt(galls.tree[-c(17,121),]), adjust = "holm")

plot(Platy_density ~ vLG_density, galls.traits[-c(17,121),])

### vLG density mechanisms. Note that using a glm didn't improve the residual plots.
scatterplotMatrix(galls.traits[ ,c("vLG_density","myricitrin__A320nm","Height", "C_N_imputed", "D_mean_smoothed","N_imputed","Total_Area")])
plot(sqrt(vLG_density) ~ Genotype, galls.traits)
plot(sqrt(vLG_density) ~ myricitrin__A320nm, galls.traits)
plot(sqrt(vLG_density) ~ Height, galls.traits)
plot(sqrt(vLG_density) ~ C_N_imputed, galls.traits)

vLG.galls.traits <- galls.traits[-c(17,121), ]
vLG.galls.traits.comp <- vLG.galls.traits[complete.cases(vLG.galls.traits[ ,c("C_N_imputed","Height")]), c("vLG_abund", "shootEst.no18","C_N_imputed","Height") ]
which(rownames(vLG.galls.traits.comp) == 79) # point 67 is a potential outlier
which(rownames(vLG.galls.traits.comp) == 69) # 57
which(rownames(vLG.galls.traits.comp) == 99) #84
which(rownames(vLG.galls.traits.comp) == 124) #104

vLG_2012.bin <- glm(vLG_abund/shootEst.no18 ~ C_N_imputed + Height, 
                    data = vLG.galls.traits.comp,
                    family = quasibinomial,
                    weights = vLG.galls.traits.comp$shootEst.no18) 
summary(vLG_2012.bin) # results are robust to removal of outliers.
plot(vLG_2012.bin)
visreg(vLG_2012.bin, scale = "response") # very low probability of finding a gall for each shoot sampled, which makes sense. I feel like I can convert this to a probability per 100 shoots sampled.


vLG_2012.lm <- lm(sqrt(vLG_density) ~ C_N_imputed + Height, vLG.galls.traits) # myricitrin was also an important predictor
summary(vLG_2012.lm)
plot(vLG_2012.lm)

visreg(vLG_2012.lm, partial = FALSE)

ggplot(vLG.galls.traits, aes(x = Height, y = vLG_density)) + geom_point() + stat_smooth(method = lm)

### rG density mechanisms
plot(sqrt(rG_density) ~ C_N_imputed, galls.traits[-121,])
plot(sqrt(rG_density) ~ Height, galls.traits[-121,])

rG.galls.traits <- galls.traits[-c(121), ]

rG_2012.bin <- glm(rG_abund/shootEst.no18 ~ C_N_imputed + Height, rG.galls.traits, family = quasibinomial, weights = rG.galls.traits$shootEst.no18)
summary(rG_2012.bin)
plot(rG_2012.bin)
#visreg(rG_2012.bin, scale = "response") # need to scale per 100 shoots sampled.


rG_2012.lm <- lm(sqrt(rG_density) ~ C_N_imputed + Height, rG.galls.traits)
summary(rG_2012.lm) # same results as quasibinomial model.
plot(rG_2012.lm)


### rsLG density mechanisms
plot(sqrt(rsLG_density) ~ luteolin.der1__A320nm, galls.traits)
plot(sqrt(rsLG_density) ~ Total_Area, galls.traits)

rsLG_2012.bin <- glm(rsLG_abund/shootEst.no18 ~ luteolin.der1__A320nm + Total_Area, galls.traits, family = quasibinomial, weights = galls.traits$shootEst.no18)
summary(rsLG_2012.bin)
plot(rsLG_2012.bin)

rsLG_2012.lm <- lm(sqrt(rsLG_density) ~ flavonOLES.PC1 + log(Trichome.No.+1) + Total_Area, galls.traits) # luteolin.der1__A320nm
summary(rsLG_2012.lm)
plot(rsLG_2012.lm)
visreg(rsLG_2012.lm)

AIC(rsLG_2012.lm)

vif(rsLG_2012.lm)


### aSG density mechanisms
plot(sqrt(aSG_density) ~ eriodictyol.7.glucoside__A270nm, galls.traits)

aSG_2012.bin <- glm(aSG_abund/shootEst.no18 ~ eriodictyol.7.glucoside__A270nm,
                    galls.traits, 
                    family = quasibinomial, 
                    weights = galls.traits$shootEst.no18)
summary(aSG_2012.bin) # results appear robust to the removal of several datapoints with high leverage (22,95,98,99)
plot(aSG_2012.bin)

aSG_2012.lm <- lm(sqrt(aSG_density) ~ eriodictyol.7.glucoside__A270nm, galls.traits) # eriodictyol.7.glucoside__A270nm, flavanonOLES.PC1
summary(aSG_2012.lm)
plot(aSG_2012.lm)


###### Gall density plots among genotypes
vLG.dens <- galls.traits[-c(17,121),] %>%
  select(Genotype, plant.position, vLG_density) %>%
  melt()

rG.dens <-  galls.traits[-c(121),] %>%
  select(Genotype, plant.position, rG_density) %>%
  melt()

rsLG.dens <-  galls.traits %>%
  select(Genotype, plant.position, rsLG_density) %>%
  melt()

aSG.dens <- galls.traits %>%
  select(Genotype, plant.position, aSG_density) %>%
  melt()

gall_dens.df <- rbind_all(list(vLG.dens, rG.dens, rsLG.dens, aSG.dens))
gall_dens.df <- mutate(gall_dens.df, 
                       variable = factor(variable, levels = c("vLG_density","rG_density","rsLG_density","aSG_density"), ordered = TRUE),
                       Genotype = as.character(Genotype))
                       

mean_gall_dens.df <- gall_dens.df %>%
  group_by(Genotype, variable) %>%
  summarise(mean_dens = mean(value)) %>%
  arrange(variable, mean_dens, Genotype) %>%
  ungroup() %>%
  mutate(Genotype = factor(Genotype, levels = Genotype[1:26], ordered = TRUE))

(gall.dens.plot <- ggplot(mean_gall_dens.df, aes(x = Genotype, y = mean_dens*100, color = variable, group = variable)) + 
  geom_line(aes(linetype = variable), size = 1) +
  geom_point(aes(shape = variable), size = 4) +
  theme_classic() +
  theme(legend.position = "top") +
  theme(axis.text = element_text(size = 18),
         axis.title = element_text(size = 20),
         axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 1),
        legend.text = element_text(size = 12, face = "italic")) + 
  scale_color_manual(values = c("#000000","#333333","#666666","#CCCCCC"),
                     name = "",
                     breaks = c("vLG_density","rG_density","rsLG_density","aSG_density"),
                     labels = c("Iteomyia","Rabdophaga-B","Pontania","Cecidomyiid")) +
  scale_shape_manual(values = c(16,17,15,8),
                     name = "",
                     breaks = c("vLG_density","rG_density","rsLG_density","aSG_density"),
                     labels = c("Iteomyia","Rabdophaga-B","Pontania","Cecidomyiid")) +
   scale_linetype_manual(values = c(1:4),
                         name = "",
                         breaks = c("vLG_density","rG_density","rsLG_density","aSG_density"),
                         labels = c("Iteomyia","Rabdophaga-B","Pontania","Cecidomyiid")) +
  ylab("No. galls per 100 shoots") +
  xlab("") +
   annotate("text", label = "(A)", size = 8, x = 1.5, y = 4.5)) # omitting xlab because it will be above another plot with the same label



###### Gather data for genetic correlations. According to Cheverud 1988, it may be best to use phenotypic correlations, rather than genetic ones... Note that many of the phenotypic correlations I detected are not retained in the genetic data...

# removed these data points because they bias the density estimates for vLG and rG, the dominants gall species. But note that they should only be removed for the galls and note the traits. Also, C, N, and C:N should be estimated from the original mean estimates. Actually, that should be the case for all of the traits, instead of the restricted sample sizes I have.
galls.ptoids.geno <- galls.traits[-c(17,121), ] %>%
  group_by(Genotype) %>%
  summarise_each(funs(mean(., na.rm = TRUE))) %>%
  select(Genotype, vLG_density:Lestodip_density)

galls.geno <- select(galls.ptoids.geno, Genotype, vLG_density:aSG_density)
traits.geno <- select(geno_level_traits, Genotype, Total_Area:flavanonOLES.PC1)

geno.corrs <- corr.test(y = galls.geno[ ,-1], x = traits.geno[ ,-1], adjust = "none")
round(geno.corrs$r, 2)


### gall correlations
gall.geno.corrs <- corr.test(galls.geno)
gall.geno.corrs$r
gall.geno.corrs$p

plot(vLG_density ~ rG_density, galls.geno)
plot(vLG_density ~ rsLG_density, galls.geno)

### vLG density geno correlations
plot(sqrt(vLG_density) ~ C_N_imputed, galls.traits.geno)
plot(vLG_density ~ sal_tannin.PC1, galls.traits.geno)
plot(vLG_density ~ D_mean_smoothed, galls.traits.geno)
plot(vLG_density ~ Height, galls.traits.geno)
plot(vLG_density ~ ampelopsin.der__A320nm, galls.traits.geno) # not promising
plot(vLG_density ~ eriodictyol.7.glucoside__A270nm, galls.traits.geno) # not promising


vLG_geno.lm <- lm(vLG_density ~ D_mean_smoothed, galls.traits.geno)
summary(vLG_geno.lm)
plot(vLG_geno.lm)

### rG density geno correlations
plot(rG_density ~ D_mean_smoothed, galls.traits.geno)

rG_geno.lm <- lm(rG_density ~ D_mean_smoothed, galls.traits.geno)
summary(rG_geno.lm)
plot(rG_geno.lm)


### rsLG density geno correlations
plot(rsLG_density ~ Total_Area, galls.traits.geno)
plot(rsLG_density ~ luteolin.der1__A320nm, galls.traits.geno)
plot(rsLG_density ~ water_content, galls.traits.geno)

rsLG_geno.lm <- lm(rsLG_density ~  water_content + luteolin.der1__A320nm, galls.traits.geno) # Total_Area, water_content, luteolin.der1 all good ones
summary(rsLG_geno.lm)
plot(rsLG_geno.lm)


### aSG density geno correlations
plot(aSG_density ~ eriodictyol.7.glucoside__A270nm, galls.traits.geno)
plot(aSG_density ~ water_content, galls.traits.geno) # looks more real

aSG_geno.lm <- lm(aSG_density ~ water_content, galls.traits.geno)
summary(aSG_geno.lm)
plot(aSG_geno.lm)

###### Nested random effect models of variation in gall size among genotypes. Note that calculations are based on adjusted repeatabilities, which do not include the variance explained by different plant positions.


### vLG
vLG.size.df <- gall_size.df %>%
  filter(vLG_total > 0)

geno.levels <- c("J","D","A","*","X","O","Y","V","Z","L","I","P","K","M","E","T","H","W","R","B","S","N","F","G") 

mean.vLG.size.by.tree <- vLG.size.df %>%
  group_by(plant.position, Genotype) %>%
  summarise(vLG.g.height = mean(g.height, na.rm = TRUE), N = n()) %>%
  ungroup() %>%
  mutate(Genotype.ord = factor(as.character(Genotype), levels = geno.levels, ordered = TRUE))

mean.vLG.size.geno.df <- mean.vLG.size.by.tree %>%
  group_by(Genotype) %>%
  summarise(vLG.mean.g.height = mean(vLG.g.height, na.rm = TRUE), ss = n(), se = sd(vLG.g.height, na.rm = TRUE)/sqrt(ss)) %>%
  arrange(vLG.mean.g.height) %>%
  ungroup() %>%
  mutate(Genotype = factor(Genotype, levels = Genotype, ordered = TRUE))

#levels(mean.vLG.size.geno.df$Genotype)

11.008636/4.830000 # 2.3-fold

(vLG.size.plot <- ggplot(mean.vLG.size.by.tree, aes(x = Genotype.ord, y = vLG.g.height)) +
  geom_point(color = "grey", size = 4, shape = 21) + 
  stat_summary(fun.y = mean, color = "black", size = 10, geom = "point", shape = "-") +
  theme_classic() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.title.x = element_text(vjust = -0.25)) + 
  xlab("Willow Genotype") +
  ylab(expression(paste(italic("Iteomyia")," Gall Size (mm)"))) +
   annotate("text", label = "(B)", size = 8, x = 1.5, y = 14))

require(gridExtra)

grid.arrange(gall.dens.plot, vLG.size.plot, ncol = 1)

#ggplot(mean.vLG.size.geno.df, aes(x = Genotype, y = vLG.mean.g.height)) +
 # geom_point(aes(size = ss)) +
  #geom_errorbar(aes(ymin = vLG.mean.g.height - se, ymax = vLG.mean.g.height + se), width = 0.2) +
  #theme_classic() +
  #ylab("Iteomyia Gall Size (mm)")

#vLG.size.df.special <- mutate(vLG.size.df, Genotype.ord = factor(as.character(Genotype), levels = levels(mean.vLG.size.geno.df$Genotype), ordered = TRUE))

pp.mean <- function(x){
    x %>%
    group_by(Genotype, plant.position) %>%
    summarise(mean = mean(x$g.height, na.rm = TRUE)) 
}
gen.mean <- function(x){
  x %>%
    group_by(Genotype, plant.position) %>%
    summarise(mean.size = mean(x$g.height, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(Genotype) %>%
    summarise(mean = mean(mean.size, na.rm = TRUE))
}

#p <- ggplot(data = vLG.size.df.special, aes(x = Genotype.ord, y = g.height, group = Genotype:plant.position)) +
 # geom_point(color = "grey") +
  #theme_classic()
#p + stat_summary(fun.y = mean, geom = "point", color = "red", size = 3)

# boxplots
ggplot(data = vLG.size.df.special, aes(x = Genotype.ord, y = g.height)) + geom_boxplot() + theme_classic()
ggplot(data = vLG.size.df, aes(x = plant.position, y = g.height)) + geom_boxplot() # based on heterogeneity, it will be important to take into account plant position

# nested random effect model
vLG.size.lmer <- lmer(sqrt(g.height) ~ 1 + (1 | Genotype) + (1 |plant.position:Genotype),
                      data = vLG.size.df)
summary(vLG.size.lmer)
vLG.size.lmer.gen <- update(vLG.size.lmer, .~. - (1 | plant.position:Genotype))
vLG.size.lmer.pp <- update(vLG.size.lmer, .~. - (1 | Genotype))

(0.0184)/(0.0184+0.1181) # H2 = 0.13
exactRLRT(m = vLG.size.lmer.gen, mA = vLG.size.lmer, m0 = vLG.size.lmer.pp)
#PBmodcomp(vLG.size.lmer, vLG.size.lmer.pp)

#vLG.complete <- vLG.size.df[complete.cases(vLG.size.df[ ,c("g.height","Genotype","plant.position")]), ]


#genotype_nested_model(data = vLG.size.df, y = sqrt(vLG.size.df$g.height)) # 0.13, p < 0.001


ggplot(data = data.frame(resid = residuals(vLG.size.lmer, type = "pearson"), fit = fitted(vLG.size.lmer)), aes(x = fit, y = resid)) + 
  geom_point() +
  geom_hline(y = 0) +
  stat_smooth() # some bias and heteroscedasticity
vLG.size.geno_ranef <- ranef(vLG.size.lmer)$Genotype$"(Intercept)"
vLG.size.pp_ranef <- ranef(vLG.size.lmer)$"plant.position:Genotype"$"(Intercept)"
ggQQ_ranef(vLG.size.geno_ranef) # looks great
ggQQ_ranef(vLG.size.pp_ranef) # looks like there is one outlier


### rG
rG.size.df <- gall_size.df %>%
  filter(rG_total > 0)

# boxplots
ggplot(data = rG.size.df, aes(x = Genotype, y = g.height)) + geom_boxplot()
ggplot(data = rG.size.df, aes(x = plant.position, y = g.height)) + geom_boxplot() # based on heterogeneity, it will be important to take into account plant position

# nested random effect model
rG.size.lmer <- lmer(sqrt(g.height) ~ 1 + (1 | Genotype/plant.position),
                     data = rG.size.df)
summary(rG.size.lmer)

rG.size.lmer.gen <- update(rG.size.lmer, .~. - (1 | plant.position:Genotype))
rG.size.lmer.pp <- update(rG.size.lmer, .~. - (1 | Genotype))

exactRLRT(m = rG.size.lmer.gen, mA = rG.size.lmer, m0 = rG.size.lmer.pp)
0.002969/(0.002969+0.06519) # H2 = 0.04
#genotype_nested_model(data = rG.size.df, y = (rG.size.df$g.height)) # 0.04, p = 0.007

ggplot(data = data.frame(resid = residuals(rG.size.lmer, type = "pearson"), fit = fitted(rG.size.lmer)), aes(x = fit, y = resid)) + 
  geom_point() +
  geom_hline(y = 0) +
  stat_smooth() # bias and heteroscedasticity
rG.size.geno_ranef <- ranef(rG.size.lmer)$Genotype$"(Intercept)"
rG.size.pp_ranef <- ranef(rG.size.lmer)$"plant.position:Genotype"$"(Intercept)"
ggQQ_ranef(rG.size.geno_ranef) # looks great
ggQQ_ranef(rG.size.pp_ranef) # looks like there is one outlier

### rsLG
rsLG.size.df <- gall_size.df %>%
  filter(rsLG_total > 0)

# boxplots
ggplot(data = rsLG.size.df, aes(x = Genotype, y = g.height)) + geom_boxplot()
ggplot(data = rsLG.size.df, aes(x = plant.position, y = g.height)) + geom_boxplot() # based on heterogeneity, it will be important to take into account plant position

# nested random effect model
rsLG.size.lmer <- lmer(sqrt(g.height) ~ 1 + (1 | Genotype) + (1 | plant.position:Genotype),
                       data = rsLG.size.df)
summary(rsLG.size.lmer)

rsLG.size.lmer.gen <- update(rsLG.size.lmer, .~. - (1 | plant.position:Genotype))
rsLG.size.lmer.pp <- update(rsLG.size.lmer, .~. - (1 | Genotype))

(0.004115)/(0.004115+0.044613) # H2 = 0.08
exactRLRT(m = rsLG.size.lmer.gen, mA = rsLG.size.lmer, m0 = rsLG.size.lmer.pp)


ggplot(data = data.frame(resid = residuals(rsLG.size.lmer, type = "pearson"), fit = fitted(rsLG.size.lmer)), aes(x = fit, y = resid)) + 
  geom_point() +
  geom_hline(y = 0) +
  stat_smooth() # some bias and heteroscedasticity
rsLG.size.geno_ranef <- ranef(rsLG.size.lmer)$Genotype$"(Intercept)"
rsLG.size.pp_ranef <- ranef(rsLG.size.lmer)$"plant.position:Genotype"$"(Intercept)"
ggQQ_ranef(rsLG.size.geno_ranef) # looks okay
ggQQ_ranef(rsLG.size.pp_ranef) # looks pretty good


### aSG. P-value doesn't seem reliable.
aSG.size.df <- gall_size.df %>%
  filter(aSG_total > 0)

# boxplots
ggplot(data = aSG.size.df, aes(x = Genotype, y = g.height)) + geom_boxplot()
ggplot(data = aSG.size.df, aes(x = plant.position, y = g.height)) + geom_boxplot() # based on heterogeneity, it will be important to take into account plant position

# nested random effect model
aSG.size.lmer <- lmer(sqrt(g.height) ~ 1 + (1 | Genotype) + (1 | plant.position:Genotype),
                      data = aSG.size.df)
summary(aSG.size.lmer)

aSG.size.lmer.gen <- update(aSG.size.lmer, .~. - (1 | plant.position:Genotype))
aSG.size.lmer.pp <- update(aSG.size.lmer, .~. - (1 | Genotype))

# H2 < 0.01
exactRLRT(m = aSG.size.lmer.gen, mA = aSG.size.lmer, m0 = aSG.size.lmer.pp)


#genotype_nested_model(data = aSG.size.df, y = sqrt(aSG.size.df$g.height)) # H2 < 0.001, p = 0.0852. Weird...

ggplot(data = data.frame(resid = residuals(aSG.size.lmer, type = "pearson"), fit = fitted(aSG.size.lmer)), aes(x = fit, y = resid)) + 
  geom_point() +
  geom_hline(y = 0) +
  stat_smooth() # some bias and heteroscedasticity
aSG.size.geno_ranef <- ranef(aSG.size.lmer)$Genotype$"(Intercept)"
aSG.size.pp_ranef <- ranef(aSG.size.lmer)$"plant.position:Genotype"$"(Intercept)"
ggQQ_ranef(aSG.size.geno_ranef) # some skew.
ggQQ_ranef(aSG.size.pp_ranef) # some definite outliers


### SG. H2-value doesn't seem reliable. Note that sample sizes were incredibly small and therefore likely not useful.
SG.size.df <- gall_size.df %>%
  filter(SG_total > 0)

# boxplots
ggplot(data = SG.size.df, aes(x = Genotype, y = g.height)) + geom_boxplot()
ggplot(data = SG.size.df, aes(x = plant.position, y = g.height)) + geom_boxplot() # based on heterogeneity, it will be important to take into account plant position

# nested random effect model
SG.size.lmer <- lmer(sqrt(g.height) ~ 1 + (1 | Genotype/plant.position),
                     data = SG.size.df)
summary(SG.size.lmer)

SG.size.lmer.gen <- update(SG.size.lmer, .~. - (1 | plant.position:Genotype))
SG.size.lmer.pp <- update(SG.size.lmer, .~. - (1 | Genotype))

0.6607/(0.6607+0.253) # H2 = 0.723
exactRLRT(m = SG.size.lmer.gen, mA = SG.size.lmer, m0 = SG.size.lmer.pp)


ggplot(data = data.frame(resid = residuals(SG.size.lmer, type = "pearson"), fit = fitted(SG.size.lmer)), aes(x = fit, y = resid)) + 
  geom_point() +
  geom_hline(y = 0) +
  stat_smooth() # some bias and heteroscedasticity
SG.size.geno_ranef <- ranef(SG.size.lmer)$Genotype$"(Intercept)"
SG.size.pp_ranef <- ranef(SG.size.lmer)$"plant.position:Genotype"$"(Intercept)"
ggQQ_ranef(SG.size.geno_ranef) # normal
ggQQ_ranef(SG.size.pp_ranef) # normal


###### Which plant traits predict variation in mean gall size among trees?

vLG.size.traits <- left_join(ungroup(mean.vLG.size.by.tree), galls.traits, by = "plant.position")

size.df <- select(vLG.size.traits, vLG.g.height, N, vLG_density)
size.trait.df <- select(vLG.size.traits, Total_Area:flavanonOLES.PC1)

vLG.size.corrs <- corr.test(y = size.df, x = size.trait.df, adjust = "none")
vLG.size.corrs$r # interesting that total area and D smoothed are now better predictors 


#vLG.size.chems <- select(vLG.size.traits, vLG.g.height, N, sal_tannin.PC1:flavanonOLES.PC1)

#vLG.size.traits.df.sub <- vLG.size.traits[complete.cases(vLG.size.traits[ ,c("vLG.g.height","salicortin__A270nm", "luteolin.der1__A320nm")]), ]

plot(vLG.g.height ~ luteolin.der1__A320nm, vLG.size.traits)
plot(vLG.g.height ~ salicortin__A270nm, vLG.size.traits)

vLG.size.trait.lm <- lm(vLG.g.height ~ sal_tannin.PC1 + flavonOLES.PC1, vLG.size.traits, weights = vLG.size.traits$N) # weighted the observations by sample size, because I have a better estimate of mean gall size of trees with larger sample size.
summary(vLG.size.trait.lm)
plot(vLG.size.trait.lm)
visreg(vLG.size.trait.lm)

AIC(vLG.size.trait.lm)

vLG.size.traits.chem.sub <- vLG.size.traits[complete.cases(vLG.size.traits[ ,c("sal_tannin.PC1", "cinn.PC1", "cinn.PC2", "flavonOLES.PC1", "flavonOLES.PC2", "flavanonOLES.PC1")]), ]
n <- lm(vLG.g.height ~ 1, vLG.size.traits.chem.sub, weights = vLG.size.traits.chem.sub$N)
f <- lm(vLG.g.height ~ sal_tannin.PC1 + cinn.PC1 + cinn.PC2 + flavonOLES.PC1 + flavonOLES.PC2 + flavanonOLES.PC1, vLG.size.traits.chem.sub, weights = vLG.size.traits.chem.sub$N)
summary(f)
test <- stepAIC(object = n, scope = formula(f), direction = "forward")
summary(test)

###### Nested random effect models of variation in gall survival among genotypes. Note that using the glmerControl(optimizer = "bobyqa") does not result in any model convergence problems and that that likelihood ratio test statistics are the same as other optimizer.

### vLG survival

# models
vLG_surv <- glmer(cbind(vLG_vLG.pupa, vLG_total - vLG_vLG.pupa) ~ 1 + (1|Genotype) + (1 | plant.position:Genotype), data = filter(gall_size.df, vLG_total > 0), family = binomial, control=glmerControl(optimizer="bobyqa"))
summary(vLG_surv)
ranef(vLG_surv)

overdisp_fun(vLG_surv) # note that the overdispersion ratio almost exactly matches that from the glmmPQL model's residual variance estimate. Therefore, rptR does calculate the residual variance correctly

vLG_surv.red <- glmer(cbind(vLG_vLG.pupa, vLG_total - vLG_vLG.pupa) ~ 1 + (1 | plant.position:Genotype), data = filter(gall_size.df, vLG_total > 0), family = binomial, control=glmerControl(optimizer="bobyqa"))


vLG_surv.PQL <- glmmPQL(cbind(vLG_vLG.pupa, vLG_total - vLG_vLG.pupa) ~ 1,
                        random = ~ 1|Genotype/plant.position,
                        data = filter(gall_size.df, vLG_total > 0), 
                        family = quasibinomial)
summary(vLG_surv.PQL)
0.6283202^2/(0.6283202^2 + (1.156049)*(pi^2)/3)

#PBmodcomp(vLG_surv, vLG_surv.red) # barely significant...

# residuals
plot(vLG_surv)


### rG survival

# models
rG_surv <- glmer(cbind(rG_rG.larv, rG_total - rG_rG.larv) ~ 1 + (1|Genotype) + (1 | plant.position:Genotype), data = filter(gall_size.df, rG_total > 0), family = binomial, control=glmerControl(optimizer="bobyqa"))
summary(rG_surv) # note that zero variance was explained, so testing whether the variance explained is greater than zero doesn't make sense.

rG_surv.red <- glmer(cbind(rG_rG.larv, rG_total - rG_rG.larv) ~ 1 + (1 | plant.position:Genotype), data = filter(gall_size.df, rG_total > 0), family = binomial, control=glmerControl(optimizer="bobyqa"))

# residuals
plot(rG_surv)

### vLG_Platy
vLG_Platy.PQL <- glmmPQL(cbind(vLG_Platy, vLG_total - vLG_Platy) ~ 1,
                        random = ~ 1|Genotype/plant.position,
                        data = filter(gall_size.df, vLG_total > 0), 
                        family = quasibinomial)
summary(vLG_Platy.PQL)
plot(vLG_Platy.PQL)
0.8783298^2/(0.8783298^2 + 1.013257^2*pi^2/3)

# vLG_Mesopol
vLG_Mesopol.PQL <- glmmPQL(cbind(vLG_Mesopol + vLG_Ptero.2, vLG_total - vLG_Mesopol - vLG_Ptero.2) ~ 1,
                         random = ~ 1|Genotype/plant.position,
                         data = filter(gall_size.df, vLG_total > 0), 
                         family = quasibinomial)
summary(vLG_Mesopol.PQL)
plot(vLG_Mesopol.PQL)
0.0006728894^2/(0.0006728894^2 + 0.849551^2*pi^2/3) # very small heritability...

# vLG_Ecto
vLG_ecto.PQL <- glmmPQL(cbind(vLG_Mesopol + vLG_Ptero.2 + vLG_Tory.mal + vLG_Tory.fem + vLG_Eulo.fem + vLG_Eulo.mal, vLG_total - vLG_Mesopol - vLG_Ptero.2 - vLG_Tory.mal - vLG_Tory.fem - vLG_Eulo.fem - vLG_Eulo.mal) ~ 1,
                        random = ~ 1|Genotype/plant.position,
                           data = filter(gall_size.df, vLG_total > 0), 
                           family = quasibinomial)
summary(vLG_ecto.PQL)

0.003322889^2/(0.003322889^2 + 1.061833^2*pi*2/3) # < 0.001


###### Random effect models of variation in gall-ptoid attack among genotypes

colSums(select(tree_level_interaxn_all_plants, aSG_aSG.larv:vLG_vLG.pupa))

# note that vLG_Platy, vLG_Tory, rG_Tory, and vLG_Mesopol make up 76% of the trophic interactions in this system.
sum(colSums(select(tree_level_interaxn_all_plants, aSG_Tory,rG_Eulo,rG_Lestodip, rG_Mesopol, rG_Platy, rG_Tory, rsLG_Eury, rsLG_Lathro, SG_Platy, vLG_Eulo, vLG_Mesopol, vLG_Mymarid, vLG_Platy, vLG_Tory)))

colSums(select(tree_level_interaxn_all_plants, Platy_abund:Mymarid_abund))
(97+75+51)/265 # Platy, Tory, and Mesopol make up 84% of parasitoid abundance.

### vLG_Platy density. Note that I tried a quasipoisson glmmPQL and the residuals still don't look great, so I've decided to go with the more straight-forward linear models.

# boxplots
plot(vLG_Platy/shootEst.no18*100 ~ Genotype, tree_level_interaxn_all_plants) # one outlier #121

tree_level_interaxn_all_plants[-121, ] %>%
  group_by(Genotype) %>%
  summarise(vLG_Platy_dens = mean(vLG_Platy/shootEst.no18*100, na.rm = T))
1.90521520/0.05460694

# analyses
vLG.platy.lmer <- lmer(sqrt(vLG_Platy/shootEst.no18) ~ (1 | Genotype), data = tree_level_interaxn_all_plants[-121, ])
summary(vLG.platy.lmer)
plot(vLG.platy.lmer)

vLG.platy.lmer.ranef <- ranef(vLG.platy.lmer)$Genotype$"(Intercept)"
ggQQ_ranef(vLG.platy.lmer.ranef)

genotype_random_model(x = tree_level_interaxn_all_plants[-121,], y = with(tree_level_interaxn_all_plants[-121,], sqrt(vLG_Platy/shootEst.no18)))
exactRLRT(vLG.platy.lmer)


### vLG_Mesopol density

# boxplots
plot(vLG_Mesopol/shootEst.no18*100 ~ Genotype, tree_level_interaxn_all_plants) 

tree_level_interaxn_all_plants %>%
  group_by(Genotype) %>%
  summarise(vLG_Mesopol_dens = mean(vLG_Mesopol/shootEst.no18*100, na.rm = T))
0.54467040/0.05164559 # 10.5

# analyses
vLG_Mesopol.lmer <- lmer(sqrt(vLG_Mesopol/shootEst.no18) ~ (1 | Genotype), data = tree_level_interaxn_all_plants)
summary(vLG_Mesopol.lmer)
plot(vLG_Mesopol.lmer)

vLG_Mesopol.lmer.ranef <- ranef(vLG_Mesopol.lmer)$Genotype$"(Intercept)"
ggQQ_ranef(vLG_Mesopol.lmer.ranef) # looks pretty good...

genotype_random_model(x = tree_level_interaxn_all_plants, y = with(tree_level_interaxn_all_plants, sqrt(vLG_Mesopol/shootEst.no18)))
exactRLRT(vLG_Mesopol.lmer)


### vLG_Tory density

# boxplots
plot(vLG_Tory/shootEst.no18*100 ~ Genotype, tree_level_interaxn_all_plants)

tree_level_interaxn_all_plants %>%
  group_by(Genotype) %>%
  summarise(vLG_Tory_dens = mean(vLG_Tory/shootEst.no18*100, na.rm = T))
0.30377293/0.05345305 # 5.7

# analyses
vLG_Tory.lmer <- lmer(sqrt(vLG_Tory/shootEst.no18) ~ (1 | Genotype), data = tree_level_interaxn_all_plants)
summary(vLG_Tory.lmer)
plot(vLG_Tory.lmer)

vLG_Tory.lmer.ranef <- ranef(vLG_Tory.lmer)$Genotype$"(Intercept)"
ggQQ_ranef(vLG_Tory.lmer.ranef) # looks pretty good...

genotype_random_model(x = tree_level_interaxn_all_plants, y = with(tree_level_interaxn_all_plants, sqrt(vLG_Tory/shootEst.no18)))
exactRLRT(vLG_Tory.lmer)


### rG_Tory density. No clear response in rG_Tory density

# boxplots
plot(rG_Tory/shootEst.no18*100 ~ Genotype, tree_level_interaxn_all_plants[-121,]) # one major outlier 

# analyses
rG_Tory.lmer <- lmer(sqrt(rG_Tory/shootEst.no18) ~ (1 | Genotype), data = tree_level_interaxn_all_plants)
summary(rG_Tory.lmer)
plot(rG_Tory.lmer)

rG_Tory.lmer.ranef <- ranef(rG_Tory.lmer)$Genotype$"(Intercept)"
ggQQ_ranef(rG_Tory.lmer.ranef) # looks pretty good...

genotype_random_model(x = tree_level_interaxn_all_plants, y = with(tree_level_interaxn_all_plants, sqrt(rG_Tory/shootEst.no18)))
exactRLRT(rG_Tory.lmer)

# plot of correlation between vLG_Tory and rG_Tory. Does this 
plot(I(rG_Tory/shootEst.no18) ~ I(vLG_Tory/shootEst.no18), tree_level_interaxn_all_plants[-121, ])


### vLG_Eulo density. Not enough individuals to probably even detect an effect if it existed.

# boxplots
plot(vLG_Eulo/shootEst.no18*100 ~ Genotype, tree_level_interaxn_all_plants) 

# analyses
vLG_Eulo.lmer <- lmer(sqrt(vLG_Eulo/shootEst.no18) ~ (1 | Genotype), data = tree_level_interaxn_all_plants)
summary(vLG_Eulo.lmer)
plot(vLG_Eulo.lmer)

vLG_Eulo.lmer.ranef <- ranef(vLG_Eulo.lmer)$Genotype$"(Intercept)"
ggQQ_ranef(vLG_Eulo.lmer.ranef) # looks pretty good...

genotype_random_model(x = tree_level_interaxn_all_plants, y = with(tree_level_interaxn_all_plants, sqrt(vLG_Eulo/shootEst.no18)))
exactRLRT(vLG_Eulo.lmer)

### rsLG_Eury density

# boxplots
plot(rsLG_Eury/shootEst.no18*100 ~ Genotype, tree_level_interaxn_all_plants) 

# analyses
rsLG_Eury.lmer <- lmer(sqrt(rsLG_Eury/shootEst.no18) ~ (1 | Genotype), data = tree_level_interaxn_all_plants)
summary(rsLG_Eury.lmer)
plot(rsLG_Eury.lmer)

rsLG_Eury.lmer.ranef <- ranef(rsLG_Eury.lmer)$Genotype$"(Intercept)"
ggQQ_ranef(rsLG_Eury.lmer.ranef) # looks pretty good...

genotype_random_model(x = tree_level_interaxn_all_plants, y = with(tree_level_interaxn_all_plants, sqrt(rsLG_Eury/shootEst.no18)))
exactRLRT(rsLG_Eury.lmer)

with(tree_level_interaxn_all_plants, cor.test(sqrt(rsLG_Eury/shootEst.no18), sqrt(rsLG_abund/shootEst.no18)))
with(tree_level_interaxn_all_plants, plot(sqrt(rsLG_Eury/shootEst.no18), sqrt(rsLG_abund/shootEst.no18)))


### Plots of phenotypic correlations in dominant species interactions

dominants.df <- with(tree_level_interaxn_all_plants[-121, ], data.frame(vLG_Platy_density = sqrt(vLG_Platy/shootEst.no18), vLG_Mesopol_density = sqrt(vLG_Mesopol/shootEst.no18), vLG_Tory_density = sqrt(vLG_Tory/shootEst.no18), rG_Tory_density = sqrt(rG_Tory/shootEst.no18))) # removed this data point because it was an outlier for a couple of interactions

scatterplotMatrix(dominants.df)


###### Correlations in 

vLG.ptoid.cors.df <- with(tree_level_interaxn_all_plants[-c(17,121), ], data.frame(vLG_Platy_density = sqrt(vLG_Platy/shootEst.no18), vLG_Mesopol_density = sqrt(vLG_Mesopol/shootEst.no18), vLG_Tory_density = sqrt(vLG_Tory/shootEst.no18), vLG_Mymarid_density = sqrt(vLG_Mymarid/shootEst.no18), vLG_Eulo_density = sqrt(vLG_Eulo/shootEst.no18), vLG_density = sqrt(vLG_abund/shootEst.no18)))

corr.test(vLG.ptoid.cors.df)
scatterplotMatrix(vLG.ptoid.cors.df)

rG.ptoid.cors.df <- with(tree_level_interaxn_all_plants[-c(121), ], data.frame(rG_Platy_density = sqrt(rG_Platy/shootEst.no18), rG_Mesopol_density = sqrt(rG_Mesopol/shootEst.no18), rG_Tory_density = sqrt(rG_Tory/shootEst.no18), rG_Eulo_density = sqrt(rG_Eulo/shootEst.no18), rG_density = sqrt(rG_abund/shootEst.no18)))

corr.test(rG.ptoid.cors.df)
scatterplotMatrix(rG.ptoid.cors.df)


##### Mixed effect models of the determinants of parasitoid attack rates

vLG.gall.mech.df <- filter(gall_mech.df, vLG_total > 0 & vLG_density < 0.10) # this restricts the data frame to only those where vLG galls were present and also removes the outlying data points for vLG_density from the vLG_density ~ Genotype analysis.

resid.D <- residuals(lm(D ~ vLG_density, vLG.gall.mech.df)) # need to figure out how to input

vLG.gall.mech.df <- mutate(vLG.gall.mech.df,
                           s.g.height = scale(g.height, center = TRUE, scale = FALSE),
                           s.vLG_density = scale(vLG_density, center = TRUE, scale = FALSE),
                           sc.g.height = scale(g.height, center = TRUE, scale = TRUE),
                           sc.vLG_density = scale(vLG_density, center = TRUE, scale = TRUE),
                           s.D = scale(D, center = TRUE, scale = FALSE))

### vLG_vLG.pupa
vLG_surv.gamm <- gamm4(cbind(vLG_vLG.pupa, vLG_total-vLG_vLG.pupa) ~ s(g.height, by = vLG_density), 
                        random = ~(1 | Genotype/plant.position),
                        data = filter(vLG.gall.mech.df, g.height < 15), 
                        family = binomial(link = "logit"))
summary(vLG_surv.gamm$gam)
summary(vLG_surv.gamm$mer)
AIC(vLG_surv.gamm$mer)

vis.gam(vLG_surv.gamm$gam, type = "response", ticktype = "detailed", theta = -45, plot.type = "persp", se = 0)
vis.gam(vLG_surv.gamm$gam, type = "response", plot.type = "contour", se = 0)

vLG_surv.glmer <- glmer(cbind(vLG_vLG.pupa, vLG_total-vLG_vLG.pupa) ~ g.height*vLG_density + (1 | Genotype/plant.position),
                       data = vLG.gall.mech.df, 
                       family = binomial(link = "logit"))
summary(vLG_surv.glmer)
visreg(vLG_surv.glmer, "g.height", by = "vLG_density", scale = "response")


AIC(vLG_surv.glmer)

### vLG_Platy

# basic gam model. The "by" model appears to be a much better fit than treating vLG_density as a continuous variable (model comparison no longer shown, judged by difference in AIC and likelihood ratio test)
vLG_Platy.gam.by <- gam(cbind(vLG_Platy, vLG_total-vLG_Platy) ~ s(g.height, by = vLG_density), # D, Total_Area, and Height are all significant. Density and log(Trichomes+1) are not.
                     data = vLG.gall.mech.df, 
                     family = binomial(link = "logit"),
                     method = "GCV.Cp")
summary(vLG_Platy.gam.by)
visreg(vLG_Platy.gam.by, "g.height", by = "vLG_density", scale = "response")
visreg2d(vLG_Platy.gam.by, xvar = "g.height", yvar = "vLG_density", scale = "response")
visreg2d(vLG_Platy.gam.by, xvar = "g.height", yvar = "vLG_density", scale = "response", plot.type = "persp")

# mixed effects gam model. Can't use visreg to visualize
vLG_Platy.gamm <- gamm4(cbind(vLG_Platy, vLG_total-vLG_Platy) ~ s(g.height, by = vLG_density), 
                        random = ~(1 | plant.position),
                        data = vLG.gall.mech.df, 
                        family = binomial(link = "logit"))
summary(vLG_Platy.gamm$gam)
summary(vLG_Platy.gamm$mer)
AIC(vLG_Platy.gamm$mer)

vis.gam(vLG_Platy.gamm$gam, theta = 45, type = "response", color = "terrain", ticktype = "detailed")
vis.gam(vLG_Platy.gamm$gam, type = "response", plot.type = "contour", color = "topo")


#consider scaling "D". Also, D may be correlated with vLG_density so I should take this into account...glmer binomial looks just as good as gamm mixed effects model.
vLG.gall.mech.df.sub <- vLG.gall.mech.df[complete.cases(select(vLG.gall.mech.df, vLG_Platy, vLG_total, s.g.height, s.vLG_density, sc.g.height, sc.vLG_density, Genotype, plant.position)), ] # focus on complete cases
vLG.gall.mech.df.sub$Genotype <- factor(as.character(vLG.gall.mech.df.sub$Genotype))
vLG.gall.mech.df.sub$plant.position <- factor(as.character(vLG.gall.mech.df.sub$plant.position))

vLG_Platy.glmer.visuals <- glmer(cbind(vLG_Platy, vLG_total-vLG_Platy) ~ g.height + vLG_density + g.height:vLG_density + (1 | plant.position), # maybe include D
                         data = vLG.gall.mech.df.sub, 
                         family = binomial(link = "logit"),
                         glmerControl(optimizer="bobyqa"))
summary(vLG_Platy.glmer.visuals) # same coefficient estimate for interaction as model below

vLG_Platy.glmer <- glmer(cbind(vLG_Platy, vLG_total-vLG_Platy) ~ s.g.height + s.vLG_density + s.g.height:s.vLG_density + (1 | Genotype/plant.position), # maybe include D, maybe exclude Genotype as random effect.
                        data = vLG.gall.mech.df.sub, 
                        family = binomial(link = "logit"),
                        glmerControl(optimizer="bobyqa"))
summary(vLG_Platy.glmer)
plot(vLG_Platy.glmer)
AIC(vLG_Platy.glmer)

vLG_Platy.glmer.red <- update(vLG_Platy.glmer, .~. - s.g.height:s.vLG_density)
summary(vLG_Platy.glmer.red)

anova(vLG_Platy.glmer.red, vLG_Platy.glmer)
#PBmodcomp(largeModel = vLG_Platy.glmer, smallModel = vLG_Platy.glmer.red, nsim = 100)

# for some reason, the visreg function is having a hard time plotting the results with a nested random effect model...
visreg(vLG_Platy.glmer.visuals, "g.height", by = "vLG_density", scale = "response")
visreg(vLG_Platy.glmer.visuals, "vLG_density", by = "g.height", scale = "response")
visreg2d(vLG_Platy.glmer.visuals, "g.height", "vLG_density", scale = "response", xlab = "Leaf gall size (mm)", ylab = "Leaf gall density", main = "Probability of Parasitism", cex.lab = 1.5, cex.main = 2)
visreg2d(vLG_Platy.glmer.visuals, "vLG_density", "g.height", scale = "response", xlab = "Leaf gall density", ylab = "Leaf gall size (mm)", main = "Probability of Parasitism", cex.lab = 1.5, cex.main = 2)
visreg2d(vLG_Platy.glmer.visuals, xvar = "g.height", yvar = "vLG_density", plot.type = "persp", scale = "response")

vLG_Platy.pred <- predict(vLG_Platy.glmer.visuals, type = "response")

vLG_Platy.plot.df <- data.frame(g.height = vLG.gall.mech.df.sub$g.height, vLG_density = vLG.gall.mech.df.sub$vLG_density, vLG_Platy.pred = vLG_Platy.pred, vLG_Mesopol.orig = with(vLG.gall.mech.df.sub, (vLG_Platy/vLG_total)))

ggplot(filter(vLG_Platy.plot.df, g.height < 8), aes(x = vLG_density, y = vLG_Platy.pred)) + geom_smooth(aes(x = vLG_density, y = vLG_Platy.pred))

### vLG_Mesopol. Very interesting...Both models are showing lower vLG_Mesopol attack at higher vLG_densities. Perhaps this is because it is being outcompted by vLG_Platy?!?
vLG_Mesopol.gamm4 <- gamm4(cbind(vLG_Mesopol + vLG_Ptero.2, vLG_total-vLG_Mesopol - vLG_Ptero.2) ~ s(g.height, vLG_density), 
                        random = ~(1 | plant.position),# Genotype
                       data = vLG.gall.mech.df.sub, 
                       family = binomial(link = "logit"))
summary(vLG_Mesopol.gamm4$gam)
summary(vLG_Mesopol.gamm4$mer)

vLG_Mesopol.pred <- predict(vLG_Mesopol.gamm4$gam, type = "response")

vLG_Mesopol.plot.df <- data.frame(g.height = vLG.gall.mech.df.sub$g.height, vLG_density = vLG.gall.mech.df.sub$vLG_density, vLG_Mesopol.pred = vLG_Mesopol.pred, vLG_Mesopol.orig = with(vLG.gall.mech.df.sub, (vLG_Mesopol + vLG_Ptero.2)/vLG_total))

ggplot(filter(vLG_Mesopol.plot.df, g.height < 8), aes(x = vLG_density, y = vLG_Mesopol.pred)) + stat_smooth(method = "glm", family = "binomial", se = FALSE) #geom_smooth(aes(x = vLG_density, y = vLG_Mesopol.pred))

vLG_PlatyMesopol.df <- data.frame(g.height = vLG.gall.mech.df.sub$g.height, vLG_density = vLG.gall.mech.df.sub$vLG_density, vLG_Mesopol.pred = vLG_Mesopol.pred, vLG_total = vLG.gall.mech.df.sub$vLG_total, vLG_Mesopol = with(vLG.gall.mech.df.sub, vLG_Mesopol + vLG_Ptero.2), vLG_Platy = vLG.gall.mech.df.sub$vLG_Platy, vLG_Platy.pred = vLG_Platy.pred)

(vLG_PlatyMesopol.df.melt <- melt(vLG_PlatyMesopol.df, id.vars = c("g.height","vLG_density")))

ggplot(vLG_PlatyMesopol.df, aes(x = g.height, y = vLG_Mesopol/vLG_total, weight = vLG_total, size = vLG_total)) + 
  stat_smooth(method = gam, formula = y ~ s(x), se = FALSE, size = 5) + 
  geom_point(shape = 21) +
  theme_classic() +
  theme(legend.position = "none") +
  ylab(expression(paste("Probability of ",italic("Mesopolobus")," parasitism"))) +
  xlab(expression(paste(italic("Iteomyia")," Gall Size (mm)")))

ggplot(vLG_PlatyMesopol.df, aes(x = g.height, y = vLG_Platy/vLG_total, weight = vLG_total, size = vLG_total)) + 
  stat_smooth(method = glm, formula = y ~ x, family = binomial, se = FALSE, size = 5) + 
  geom_point(shape = 21) +
  theme_classic() +
  theme(legend.position = "none") +
  ylab(expression(paste("Probability of ",italic("Platygaster")," parasitism"))) +
  xlab(expression(paste(italic("Iteomyia")," Gall Size (mm)")))

ggplot(filter(vLG_PlatyMesopol.df, g.height < 8), aes(x = vLG_density*100, y = vLG_Mesopol/vLG_total, weight = vLG_total, size = vLG_total)) + 
  stat_smooth(method = glm, family = binomial, se = FALSE, size = 5, color = "black") + 
  geom_point(shape = 21, color = "grey") +
  theme_classic() +
  theme(legend.position = "none") +
  ylab(expression(paste("Probability of ",italic("Mesopolobus")," parasitism"))) +
  xlab(expression(paste("No. ",italic("Iteomyia")," galls per 100 shoots")))

ggplot(filter(vLG_PlatyMesopol.df, g.height < 8), aes(x = vLG_density*100, y = vLG_Platy/vLG_total, weight = vLG_total, size = vLG_total)) + 
  stat_smooth(method = glm, family = binomial, se = FALSE, size = 5, color = "black") + 
  geom_point(shape = 21, color = "grey") +
  theme_classic() +
  theme(legend.position = "none") +
  ylab(expression(paste("Probability of ",italic("Platygaster")," parasitism"))) +
  xlab(expression(paste("No. ",italic("Iteomyia")," galls per 100 shoots")))



vis.gam(vLG_Mesopol.gamm4$gam, theta = 45, type = "response")
vis.gam(vLG_Mesopol.gamm4$gam, view = c("vLG_density", "g.height"), type = "response", plot.type = "contour")
plot(vLG_Mesopol.gamm4$mer)

vLG_Mesopol.gam <- gam(cbind(vLG_Mesopol + vLG_Ptero.2, vLG_total-vLG_Mesopol - vLG_Ptero.2) ~ s(g.height, vLG_density), 
                       data = filter(gall_mech.df, vLG_total > 0), 
                         family = binomial(link = "logit"),
                         method = "GCV.Cp")
summary(vLG_Mesopol.gam)
visreg(vLG_Mesopol.gam,"vLG_density", by = "g.height", scale = "response")
visreg(vLG_Mesopol.gam,"g.height", by = "vLG_density", scale = "response")
visreg2d(vLG_Mesopol.gam, xvar = "g.height", yvar = "vLG_density", scale = "response", plot.type = "persp")
vis.gam(vLG_Mesopol.gam, theta = 45, type = "response")
vis.gam(vLG_Mesopol.gam, type = "response", plot.type = "contour")


### vLG_Tory. 

vLG_Tory.gamm4 <- gamm4(cbind(vLG_Tory.mal + vLG_Tory.fem, vLG_total - vLG_Tory.mal- vLG_Tory.fem) ~ s(g.height) + s(vLG_density), #  didn't pick up anything new when I tried to tease apart the sexes.
                    random = ~(1 | Genotype/plant.position),
                    data = vLG.gall.mech.df.sub, 
                    family = binomial(link = "logit"))
summary(vLG_Tory.gamm4$gam) # no interactive effects
vis.gam(vLG_Tory.gamm4$gam, theta = 45, type = "response", plot.type = "persp", ticktype = "detailed", zlim = c(0,1))
vis.gam(vLG_Tory.gamm4$gam, type = "response", plot.type = "contour", zlim = c(0,1))

# regular gam model did show an interactions between vLG density and gall height, but it is pretty weak. Also, this interaction doesn't hold up in the mixed model.
vLG_Tory.gam <- gam(cbind(vLG_Tory.fem + vLG_Tory.mal, vLG_total - vLG_Tory.fem - vLG_Tory.mal) ~ s(g.height, vLG_density), 
                       data = vLG.gall.mech.df.sub, 
                       family = binomial(link = "logit"))
summary(vLG_Tory.gam)
visreg(vLG_Tory.gam, "g.height", by = "vLG_density", scale = "response")

vLG_Tory.glmer <- glmer(cbind(vLG_Tory.fem + vLG_Tory.mal, vLG_total - vLG_Tory.fem - vLG_Tory.mal) ~ s.g.height + s.vLG_density + (1|Genotype/plant.position), 
                    data = vLG.gall.mech.df.sub, 
                    family = binomial(link = "logit"),
                    glmerControl(optimizer="bobyqa"))
summary(vLG_Tory.glmer)
plot(vLG_Tory.glmer)

vLG_Tory.glmer.red <- update(vLG_Tory.glmer, .~. - s.vLG_density)
vLG_Tory.glmer.red2 <- update(vLG_Tory.glmer, .~. - s.g.height)
vLG_Tory.glmer.red3 <- update(vLG_Tory.glmer, .~. - s.g.height - s.vLG_density)
summary(vLG_Tory.glmer.red)

anova(vLG_Tory.glmer.red, vLG_Tory.glmer)
anova(vLG_Tory.glmer.red2, vLG_Tory.glmer) # suggests only vLG density is an important predictor
anova(vLG_Tory.glmer.red2, vLG_Tory.glmer.red3) # again suggests only vLG density is an important predictor

summary(update(vLG_Tory.glmer.red2, .~. + rG_density)) # no effect of rG_density or an interaction with vLG_density

vLG_Tory.glmer.vis <- update(vLG_Tory.glmer, .~. - s.vLG_density - s.g.height + vLG_density)
summary(vLG_Tory.glmer.vis)
visreg(vLG_Tory.glmer.vis, "vLG_density", scale = "response")

vLG_Tory.pred <- predict(vLG_Tory.glmer.red2, type = "response")

vLG_Tory.glm <- glm(cbind(vLG_Tory, vLG_abund - vLG_Tory) ~ I(vLG_abund/shootEst.no18),
                    data = tree_level_interaxn_all_plants[-c(17,121), ],
                    family = binomial(link = "logit"))
vLG_Tory.glm.null <- glm(cbind(vLG_Tory, vLG_abund - vLG_Tory) ~ 1,
                    data = tree_level_interaxn_all_plants[-c(17,121), ],
                    family = binomial(link = "logit"))
summary(vLG_Tory.glm)
anova(vLG_Tory.glm.null, vLG_Tory.glm, test = "Chisq")
plot(vLG_Tory.glm)
visreg(vLG_Tory.glm, scale = "response") # need to fix with shootEst

vLG_PlatyMesopol.df <- data.frame(g.height = vLG.gall.mech.df.sub$g.height, vLG_density = vLG.gall.mech.df.sub$vLG_density, vLG_Mesopol.pred = vLG_Mesopol.pred, vLG_total = vLG.gall.mech.df.sub$vLG_total, vLG_Mesopol = with(vLG.gall.mech.df.sub, vLG_Mesopol + vLG_Ptero.2), vLG_Platy = vLG.gall.mech.df.sub$vLG_Platy, vLG_Platy.pred = vLG_Platy.pred, vLG_Tory = with(vLG.gall.mech.df.sub, vLG_Tory.fem + vLG_Tory.mal), vLG_Tory.pred = vLG_Tory.pred)

#(vLG_PlatyMesopol.df.melt <- melt(vLG_PlatyMesopol.df, id.vars = c("g.height","vLG_density")))

meso.size.p <- ggplot(vLG_PlatyMesopol.df, aes(x = g.height, y = vLG_Mesopol/vLG_total, weight = vLG_total, size = vLG_total)) + 
  stat_smooth(aes(y = vLG_Mesopol.pred), method = gam, formula = y ~ s(x), se = FALSE, size = 3, color = "black") + 
  geom_point(shape = 21, color = "grey") +
  theme_classic() +
  theme(legend.position = "none") +
  ylab(expression(paste(italic("Mesopolobus")," - % parasitized"))) +
  xlab(expression(paste(italic("Iteomyia")," Gall Size (mm)")))

platy.size.p <- ggplot(vLG_PlatyMesopol.df, aes(x = g.height, y = vLG_Platy/vLG_total, size = vLG_total)) + 
  stat_smooth(aes(y = vLG_Platy.pred), method = glm, formula = y ~ x, family = binomial, se = FALSE, size = 3, color = "black") + 
  geom_point(shape = 21, color = "grey") +
  theme_classic() +
  theme(legend.position = "none") +
  ylab(expression(paste(italic("Platygaster")," - % parasitized"))) +
  xlab(expression(paste(italic("Iteomyia")," Gall Size (mm)")))

meso.dens.p <- ggplot(filter(vLG_PlatyMesopol.df, g.height < 8), aes(x = vLG_density*100, y = vLG_Mesopol/vLG_total, size = vLG_total)) + 
  stat_smooth(aes(y = vLG_Mesopol.pred), method = glm, family = binomial, se = FALSE, size = 3, color = "black") + 
  geom_point(shape = 21, color = "grey") +
  theme_classic() +
  theme(legend.position = "none") +
  ylab(expression(paste(italic("Mesopolobus")," - % small galls parasitized (< 8 mm)"))) +
  xlab(expression(paste("No. ",italic("Iteomyia")," galls per 100 shoots")))

platy.dens.p <- ggplot(filter(vLG_PlatyMesopol.df, g.height < 8), aes(x = vLG_density*100, y = vLG_Platy/vLG_total, weight = vLG_total, size = vLG_total)) + 
  stat_smooth(aes(y = vLG_Platy.pred), method = glm, family = binomial, se = FALSE, size = 3, color = "black") + 
  geom_point(shape = 21, color = "grey") +
  theme_classic() +
  theme(legend.position = "none") +
  ylab(expression(paste(italic("Platygaster")," - % small galls parasitized (< 8 mm)"))) +
  xlab(expression(paste("No. ",italic("Iteomyia")," galls per 100 shoots")))

tory.dens.p <- ggplot(vLG_PlatyMesopol.df, aes(x = vLG_density*100, y = vLG_Tory/vLG_total, size = vLG_total)) + 
  stat_smooth(aes(y = vLG_Tory.pred), method = glm, family = binomial, se = FALSE, size = 3, color = "black") + 
  geom_point(shape = 21, color = "grey") +
  theme_classic() +
  theme(legend.position = "none") +
  ylab(expression(paste(italic("Torymus")," - % parasitized"))) +
  xlab(expression(paste("No. ",italic("Iteomyia")," galls per 100 shoots")))

grid.arrange(platy.size.p, platy.dens.p, tory.dens.p, meso.size.p, meso.dens.p, ncol = 3)

with(filter(vLG_PlatyMesopol.df, g.height < 4.3), mean(vLG_Platy.pred))/with(filter(vLG_PlatyMesopol.df, g.height > 14), mean(vLG_Platy.pred)) # Platygaster attack decreased 18-fold from small to large galls.

with(filter(vLG_PlatyMesopol.df, g.height < 8 & vLG_density > 0.08), mean(vLG_Platy.pred))/with(filter(vLG_PlatyMesopol.df, g.height < 8 & vLG_density < 0.01), mean(vLG_Platy.pred))

with(filter(vLG_PlatyMesopol.df, g.height < 8 & vLG_density < 0.01), mean(vLG_Mesopol.pred))/with(filter(vLG_PlatyMesopol.df, g.height < 8 & vLG_density > 0.08), mean(vLG_Mesopol.pred))

with(filter(vLG_PlatyMesopol.df, g.height < 8 & vLG_density < 0.01), mean(vLG_Tory.pred))/with(filter(vLG_PlatyMesopol.df, g.height < 8 & vLG_density > 0.06), mean(vLG_Tory.pred))

#### rG_Tory
rG.gall.mech.df <- filter(gall_mech.df, rG_total > 0 & rG_density < 0.5) # this restricts the data frame to only those where rG galls were present and also removes the outlying data points for rG_density from the rG_density ~ Genotype analysis.

rG_Tory.glm <- glm(cbind(rG_Tory.fem + rG_Tory.mal, rG_total - rG_Tory.fem - rG_Tory.mal) ~ rG_density, 
                    data = rG.gall.mech.df, 
                    family = binomial(link = "logit"))
summary(rG_Tory.glm) # no relationship
visreg(rG_Tory.glm, scale = "response")

rG_Tory.glm2 <- glm(cbind(rG_Tory, rG_abund - rG_Tory) ~ I(rG_abund/shootEst.no18), 
                   data = tree_level_interaxn_all_plants[-121, ], 
                   family = binomial(link = "logit"))
summary(rG_Tory.glm2) # no relationship

##### Plots for dominant interactions
par(mfrow = c(2,2),
    mar = c(1.5,1.5,1.5,1.5))
vis.gam(vLG_Platy.gamm$gam, theta = -45, type = "response", plot.type = "persp", ticktype = "detailed", zlim = c(0,1), xlab = "Gall Size (mm)", ylab = "Iteomyia density", zlab = "Probability of Parasitism", main = "Platygaster")
vis.gam(vLG_Mesopol.gamm4$gam, theta = -45, type = "response", plot.type = "persp", ticktype = "detailed", zlim = c(0,0.17), xlab = "Gall Size (mm)", ylab = "Iteomyia density", zlab = "Probability of Parasitism", main = "Mesopolobus")
vis.gam(vLG_Tory.gamm4$gam, theta = -45, type = "response", plot.type = "persp", zlim = c(0,0.11), ticktype = "detailed", xlab = "Gall Size (mm)", ylab = "Iteomyia density", zlab = "Probability of Parasitism", main = "Torymus")

################ EVERYTHING BELOW THIS IS OLD
tree_level_interaxn_all_plants <- read.csv("~/Documents/Genotype_Networks/data/tree_level_interaxn_all_plants.csv")

vLG <- tree_level_interaxn_all_plants$vLG_abund
log.shoot.count <- log(tree_level_interaxn_all_plants$shootEst.no18)
Genotype <- tree_level_interaxn_all_plants$Genotype

test <- rpt.poisGLMM.multi(y = vLG, groups = Genotype) # generates and error but may be okay.

plot(vLG/tree_level_interaxn_all_plants$shootEst.no18 ~ Genotype)
which(vLG/tree_level_interaxn_all_plants$shootEst.no18 > 0.1)
vLG.mod.ests <- pqlglmm.pois.model(vLG[-c(17,121)], Genotype[-c(17,121)], "log", returnR=FALSE)

t <- log.shoot.count[-c(17,121)]
glmmPQL(vLG[-c(17,121)] ~ 1 + offset(t), random = ~1|Genotype[-c(17,121)], family = quasipoisson("log"))

test <- rpt.poisGLMM.multi(y = vLG, groups = Genotype, link = "log", nboot = 100, npermut = 100)

nboot <- 100
k <- levels(Genotype)
N <- length(vLG)

vLG.R.boot   <- replicate(nboot, bootstr(vLG, log.shoot.count, Genotype, k, N, vLG.mod.ests$beta0, vLG.mod.ests$var.a, vLG.mod.ests$omega, "log"), simplify=TRUE)
vLG.R.boot   <- list(R.link = as.numeric(unlist(vLG.R.boot["R.link",])), R.org = as.numeric(unlist(vLG.R.boot["R.org",]))) 	

R <-vLG.mod.ests <- pqlglmm.pois.model.offset(vLG[-c(17,121)], offset.trans = log.shoot.count[-c(17,21)], Genotype[-c(17,121)], "log")

npermut = 100
R.permut <- replicate(npermut-1, permut(y = vLG, offset.trans = log.shoot.count, groups = Genotype, N = N, link = "log"), simplify=TRUE)
R.permut <- list(R.link = c(R$R.link, unlist(R.permut["R.link",])), R.org = c(R$R.org, unlist(R.permut["R.org",])))
P.link   <- sum(R.permut$R.link >= R$R.link) / npermut
P.org    <- sum(R.permut$R.org >= R$R.org) / npermut

m1 <- glm(vLG ~ 1)
m2 <- glm(vLG ~ offset(log.shoot.count))

bootstr(vLG, log.shoot.count, Genotype, k, N, vLG.mod.ests$beta0, vLG.mod.ests$var.a, vLG.mod.ests$omega, "log")
#
mu/(vLG.mod.ests$omega-1)
# work on bootstrap
groupMeans <- rnorm(26, 0, sqrt(vLG.mod.ests$var.a))
mu <- exp(vLG.mod.ests$beta0 + groupMeans[Genotype])
y.boot <- rnbinom(N, size=(mu/(vLG.mod.ests$omega-1)), mu=mu)
pqlglmm.pois.model.offset(y.boot, offset.trans = log.shoot.count, groups = Genotype, link = "log")


####### TESTING OUT MIXED EFFECT MODEL ON VLG SURVIVAL
tree_level_interaxn_all_plants <- read.csv("~/Documents/Genotype_Networks/data/tree_level_interaxn_all_plants.csv")
vLGfocus.tree <- tree_level_interaxn_all_plants %>%
  mutate(vLG_density_tree = vLG_abund/shootEst.no18*200) %>%
  select(Genotype, plant.position, vLG_density_tree)

vLG_model.df <- left_join(gall_parasitoid_size_niche_data, vLGfocus.tree, by = "plant.position")
vLG_model.df <- filter(vLG_model.df, vLG_density_tree > 0 & vLG_total > 0)
vLG_model.df <- mutate(vLG_model.df, c.vLG_density = scale(vLG_density_tree, center = TRUE, scale = FALSE),
                       c.g.height = scale(g.height, center = TRUE, scale = FALSE),
                       vLG_density_cut2 = cut(vLG_density_tree, c(0, 4.054, 33.7), labels = c("low density","high density")))

# consider removing the plants that I removed from the density estimations for the different genotypes...
library(lme4) # FOLLOW 10 STEP PROCESS BEFORE DECIDING ON THIS FINAL MODEL!!!
library(lattice)
dotchart(vLG_model.df$vLG_vLG.pupa, groups = factor(vLG_model.df$Genotype))
plot(vLG_vLG.pupa ~ vLG_density_tree, vLG_model.df)
pairs(cbind(vLG_model.df$vLG_density_tree, vLG_model.df$g.height)) # no clear correlation among predictor variables
table(vLG_model.df$vLG_vLG.pupa)
table(vLG_model.df$vLG_Platy)
dotchart(table(vLG_model.df$plant.position))
plot(table(vLG_model.df$plant.position) ~ vLG_model.df$vLG_density_tree)
pp.3plus.vLG <- names(which(table(vLG_model.df$plant.position) > 2))

vLG_glmer <- glmer(cbind(vLG_vLG.pupa, vLG_total - vLG_vLG.pupa) ~ c.vLG_density*c.g.height + (1 | Genotype/plant.position), data = filter(vLG_model.df, plant.position %in% pp.3plus.vLG), family = "binomial")
summary(vLG_glmer)
plot(vLG_glmer)
qqnorm(vLG_glmer)
vLG_glmer_resid <- residuals(vLG_glmer, type = "pearson")

# don't have to check the random effects, because I have already accounted for the correlation structure! I DO NEED TO CHECK AND MAKE SURE THE EXPLANATORY VARIABLES ARE OKAY.
plot(vLG_glmer_resid ~ Genotype, vLG_model.df[names(vLG_glmer_resid), ]) # note heterogeneity in variance of residuals among genotypes
plot(vLG_glmer_resid ~ factor(plant.position), vLG_model.df[names(vLG_glmer_resid), ]) # note heterogeneity in variance of residuals among genotypes
plot(vLG_glmer_resid ~ c.g.height, vLG_model.df[names(vLG_glmer_resid), ])
plot(vLG_glmer_resid ~ c.vLG_density, vLG_model.df[names(vLG_glmer_resid), ])

vLG_glm <- glm(cbind(vLG_vLG.pupa, vLG_total - vLG_vLG.pupa) ~ vLG_density_tree * g.height, data = vLG_model.df, family = "binomial")
summary(vLG_glm)
plot(vLG_glm)

vLG_glm.resids <- residuals(vLG_glm, type = "pearson")
plot(vLG_glm.resids ~ Genotype, vLG_model.df[names(vLG_glm.resids), ])
plot(vLG_glm.resids ~ factor(plant.position), vLG_model.df[names(vLG_glm.resids), ])

visreg(vLG_glm, "g.height", by = "vLG_density_cut2", scale = "response")

### vLG_Platy glmer. Consider doing non-linear term
vLG_Platy_glmer <- glmer(cbind(vLG_Platy, vLG_total - vLG_Platy) ~ c.vLG_density*c.g.height + (1 | Genotype/plant.position), data = vLG_model.df, family = "binomial")
summary(vLG_Platy_glmer)
plot(vLG_Platy_glmer)

vLG_Platy_glm <- glm(cbind(vLG_Platy, vLG_total - vLG_Platy) ~ vLG_density_tree * g.height, data = vLG_model.df, family = "binomial") 
summary(vLG_Platy_glm)
plot(vLG_Platy_glm)

visreg(vLG_Platy_glm, "g.height", by = "vLG_density_tree", scale = "response")



### vLG_Mesopol glmer

vLG_Mesopol_glmer <- glmer(cbind(vLG_Mesopol + vLG_Ptero.2, vLG_total - vLG_Mesopol - vLG_Ptero.2) ~ c.vLG_density*c.g.height + (1 | Genotype/plant.position), data = vLG_model.df, family = "binomial")
summary(vLG_Mesopol_glmer)
plot(vLG_Mesopol_glmer)

vLG_Mesopol_gamm4 <- gamm4(cbind(vLG_Mesopol + vLG_Ptero.2, vLG_total - vLG_Mesopol - vLG_Ptero.2) ~  s(g.height), 
                           random = ~ (1 | plant.position), 
                           data = vLG_model.df, family = "binomial")
summary(vLG_Mesopol_gamm4$mer)
summary(vLG_Mesopol_gamm4$gam)
plot(vLG_Mesopol_gamm4$mer)

plot(vLG_Mesopol_gamm4$gam, se = 1, seWithMean = TRUE, rug = FALSE, shift = mean(predict(vLG_Mesopol_gamm4$gam)),
     trans = function(x){exp(x)/(1+exp(x))})

vLG_Mesopol_glm <- gam(cbind(vLG_Mesopol, vLG_total - vLG_Mesopol) ~ s(g.height), data = vLG_model.df, family = "binomial", method = "GCV.Cp") 
summary(vLG_Mesopol_glm)
plot(vLG_Mesopol_glm)

visreg(vLG_Mesopol_glm, "g.height", scale = "response")