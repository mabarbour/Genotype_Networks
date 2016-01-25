

library(dplyr)
library(reshape2)
library(visreg)
library(mgcv)
library(ggplot2)

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



#library(arm)
#rpt.remlLMM.adj(formula = g.height ~ 1 + (1 | Genotype/plant.position), 
 #               data = vLG.size.df, 
  #              grname = "Genotype",
   #             nboot = 10, 
    #            npermut = 10)


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

### Use cubic splines to identify gall attack thresholds for different parasitoid species

# rG. No evidence for a size affect rG survival so I did not analyze gall size for this species.
rG_df <- filter(gall_parasitoid_size_niche_data, rG_total > 0)

ggplot(rG_df, aes(y = rG_rG.larv/rG_total, x = g.height)) + geom_point() + stat_smooth(method = "glm", family = binomial, aes(weight = rG_df$rG_total)) # no clear effect.

rG.larv.glm <- glm(cbind(rG_rG.larv, rG_total-rG_rG.larv) ~ g.height, rG_df, family = "binomial")
summary(rG.larv.glm)

rG_larv.thres <- gam(cbind(rG_rG.larv, rG_total - rG_rG.larv) ~ s(g.height), 
                      rG_df,
                      family = binomial(link = "logit"),
                      method = "GCV.Cp")
summary(rG_larv.thres)
plot(rG_larv.thres, se = 1, seWithMean = TRUE, rug = FALSE, shift = mean(predict(rG_larv.thres)),
     trans = function(x){exp(x)/(1+exp(x))})

rG_Tory.thres <- gam(cbind(rG_Tory.fem + rG_Tory.mal, rG_total - rG_Tory.fem + rG_Tory.mal) ~ s(g.height), 
                     rG_df,
                     family = binomial(link = "logit"),
                     method = "GCV.Cp")
summary(rG_Tory.thres)
plot(rG_Tory.thres, se = 1, seWithMean = TRUE, rug = FALSE, shift = mean(predict(rG_Tory.thres)),
     trans = function(x){exp(x)/(1+exp(x))})

### vLG. Evidence for gall size affecting both vLG survival, as well as probability of Platy and Mesopolobus attack 

vLG_df <- filter(gall_parasitoid_size_niche_data, vLG_total > 0)

hist(vLG_df$g.height)
median(vLG_df$g.height, na.rm = T) # 8.36
mean(vLG_df$g.height, na.rm = T) # 8.62

# vLG pupa survival
ggplot(data = vLG_df, aes(x = g.height, y = vLG_vLG.pupa/vLG_total)) + 
  geom_point(pch = 1, size = 5) + 
  stat_smooth(method = "glm", family = binomial, aes(weight = vLG_df$vLG_total), size = 2, se = TRUE) + 
  theme_classic() +
  theme(axis.text = element_text(size = 18),
        axis.title.x = element_text(size = 22, vjust = -0.5),
        axis.title.y = element_text(size = 22, vjust = 1)) +
  ylab("Probability of Iteomyia survival") + 
  xlab("Gall height (mm)")
  
vLG.pupa_glm <- glm(cbind(vLG_vLG.pupa, vLG_total - vLG_vLG.pupa) ~ g.height, vLG_df, family = binomial(link = "logit"))
summary(vLG.pupa_glm)
plot(vLG.pupa_glm)
visreg(vLG.pupa_glm, scale = "response")

vLG_pupa.thres <- gam(cbind(vLG_vLG.pupa, vLG_total - vLG_vLG.pupa) ~ s(g.height), 
                      filter(gall_parasitoid_size_niche_data, vLG_total > 0),
                      family = binomial(link = "logit"),
                      method = "GCV.Cp")
summary(vLG_pupa.thres)
plot(vLG_pupa.thres, se = 1, seWithMean = TRUE, rug = FALSE, shift = mean(predict(vLG_pupa.thres)),
     trans = function(x){exp(x)/(1+exp(x))})



vLG_Platy.thres <- gam(cbind(vLG_Platy, vLG_total-vLG_Platy) ~ s(g.height), 
                       filter(vLG_model.df, vLG_density_cut2 == "high density"), 
                       family = binomial(link = "logit"),
                       method = "GCV.Cp")
summary(vLG_Platy.thres)
plot(vLG_Platy.thres, se = 1.96, seWithMean = TRUE, rug = FALSE, shift = mean(predict(vLG_Platy.thres)),
     trans = function(x){exp(x)/(1+exp(x))})

vLG_Mesopol.thres <- gam(cbind(vLG_Mesopol + vLG_Ptero.2, vLG_total-vLG_Mesopol - vLG_Ptero.2) ~ s(g.height), 
                         vLG_model.df, 
                         family = binomial(link = "logit"),
                         method = "GCV.Cp")
summary(vLG_Mesopol.thres)
plot(vLG_Mesopol.thres, se = 1, seWithMean = TRUE, rug = FALSE, shift = mean(predict(vLG_Mesopol.thres)),
     trans = function(x){exp(x)/(1+exp(x))})
gam.check(vLG_Mesopol.thres)



### Gall sizes at which different gall occupants are found

gall_parasitoid_size_niche_data %>%
  filter(vLG_vLG.pupa > 0) %>%
  summarise(vLG.pupa_niche = mean(g.height, na.rm = T), vLG.pupa_sd = sd(g.height, na.rm = TRUE), N = n())

gall_parasitoid_size_niche_data %>%
  filter(vLG_Platy > 0) %>%
  summarise(vLG_Platy_niche = mean(g.height, na.rm = T), vLG_Platy_sd = sd(g.height, na.rm = TRUE), N = n())

gall_parasitoid_size_niche_data %>%
  filter(vLG_Mesopol + vLG_Ptero.2 > 0) %>%
  summarise(vLG_Mesopol_niche = mean(g.height, na.rm = TRUE), vLG_Mesopol_sd = sd(g.height, na.rm = TRUE), N = n())

gall_parasitoid_size_niche_data %>%
  filter(vLG_Tory.fem + vLG_Tory.mal > 0) %>%
  summarise(vLG_Tory_niche = mean(g.height, na.rm = TRUE), vLG_Tory_sd = sd(g.height, na.rm = TRUE), N = n())

gall_parasitoid_size_niche_data %>%
  filter(vLG_Eulo.mal + vLG_Eulo.fem > 0) %>%
  summarise(vLG_Eulo_niche = mean(g.height, na.rm = T), vLG_Eulo_sd = sd(g.height, na.rm = T), N = n())


### Deciding which measure of gall size to use. I have decided to use gall height since it is a direct measurement from my data and appears to do a slightly better job at predicting the probability of attack (AIC and p-value for vLG survival)
plot(log(gall.volume.height.width.cm) ~ log(gall.volume.cm), gall_size_filter_add) # virtually the same, so I decided to use gall.volume.height.width.cm since I have more data for it.
plot(gall.volume.height.width.cm ~ g.height, filter(gall_size_filter_add, gall.sp == "vLG")) # I wonder whether there is anything important to this curvilinear relationship.

vLG_size_prop_df <- gall_size_filter_add %>%
  group_by(Genotype, plant.position, gall.sp, gall.id.nest) %>%
  filter(g.height > 0) %>% # removes NA data
  summarise(g.height = mean(g.height), N = n()) %>%
  ungroup() %>%
  group_by(Genotype, plant.position, gall.sp) %>%
  filter(gall.sp == "vLG") %>%
  mutate(g.height.cut = ifelse(g.height > 8.62, "large", "small"))

vLG_mean_g.height.df <- summarise(vLG_size_prop_df, g.height.mean = mean(g.height))
vLG_mean_g.height.df <- vLG_mean_g.height.df %>%
  ungroup() %>%
  mutate(plant.position = as.character(plant.position)) %>%
  select(Genotype, plant.position, g.height.mean)
  

vLG_g.height_df <- as.data.frame.matrix(with(vLG_size_prop_df, table(plant.position, g.height.cut))) 
vLG_g.height_df <- mutate(vLG_g.height_df, plant.position = rownames(vLG_g.height_df))


vLG_size_df <- left_join(vLG_mean_g.height.df, vLG_g.height_df)


### How does mean gall size vary among genotypes?
vLG.size.geno.sum <- vLG_size_df %>%
  group_by(Genotype) %>%
  summarise(vLG.size.mean = mean(g.height.mean), 
            vLG.size.se = sd(g.height.mean)/sqrt(length(g.height.mean))) %>%
  arrange(vLG.size.mean)

vLG.size.geno.sum <- mutate(vLG.size.geno.sum, Geno.ord = factor(Genotype, levels = as.character(Genotype), ordered = TRUE))
levels(vLG.size.geno.sum$Geno.ord)[4] <- "C" # for aesthetics

ggplot(vLG.size.geno.sum, aes(x = Geno.ord, y = vLG.size.mean)) + 
  geom_errorbar(aes(ymin = vLG.size.mean - vLG.size.se, ymax = vLG.size.mean + vLG.size.se), width = 0.2) +
  geom_point(shape = 21, fill = "steelblue", size = 10) + 
  theme + 
  xlab ("Genotype") + 
  scale_y_continuous(breaks = seq(0,14,2), limits = c(0,14)) +
  ylab ("Gall size (mm)")
ggsave("~/Documents/Genotype_Networks/vLG_gall_height_variation.png", width = 11.5, height = 8, units = "in", dpi = 300)


g.height.vLG.lm <- lm(log(g.height.mean) ~ Genotype, vLG_size_df[-c(5,44,48), ]) # remove genotypes with only one data point. Doesn't affect qualitative results.
summary(g.height.vLG.lm)
anova(g.height.vLG.lm)
plot(g.height.vLG.lm)

plot(small/(large + small) ~ Genotype, vLG_size_df)
g.height.prop.vLG <- glm(cbind(small, large) ~ Genotype, vLG_size_df[-c(5,44,48),], family = binomial(link = "logit"))
summary(g.height.prop.vLG)
plot(g.height.prop.vLG)
anova(g.height.prop.vLG, test = "Chisq")

### Genotype level data
vLG_gall_size_geno.df <- vLG_size_df %>%
  group_by(Genotype) %>%
  summarise(g.height.vLG.mean = mean(g.height.mean), large.vLG.count = sum(large), small.vLG.count = sum(small))

write.csv(vLG_gall_size_geno.df, "~/Documents/Genotype_Networks/data/vLG_gall_size_geno.df.csv")

#### everything below this is old!!!
tree_level_interaxn_all_plants <- read.csv('~/Documents/Genotype_Networks/data/tree_level_interaxn_all_plants.csv')
t <- with(gall_size_prop_df, table(plant.position, g.height.cut))
t <- as.data.frame.matrix(t)
t <- mutate(t, plant.position = rownames(t))
#vLG_size_class <- data.frame(with(gall_size_prop_df, table(plant.position, g.height.cut)))
#vLG_size_class <- mutate(vLG_size_class, plant.position = as.character(plant.position))
tree_level_interaxn_all_plants$plant.position <- as.character(tree_level_interaxn_all_plants$plant.position)

vLG_size_class_genotype <- left_join(t, tree_level_interaxn_all_plants)

plot(small/(large+small) ~ Genotype, vLG_size_class_genotype)
vLG_class_glm <- glm(cbind(small, large) ~ Genotype, data = vLG_size_class_genotype[-c(7,63,69),], family = "binomial")
summary(vLG_class_glm)
anova(vLG_class_glm, test = "Chi")
plot(vLG_class_glm)
visreg(vLG_class_glm, scale = "response")

plot(small ~ Genotype, vLG_size_class_genotype)
plot(small ~ large, vLG_size_class_genotype)

vLG_small_glm <- glm(cbind(small,large) ~ offset(log(shootEst.no18)) + Genotype, data = vLG_size_class_genotype[-c(7,63,69),], family = "binomial")
summary(vLG_small_glm)
anova(vLG_small_glm, test= "Chi")
plot(vLG_small_glm)
  

gall_size_data <- gall_size_filter_add %>%
  group_by(Genotype, plant.position, gall.sp, gall.id.nest) %>%
  filter(gall.volume.height.width.cm > 0) %>% # removes NA data
  summarise(gall.size = mean(gall.volume.height.width.cm), N = n()) %>%
  ungroup() %>%
  group_by(Genotype, plant.position, gall.sp) %>%
  summarise(gall.size.mean = mean(gall.size), gall.size.sd = sd(gall.size), N = n()) # note that all of these gall sizes either contained a surviving gall larva or a parasitoid.
  

#### Gall size data set

# vLG size

plot(gall.size.sd ~ gall.size.mean, filter(gall_size_data, gall.sp == "vLG")) # pretty tight correlation
plot(gall.size.mean ~ Genotype, filter(gall_size_data, gall.sp == "vLG")) # a lot of variance for some of the gall species.

vLG_gall.vol.lm <- lm(log(gall.size.mean) ~ Genotype, filter(gall_size_data, gall.sp == "vLG")) # had to remove these replicates because they had a large leverage on the outcome and there was only replicate per genotype. Same results whether or not a "weight" the model with the number of gall sizes sampled
summary(vLG_gall.vol.lm)
anova(vLG_gall.vol.lm) # results hold after removing points with high leverage
#plot(vLG_gall.vol.lm) # variance is "squished" on either side...

vLG_size_visreg <- visreg(vLG_gall.vol.lm, trans = exp, ylab = "vLG volume")
vLG_size_visreg$Genotype$y$fit
vLG_size_summary <- gall_size_data %>%
  filter(gall.sp == "vLG") %>%
  group_by(Genotype) %>%
  summarise(vLG_volume_mean = mean(gall.size.mean), vLG_genotype_size_N = n()) %>%
  mutate(vLG_volume_modelfit = vLG_size_visreg$Genotype$y$fit)


# rG
plot(gall.size.mean ~ Genotype, filter(gall_size_data, gall.sp == "rG"))

rG_gall.vol.lm <- lm(log(gall.size.mean) ~ Genotype, filter(gall_size_data, gall.sp == "rG"))
summary(rG_gall.vol.lm)
#plot(rG_gall.vol.lm)
anova(rG_gall.vol.lm) # analysis still holds when removing results with leverage of 1.
rG_size_visreg <- visreg(rG_gall.vol.lm, "Genotype", trans=exp, ylab = "rG Gall Volume")

rG_size_summary <- gall_size_data %>%
  filter(gall.sp == "rG") %>%
  group_by(Genotype) %>%
  summarise(rG_volume_mean = mean(gall.size.mean), rG_genotype_size_N = n()) %>%
  mutate(rG_volume_modelfit = rG_size_visreg$Genotype$y$fit)
  

# rsLG. Not including rsLG size as a factor since it didn't significantly vary among the genotypes.
plot(gall.size.mean ~ Genotype, filter(gall_size_data, gall.sp == "rsLG")[-c(1,2,22,26,37,42,43), ])

rsLG_gall.vol.lm <- lm(log(gall.size.mean) ~ Genotype, filter(gall_size_data, gall.sp == "rsLG")[-c(1,2,22,26,37,42,43), ])
summary(rsLG_gall.vol.lm)
plot(rsLG_gall.vol.lm)
anova(rsLG_gall.vol.lm) # still not significant even after data points are removed.
visreg(rsLG_gall.vol.lm, "Genotype", trans=exp, ylab = "rsLG Gall Volume")


# aSG. No difference in size, therefore I did not summarize this data.
plot(gall.size.mean ~ Genotype, filter(gall_size_data, gall.sp == "aSG"))

aSG_gall.vol.lm <- lm(log(gall.size.mean) ~ Genotype, filter(gall_size_data, gall.sp == "aSG")[-c(1,2,10,14,19,20,21,26,27,28,29), ])
summary(aSG_gall.vol.lm)
#plot(aSG_gall.vol.lm)
anova(aSG_gall.vol.lm) # still not significant after removing points with high leverage.
visreg(aSG_gall.vol.lm, "Genotype", trans=exp, ylab = "aSG Gall Volume")

# SG
plot(gall.size.mean ~ Genotype, filter(gall_size_data, gall.sp == "SG")) # not a sufficient number of reps per genotype to run the analysis (only one genotype with 2 reps)

# Join together gall size summary data for vLG and rG
anti_join(vLG_size_summary, rG_size_summary)

gall_size_genotype_summary_df <- left_join(vLG_size_summary, rG_size_summary) # vLG has a more encompassing number of genotypes, which is why I did left_join.

write.csv(gall_size_genotype_summary_df, "~/Documents/Genotype_Networks/data/gall_size_genotype_summary_df.csv")
