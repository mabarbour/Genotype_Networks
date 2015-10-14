
## load required libraries
require(visreg)
require(ggplot2)
require(grid)
require(dplyr)
require(piecewiseSEM)
require(semPlot)

## upload results from simulation
#all.measures.500.4 <- read.csv("~/Documents/Genotype_Networks/data/food web complexity simulation 500sims 4reps.csv")
#all.measures.1000.4 <- read.csv("~/Documents/Genotype_Networks/data/food web complexity simulation 1000sims 4reps.csv")
#sample.rare.100 <- read.csv("~/Documents/Genotype_Networks/data/food web complexity 100 sample-based rarefaction.csv")
#all.measures.1000.1.3 <- read.csv("~/Documents/Genotype_Networks/data/food web complexity 1000sims 1-3reps.csv")
#all.measures.2500.4 <- read.csv("~/Documents/Genotype_Networks/data/food web complexity simulation 2500sims 4reps.csv")
#all.measures.mono <- read.csv("~/Documents/Genotype_Networks/data/food web complexity simulation monocultures 1000sims 1-4reps.csv")

all.measures <- read.csv("~/Documents/Genotype_Networks/data/food web complexity simulation 40 reps of 100 sims.csv")
dim(all.measures)[1] # 2221 unique simulations. 2425 simulations originally run
table(all.measures$genotypes.sampled)

alpha.summary <- all.measures %>% 
  filter(genotypes.sampled == 1) %>%
  group_by(df.sim.number) %>%
  summarise(alpha.mean = mean(total_complexity),
            alpha.min = min(total_complexity)) %>%
  mutate(alpha.contrib = alpha.mean - alpha.min)

gamma.summary <- all.measures %>% 
  filter(genotypes.sampled == 25) %>%
  group_by(df.sim.number) %>%
  summarise(gamma = mean(total_complexity))

beta.summary <- left_join(alpha.summary, gamma.summary, by = "df.sim.number") %>%
  mutate(beta.contrib = gamma - alpha.mean,
         total.contrib = gamma - alpha.min,
         beta.prop.contrib = beta.contrib/total.contrib,
         alpha.prop.contrib = alpha.contrib/total.contrib,
         beta.magnitude = gamma/alpha.mean,
         gamma.magnitude = gamma/alpha.min)

final.summary <- beta.summary %>%
  filter(alpha.min > 0) %>%
  summarise_each(funs(mean, sd)) %>%
  select(-df.sim.number_mean, -df.sim.number_sd) %>%
  c()

summary.data.plots <- all.measures %>%
  group_by(df.sim.number, sim.number, genotypes.sampled) %>%
  summarise(total_complexity_mean = mean(total_complexity, na.rm = TRUE)) %>%
  group_by(df.sim.number, genotypes.sampled) %>%
  summarise(total_complexity_mean.mean = mean(total_complexity_mean, na.rm = TRUE)) %>%
  group_by(genotypes.sampled) %>%
  summarise(total_complexity_mean.mean.mean = mean(total_complexity_mean.mean, na.rm = TRUE))
summary.data.plots$total_complexity_mean.mean.mean[25]/summary.data.plots$total_complexity_mean.mean.mean[1]

## piecewise SEM of interactions
ggplot(all.measures, aes(x = genotypes.sampled, y = total_complexity, group = df.sim.number)) +
  geom_point(shape = 1, color = "grey") +
  geom_smooth(se = FALSE) +
  #geom_line(data = sample.rare.100.sum, aes(x = reps.sampled, y = mean_complexity)) +
  theme_classic()

## correlations
car::scatterplotMatrix(select(all.measures.1000.4, genotype.richness, total.gall.abund, total_complexity))

## determine sensitive of genetic variation coefficient to number of simulations
lm.500.4 <- lm(total_complexity ~ log(genotype.richness), all.measures.500.4)
lm.1000.4 <- lm(total_complexity ~ log(genotype.richness), all.measures.1000.4)
lm.2500.4 <- lm(total_complexity ~ log(genotype.richness), all.measures.2500.4)

# coefficients for intercept and log(genotype.richness) appear to stabilize after 1000 simulations for each level of genetic variation.
coef(lm.500.4)
coef(lm.1000.4)
coef(lm.2500.4)

# model with covariates
ggplot(all.measures.1000.4, aes(x = genotype.richness, y = total_complexity)) +
  geom_point() +
  geom_smooth() 

## generalized additive models
library(gam)
all.measures.1000.4.NAone <- all.measures.1000.4
all.measures.1000.4.NAone$total_complexity[is.na(all.measures.1000.4$total_complexity)] <- 1
summary(all.measures.1000.4.NAone$total_complexity)

gam.1000.4.cov <- gam(total_complexity ~ s(genotype.richness) + 
                        s(total.gall.abund),
                      data = all.measures.1000.4.NAone)
summary(gam.1000.4.cov)
#visreg(gam.1000.4.cov)
visreg(gam.1000.4.cov, xvar = "genotype.richness", jitter = TRUE,
       ylab = expression(paste("Food-web complexity (",italic(LD[q]),")")),
       cond = list(total.gall.abund = mean(all.measures.1000.4.NAone$total.gall.abund)),
       xlab = "No. of willow genotypes",
       band = FALSE,
       partial = TRUE,
       points.par = list(cex = 1.1, pch = 21, col = "grey"),
       line.par = list(lwd = 5, col = "steel blue"),
       xaxt = "n")
axis(side = 1, c(1,5,10,15,20,25))


pred.1geno <- data.frame(genotype.richness = 1,
                         total.gall.abund = mean(all.measures.1000.4.NAone$total.gall.abund))
pred.8geno <- data.frame(genotype.richness = 8,
                         total.gall.abund = mean(all.measures.1000.4.NAone$total.gall.abund))
pred.25geno <- data.frame(genotype.richness = 25,
                          total.gall.abund = mean(all.measures.1000.4.NAone$total.gall.abund))
predict(gam.1000.4.cov, pred.25geno)/predict(gam.1000.4.cov, pred.1geno) # 13% increase in food-web complexity
predict(gam.1000.4.cov, pred.8geno)/predict(gam.1000.4.cov, pred.1geno) # 11% increase in food-web complexity from 1 to 8 genotypes

summary(filter(all.measures.1000.4, genotype.richness == 1)$total_complexity)
summary(filter(all.measures.1000.4, genotype.richness == 25)$total_complexity)

#plot(gam.1000.4.cov, residuals = TRUE)

#QuantPsyc::lm.beta(lm.1000.4.cov)
#car::vif(lm.1000.4.cov)
#visreg(lm.1000.4.cov)
#predict(lm.1000.4.cov, predict.data)
#2.133265/1.82181 # 17%
#visreg(lm.1000.4.cov, xvar = "genotype.richness", jitter = TRUE, col = "grey", pch = 21)

## plot relationship between willow replicates and food-web complexity for monocultures. Seems to plateau as long as I remove NAs from food-web calculation
ggplot(all.measures.mono,
       aes(x = reps.sampled, y = total_complexity)) +
  geom_jitter(shape = 1, color = "grey", position = position_jitter(width = 0.15, height = NULL), size = 1) + 
  stat_summary(fun.y = mean, geom = "line", color = "steelblue", size = 2) +
  stat_summary(fun.y = mean, geom = "point", color = "black", size = 3) +
  xlab("No. of willow replicates\nfor genotype monocultures") + 
  ylab(bquote('Food-web complexity ('*italic(LD[q])*')')) +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 9),#10
        axis.text.x = element_text(size = 9),#10
        axis.title.x = element_text(size = 11, vjust = 0.1),#vjust = 0.1, 
        axis.title.y = element_text(size = 11, vjust = 0.5),#vjust = 1, 
        panel.grid = element_blank()) 

# setup multipanel plot for combining ggplot and base graphics (ordination from "link_composition_plots")
#A.total <- data.frame(x = 0.75, y = 2.5, labels = "(A)")

## plot the relationship between genetic variation and total food web complexity (weighted linkage density)
total <- ggplot(all.measures.2500.4, 
       aes(x = genotype.richness, y = total_complexity)) + 
  geom_jitter(shape = 1, color = "grey", position = position_jitter(width = 0.15, height = NULL), size = 1) + 
  stat_summary(fun.y = mean, geom = "point", color = "steelblue", size = 3) + 
  #geom_hline(yintercept = single.complexity.max, linetype = "dashed") +
  xlab("No. of willow genotypes") + 
  ylab(bquote('Food-web complexity ('*italic(LD[q])*')')) +
  scale_x_continuous(limits = c(0.5,25), 
                     breaks = c(1,5,10,15,20,25)) +
  scale_y_continuous(limits = c(1,2.5), #c(1,2.4)
                     breaks = seq(1, 2.5, by = 0.5)) + #seq(1, 2.5, by = 0.25)
  #geom_text(data = A.total, aes(x = x, y = y, label = labels), 
  #          inherit.aes = FALSE, size = 4) +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 9),#10
        axis.text.x = element_text(size = 9),#10
        axis.title.x = element_text(size = 11, vjust = 0.1),#vjust = 0.1, 
        axis.title.y = element_text(size = 11, vjust = 0.5),#vjust = 1, 
        panel.grid = element_blank()) 


# then plot ordination in "link_composition_plots.R" and save the figure as a pdf, portrait, 3.42" x 6"

## plot the relationship between genetic variation and willow-gall complexity (weighted linkage density)
# weird pattern emerging...
ggplot(all.measures.2500.4, 
       aes(x = genotype.richness, y = link.density.plant_gall)) + 
  geom_jitter(shape = 1, color = "grey", position = position_jitter(width = 0.15, height = NULL)) + 
  stat_summary(fun.y = mean, geom = "point", color = "steelblue", size = 8) + 
  xlab("Genetic variation (no. of genotypes)") + ylab("Food web complexity index") +
  scale_x_continuous(limits = c(0,25), breaks = 1:25) +
  scale_y_continuous(limits = c(1,2.5), breaks = seq(1,2.5, by = 0.25)) +
  theme_bw() + 
  theme(axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, vjust = 0.1),
        axis.title.y = element_text(size = 16, vjust = 1),
        panel.grid = element_blank())

## plot the relationship between genetic variation and gall-ptoid complexity (weighted linkage density)
ggplot(all.measures.2500.4, 
       aes(x = genotype.richness, y = link.density.gall_ptoid)) + 
  geom_jitter(shape = 1, color = "grey", position = position_jitter(width = 0.15, height = NULL)) + 
  stat_summary(fun.y = mean, geom = "point", color = "steelblue", size = 8) + 
  xlab("Genetic variation (no. of genotypes)") + ylab("Food web complexity index") +
  scale_x_continuous(limits = c(0,25), breaks = 1:25) +
  scale_y_continuous(limits = c(1,2.85), breaks = seq(1,2.75, by = 0.25)) +
  theme_bw() + 
  theme(axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, vjust = 0.1),
        axis.title.y = element_text(size = 16, vjust = 1),
        panel.grid = element_blank())

## calculate the mean food web complexity for the different levels of genetic variation
web.measures.summary <- all.measures.2500.4 %>%
  group_by(genotype.richness) %>%
  summarise(mean.complexity = mean(total_complexity, na.rm = TRUE))
with(web.measures.summary, max(mean.complexity)/min(mean.complexity)) # 52% increase in food-web complexity
#with(web.measures.summary, max(mean.complexity)/single.complexity.max) # 8% higher due to complimentarity effects

#max(all.measures$total_complexity, na.rm = TRUE)/single.complexity.max # up to 15% higher due to complimentarity effects.

## linear model with just the mean data. 
#web.sum.lm <- lm(log(mean.complexity) ~ log(Number.of.Genotypes), data = web.measures.summary)
#summary(web.sum.lm) # Need to understand interpretation here.
#visreg(web.sum.lm) # not quite there

## non-linear model with just the mean data
#require(mgcv)
#web.sum.gam <- gam(mean.complexity ~ s(Number.of.Genotypes), data = web.measures.summary)
#summary(web.sum.gam)
#visreg(web.sum.gam)


## linear model with all of the data. 
#web.lm <- lm(log(total_complexity) ~ log(Number.of.Genotypes), data = all.measures)
#summary(web.lm)
#visreg(web.lm)

## non-linear model with all of the data
#web.gam <- gam(total_complexity ~ s(Number.of.Genotypes), data = all.measures)
#summary(web.gam)
#visreg(web.gam)

## null model expectation for food-web complexity
#mean(filter(all.measures, Number.of.Genotypes == 1)$total_complexity, na.rm = TRUE) # does this represent the expectation for food-web complexity under the sampling effect?

## piecewise Structural Equation Model
all.measures.1000.4.sub <- all.measures.1000.4 %>%
  filter(total.gall_ptoid.abund > 0) %>%
  mutate(log.gall.abund = log(total.gall.abund),
         log.gall_ptoid.abund = log(total.gall_ptoid.abund),
         log.genotype.rich = log(genotype.richness))

## full mediation model
plot(log.gall.abund ~ log.genotype.rich, all.measures.1000.4.sub)
lm.gall.abund <- lm(log.gall.abund ~ log.genotype.rich + I(log.genotype.rich^2), 
                    na.action = na.omit, 
                    data = all.measures.1000.4.sub)
summary(lm.gall.abund)

plot(log.gall_ptoid.abund ~ log.gall.abund, all.measures.1000.4.sub)
lm.gall_ptoid.abund <- lm(log.gall_ptoid.abund ~ log.gall.abund +
                            I(log.gall.abund^2), 
                          na.action = na.omit, 
                          data = all.measures.1000.4.sub)
summary(lm.gall_ptoid.abund)

plot(total_complexity ~ log.gall_ptoid.abund, all.measures.1000.4.sub)
plot(total_complexity ~ log.gall.abund, all.measures.1000.4.sub)
lm.total_complexity <- lm(total_complexity ~ log.gall_ptoid.abund +
                            I(log.gall_ptoid.abund^2) +
                            log.gall.abund + I(log.gall.abund^2),
                          na.action = na.omit, 
                          data = all.measures.1000.4.sub)
summary(lm.total_complexity)

LDq.fullmed <- list(lm.gall.abund, lm.gall_ptoid.abund, 
                   lm.total_complexity)
sem.fit(LDq.fullmed, data = all.measures.1000.4.sub) # suggests that the primary pathway missing is between genotype.richness and total complexity.
#sem.model.fits(LDq.fullmed)
sem.coefs(LDq.fullmed, all.measures.1000.4.sub, standardize = "scale")
#lav.LDq.fullmed <- sem.lavaan(LDq.fullmed, all.measures.2500.4)

## new model, full mediation plus pathway between gentoype.richness and total_complexity as well as total.gall.abund and total_complexity
plot(total_complexity ~ log.genotype.rich, all.measures.1000.4.sub)
lm.total_complexity.geno <- lm(total_complexity ~ 
                                 log.gall_ptoid.abund + 
                                 log.gall.abund + 
                                 log.genotype.rich,
                               na.action = na.omit, data = all.measures.1000.4.sub)
summary(lm.total_complexity.geno)

#plot(log.gall_ptoid.abund ~ log.genotype.rich, all.measures.1000.4.sub)
#lm.gall_ptoid.geno <- lm(log.gall_ptoid.abund ~ log.gall.abund + 
 #                                log.genotype.rich,
  #                             na.action = na.omit, data = all.measures.1000.4.sub)
#summary(lm.gall_ptoid.geno)

LDq.fullmed.geno <- list(lm.gall.abund, lm.gall_ptoid.abund, 
                   lm.total_complexity.geno)
sem.fit(LDq.fullmed.geno, data = all.measures.1000.4.sub) # AIC = 24
sem.model.fits(LDq.fullmed.geno)
sem.coefs(LDq.fullmed.geno, all.measures.1000.4.sub, standardize = "scale")

## plot structural equation model
lav.LDq.fullmed.geno <- sem.lavaan(LDq.fullmed.geno, all.measures.1000.4.sub)

labels = c("Total gall\nabundance", "Frequency of\ngall-parasitoid\ninteractions", "Food-web\ncomplexity (LDq)", "No. of willow\ngenotypes")
layout.LD = matrix(c(-0.25, -0.5,
                     0.25, -0.5,
                     0.5, 0.5,
                     -0.5, 0.5), 
                   ncol = 2, byrow = TRUE)

semPaths(lav.LDq.fullmed.geno, 
         what = "std", 
         layout = layout.LD, 
         residuals = FALSE,
         nCharNodes = 0,
         nodeLabels = labels,
         sizeMan = 10,
         edge.label.cex = 1.5,
         label.prop = 1.9,
         posCol = c("blue","red"),
         fade = FALSE)

## old no longer going to use probably 
## plot updated Gravel simulation
#all.measures.4 <- filter(all.measures, reps.sampled == 4)
#all.measures.NAzero <- all.measures
#all.measures.NAzero[is.na(all.measures.NAzero)] <- 0

#all.measures.NAone <- all.measures
#all.measures.NAone[is.na(all.measures.NAone)] <- 1

# check plant sampling number
ggplot(all.measures.2500.4,
       aes(x = genotype.richness, y = total_complexity)) +
  geom_point(color = "grey", shape = 1) +
  geom_smooth(size = 2)

ggplot(all.measures,#,filter( plants_sampled %in% c(3,4,6,9,9,10,12)),#.NAone,
       aes(x = plants_sampled, y = total_complexity, group = genotype.richness, color = genotype.richness)) + # genotype.richness
  #facet_wrap(~reps.sampled, nrow = 2) +
  geom_point(color = "grey", shape = 1) +
  geom_smooth(se = FALSE)

plot(log(total_complexity) ~ log(plants_sampled), all.measures)
plot(log(total_complexity) ~ log(genotype.richness), all.measures)
#plot(total_complexity ~ functional.divergence, all.measures)
plot(plants_sampled ~ genotype.richness, all.measures)
#plot(functional.divergence ~ plants_sampled, all.measures)
#plot(functional.divergence ~ genotype.richness, all.measures)

lm.1 <- lm(total_complexity ~ log(plants_sampled) + log(genotype.richness), all.measures) #+ log(genotype.richness)
summary(lm.1)
#plot(lm.1)
visreg(lm.1)
#visreg(lm.1, xvar = "genotype.richness", by = "plants_sampled")
car::vif(lm.1) # still a bit of collinearity...but no collinearity between functional divergence and plants_sampled...
#visreg(lm.1, xvar = "genotype.richness", by = "plants_sampled")


