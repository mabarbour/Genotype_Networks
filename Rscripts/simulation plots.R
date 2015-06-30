## load required libraries
require(visreg)
require(ggplot2)
require(grid)
require(dplyr)
require(piecewiseSEM)
require(semPlot)

## upload results from simulation
all.measures <- read.csv("~/Documents/Genotype_Networks/data/food web complexity simulation Gravel.csv")
#all.measures <- read.csv("~/Documents/Genotype_Networks/data/food web complexity simulation.csv")
#dim(all.measures)[1] # 2221 unique simulations. 2425 simulations originally run
#table(all.measures$Number.of.Genotypes)

#single.complexity.max <- max(filter(all.measures, Number.of.Genotypes == 1)$total_complexity, na.rm = TRUE)

## plot updated Gravel simulation
#all.measures.4 <- filter(all.measures, reps.sampled == 4)
all.measures.NAzero <- all.measures
all.measures.NAzero[is.na(all.measures.NAzero)] <- 0

all.measures.NAone <- all.measures
all.measures.NAone[is.na(all.measures.NAone)] <- 1

# check plant sampling number
ggplot(all.measures,#,filter( plants_sampled %in% c(3,4,6,9,9,10,12)),#.NAone,
       aes(x = genotype.richness, y = total_complexity, group = reps.sampled, color = reps.sampled)) + # genotype.richness
  #facet_wrap(~reps.sampled, nrow = 2) +
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


# setup multipanel plot for combining ggplot and base graphics (ordination from "link_composition_plots")
#A.total <- data.frame(x = 0.75, y = 2.5, labels = "(A)")

## plot the relationship between genetic variation and total food web complexity (weighted linkage density)
total <- ggplot(all.measures, 
       aes(x = Number.of.Genotypes, y = total_complexity)) + 
  geom_jitter(shape = 1, color = "grey", position = position_jitter(width = 0.15, height = NULL), size = 1) + 
  stat_summary(fun.y = mean, geom = "point", color = "steelblue", size = 3) + 
  geom_hline(yintercept = single.complexity.max, linetype = "dashed") +
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
ggplot(all.measures, 
       aes(x = Number.of.Genotypes, y = link.density.plant_gall)) + 
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
ggplot(all.measures, 
       aes(x = Number.of.Genotypes, y = link.density.gall_ptoid)) + 
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
web.measures.summary <- all.measures %>%
  group_by(Number.of.Genotypes) %>%
  summarise(mean.complexity = mean(total_complexity, na.rm = TRUE))
with(web.measures.summary, max(mean.complexity)/min(mean.complexity))
with(web.measures.summary, max(mean.complexity)/single.complexity.max) # 8% higher due to complimentarity effects

max(all.measures$total_complexity, na.rm = TRUE)/single.complexity.max # up to 15% higher due to complimentarity effects.

## linear model with just the mean data. 
web.sum.lm <- lm(log(mean.complexity) ~ log(Number.of.Genotypes), data = web.measures.summary)
summary(web.sum.lm) # Need to understand interpretation here.
visreg(web.sum.lm) # not quite there

## non-linear model with just the mean data
require(mgcv)
web.sum.gam <- gam(mean.complexity ~ s(Number.of.Genotypes), data = web.measures.summary)
summary(web.sum.gam)
visreg(web.sum.gam)


## linear model with all of the data. 
web.lm <- lm(log(total_complexity) ~ log(Number.of.Genotypes), data = all.measures)
summary(web.lm)
visreg(web.lm)

## non-linear model with all of the data
web.gam <- gam(total_complexity ~ s(Number.of.Genotypes), data = all.measures)
summary(web.gam)
visreg(web.gam)

## null model expectation for food-web complexity
mean(filter(all.measures, Number.of.Genotypes == 1)$total_complexity, na.rm = TRUE) # does this represent the expectation for food-web complexity under the sampling effect?

## piecewise Structural Equation Model

## full mediation model
lm.gall.abund <- lm(total.gall.abund ~ genotype.richness, 
                    na.action = na.omit, data = all.measures)
summary(lm.gall.abund)


lm.gall_ptoid.abund <- lm(total.gall_ptoid.abund ~ total.gall.abund, 
                          na.action = na.omit, data = all.measures)
summary(lm.gall_ptoid.abund)


lm.total_complexity <- lm(total_complexity ~ total.gall_ptoid.abund,
                          na.action = na.omit, data = all.measures)
summary(lm.total_complexity)

LDq.fullmed = list(lm.gall.abund, lm.gall_ptoid.abund, 
                   lm.total_complexity)
sem.fit(LDq.fullmed, data = all.measures) # suggests that the primary pathway missing is between genotype.richness and total complexity.
sem.model.fits(LDq.fullmed)
sem.coefs(LDq.fullmed, all.measures, standardize = "scale")
lav.LDq.fullmed <- sem.lavaan(LDq.fullmed, all.measures)

## new model, full mediation plus pathway between gentoype.richness and total_complexity
lm.total_complexity.geno <- lm(total_complexity ~ total.gall_ptoid.abund + genotype.richness,
                               na.action = na.omit, data = all.measures)
summary(lm.total_complexity.geno)

LDq.fullmed.geno = list(lm.gall.abund, lm.gall_ptoid.abund, 
                   lm.total_complexity.geno)
sem.fit(LDq.fullmed.geno, data = all.measures)
sem.model.fits(LDq.fullmed.geno)
sem.coefs(LDq.fullmed.geno, all.measures, standardize = "scale")

## plot structural equation model
lav.LDq.fullmed.geno <- sem.lavaan(LDq.fullmed.geno, all.measures)

labels = c("Total gall\nabundance", "Frequency of\ngall-parasitoid\ninteractions", "Weighted\nlinkage density", "No. of willow\ngenotypes")

semPaths(lav.LDq.fullmed.geno, 
         what = "std", 
         layout = "tree", 
         rotation = 2,
         nCharNodes = 0,
         nodeLabels = labels,
         sizeMan = 10,
         edge.label.cex = 1.5,
         label.prop = 1.5,
         posCol = c("blue","red"))
