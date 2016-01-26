#############################################################
#  Description: This script plots the simulation data in Fig. 6 and reproduces the output in Tables S4 and S5 of supplementary material for the manuscript, "Genetic specificity of a plant-insect food web: Implications for linking genetic variation to food-web complexity"
#  Code author: Matthew A. Barbour
#  Email: barbour@zoology.ubc.ca
#############################################################

## load required libraries ----
library(ggplot2)
library(grid)
library(dplyr)
library(tidyr)

## upload results from full simulation ----
all.measures <- read.csv("manuscript/Dryad_data_copy/simulation_data_output/food web complexity simulation 50 reps of 200 sims 4 reps.csv")
dim(all.measures)[1] # 216,102 unique simulations. 250,00 simulations originally run

## upload results from monoculture simulations ----
mono1 <- read.csv("manuscript/Dryad_data_copy/simulation_data_output/food web complexity simulation 1000 monosims 1 reps.csv")
mono2 <- read.csv("manuscript/Dryad_data_copy/simulation_data_output/food web complexity simulation 1000 monosims 2 reps.csv")
mono3 <- read.csv("manuscript/Dryad_data_copy/simulation_data_output/food web complexity simulation 1000 monosims 3 reps.csv")
mono4 <- read.csv("manuscript/Dryad_data_copy/simulation_data_output/food web complexity simulation 1000 monosims 4 reps.csv")

monos <- bind_rows(list(mono1, mono2, mono3, mono4)) 

## calculate the mean food web complexity for the different levels of genetic variation for each randomly sampled dataset in the simulation. ----
web.measures.summary <- all.measures %>% # all.measures
  group_by(df.sim.number, genotypes.sampled) %>%
  summarise(mean.complexity = mean(total_complexity, na.rm = TRUE))

## average complexity of 25-genotype (100-plant) food web ----
mean.complexity.25 <- mean(filter(web.measures.summary, genotypes.sampled == 25)$mean.complexity)

## Table S5 ----

## Asymptotic model
# identify starting values for Micheaelis-Menton function with a linear model. Subtracting 1 shifts the function to the right. I did this because the starting value is 1 genotype. Slope and intercept will correspond to the parameters a and d, respectively for the non-linear model
b <-2 # increased progressively from a value of 1 to maximize R2. 
full.init.lm <- lm(mean.complexity ~ I((genotypes.sampled - 1)/(genotypes.sampled - 1 + b)), web.measures.summary)
summary(full.init.lm) # R2 = 0.9573

# non-linear model
full.nls <- nls(mean.complexity ~ a*(genotypes.sampled - 1)/(genotypes.sampled - 1 + b) + d, web.measures.summary,
                start = list(a = coef(full.init.lm)[2], 
                             b = b, 
                             d = coef(full.init.lm)[1]))
summary(full.nls)
full.nls.predict.df <- data.frame(mean.complexity = predict(full.nls, newdata = data.frame(genotypes.sampled = 1:25)),
                                  genotypes.sampled = 1:25)
max.full.nls <- max(predict(full.nls, newdata = data.frame(genotypes.sampled = 1:25))) # predict food-web complexity of 2.2
(max.full.nls - mean.complexity.25)/mean.complexity.25 # deviated by less than a 10th of 1% of observed value (over predicted)

# rerun non-linear model as linear model to calculate R2
full.nls.lm <- lm(mean.complexity ~ I((genotypes.sampled - 1)/(genotypes.sampled - 1 + coef(full.nls)[2])), web.measures.summary)
summary(full.nls.lm) # R2 = 0.96
plot(full.nls.lm) # residuals look good
hist(residuals(full.nls.lm))

## Non-asymptotic models

# log-linear
full.lm.loglin <- lm(mean.complexity ~ log(genotypes.sampled), web.measures.summary)
summary(full.lm.loglin)
max(predict(full.lm.loglin, newdata = data.frame(genotypes.sampled = 1:25)))
(2.262 - 2.209)/2.209
visreg::visreg(full.lm.loglin)
plot(full.lm.loglin) # horrible residuals

# log-log
full.lm.loglog <- lm(log(mean.complexity) ~ log(genotypes.sampled), web.measures.summary)
summary(full.lm.loglog)
exp(max(predict(full.lm.loglog, newdata = data.frame(genotypes.sampled = 1:25))))
(2.277 - 2.209)/2.209
visreg::visreg(full.lm.loglog)
plot(full.lm.loglog) # horrible residuals again

## Table S4 ----

## Asymptotic model
# find initial starting values for monoculture samples using a linear model
b <- 4 # progressively increased from 1 to maximize R2.
monos.init.lm <- lm(total_complexity.mean ~ I((plants.sampled - 1)/(plants.sampled - 1 + b)), monos.summary)
summary(monos.init.lm) # R2 = 0.8778

# non-linear model for monocultures. Used the coefficients from the full model as initial starting values so that the methods were comparable with the subsampling of 4 genotypes. Note however, that the coefficient estimates were robust to different starting value (e.g. based on prior iterations with linear model)
mono.nls <- nls(total_complexity.mean ~ a*(plants.sampled - 1)/(plants.sampled - 1 + b) + d, monos.summary,
                start = list(a = coef(full.nls)["a"], 
                             b = coef(full.nls)["b"], 
                             d = coef(full.nls)["d"]))
summary(mono.nls)
monos.nls.predict <- data.frame(total_complexity = predict(mono.nls, newdata = data.frame(plants.sampled = 1:100)),
                                plants.sampled = 1:100)
max.mono.nls <- max(monos.nls.predict$total_complexity) # predict food-web complexity of 1.84 at 100 plants (comparable to 25 genotypes with 4 plants for each genotype)

# rerun non-linear model as linear model to calculate R2
mono.nls.lm <- lm(total_complexity.mean ~ I((plants.sampled - 1)/(plants.sampled - 1 + coef(mono.nls)[2])), monos.summary)
summary(mono.nls.lm) # R2 = 0.885

## Non-asymptotic models

# log-log
mono.loglog.lm <- lm(log(total_complexity.mean) ~ log(plants.sampled), monos.summary)
summary(mono.loglog.lm) # R2 = 0.881
exp(predict(mono.loglog.lm, newdata = data.frame(plants.sampled = 1:100))) # 2.45

# log-linear
mono.loglin.lm <- lm(total_complexity.mean ~ log(plants.sampled), monos.summary)
summary(mono.loglin.lm) # R2 = 0.884
predict(mono.loglin.lm, newdata = data.frame(plants.sampled = 1:100)) # 2.17

## Determining accuracy of asymptotic model ----
# non-linear model based on first 4 data points. Used values from full non-linear model as initial values in this model
full.nls.4pts <- nls(mean.complexity ~ a*(genotypes.sampled - 1)/(genotypes.sampled - 1 + b) + d, 
                     filter(web.measures.summary, genotypes.sampled < 5),
                     start = list(a = coef(full.nls)["a"], 
                                  b = coef(full.nls)["b"], 
                                  d = coef(full.nls)["d"]))
summary(full.nls.4pts)
max.4pts.nls <- max(predict(full.nls.4pts, newdata = data.frame(genotypes.sampled = 1:25)))
(max.4pts.nls - mean.complexity.25)/max.4pts.nls # error is still less than 10th of 1%. used to overestimates by about 2%

pred.4pts.nls.100genos <- max(predict(full.nls.4pts, newdata = data.frame(genotypes.sampled = 1:100)))

## Food-web complexity increased by 20% due to genotypes hosting distinct sets of trophic interactions ----
mean.complexity.25/max.mono.nls # 20% increase in food-web complexity

## Fig. 6: plot of simulation data ----

# plot the relationship between genetic variation and total food web complexity (weighted linkage density)
total.labels <- paste(c(" ",1,5,10,15,20,25),
                      c("(1) "," (4)","(20)","(40)","(60)","(80)","(100)"), sep = "\n") # x-axis labels

total <- ggplot(web.measures.summary,#all.measures, 
       aes(x = genotypes.sampled, 
           y = mean.complexity,
           group = genotypes.sampled)) + 
  geom_jitter(color = "grey", shape = 1, position = position_jitter(width = 0.25, height = NULL), size = 1) +
  stat_summary(fun.y = mean, geom = "point", color = "steelblue",
                           shape = 20, size = 3) +
  geom_point(data = monos.summary.plot, 
             aes(x = plants.sampled/4, 
                 y = total_complexity.mean), 
             color = "black", size = 1.5,
             inherit.aes = FALSE) +
  geom_line(data = monos.nls.predict, 
            aes(x = plants.sampled/4, y = total_complexity),
            linetype = "dashed",
            color = "black",
            inherit.aes = FALSE) +
  xlab("No. of willow genotypes (no. of plants)") + 
  ylab(bquote('Food-web complexity ('*italic(LD[q])*')')) +
  scale_x_continuous(limits = c(0.1,25.5), 
                     breaks = c(0.25,1,5,10,15,20,25),
                     labels = total.labels) +
  scale_y_continuous(limits = c(1,2.45), #c(1,2.4)
                     breaks = seq(1, 2.25, by = 0.25)) +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 6),#9
        axis.title.x = element_text(size = 11, vjust = 0.1),
        axis.title.y = element_text(size = 11, vjust = 0.5),
        panel.grid = element_blank()) 

# then plot ordination in "link_composition_plots.R" to complete Fig. 6

