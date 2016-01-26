
## load required libraries
#require(visreg)
require(ggplot2)
require(grid)
require(dplyr)
require(tidyr)


## upload results from full simulation
all.measures <- read.csv("manuscript/Dryad_data_copy/simulation_data_output/food web complexity simulation 50 reps of 200 sims 4 reps.csv")
dim(all.measures)[1] # 216,102 unique simulations. 250,00 simulations originally run

## upload results from monoculture simulations
mono1 <- read.csv("manuscript/Dryad_data_copy/simulation_data_output/food web complexity simulation 1000 monosims 1 reps.csv")
mono2 <- read.csv("manuscript/Dryad_data_copy/simulation_data_output/food web complexity simulation 1000 monosims 2 reps.csv")
mono3 <- read.csv("manuscript/Dryad_data_copy/simulation_data_output/food web complexity simulation 1000 monosims 3 reps.csv")
mono4 <- read.csv("manuscript/Dryad_data_copy/simulation_data_output/food web complexity simulation 1000 monosims 4 reps.csv")

monos <- bind_rows(list(mono1, mono2, mono3, mono4)) #%>%

ggplot(monos, aes(x = plants.sampled, y = total_complexity, group = sim.number)) +
  #geom_jitter(shape = 1, color = "grey", position = position_jitter(width = 0.15, height = NULL)) + 
  stat_summary(fun.y = mean, geom = "point", color = "grey", size = 2, shape = 1) #+
  #geom_smooth(method = "lm", formula = y ~ log(x), color = "black")

monos.summary <- monos %>% group_by(sim.number, plants.sampled) %>% summarise(total_complexity.mean = mean(total_complexity, na.rm = TRUE))

monos.summary.plot <- monos.summary %>% group_by(plants.sampled) %>%
  summarise(total_complexity.mean = mean(total_complexity.mean))

## summary stats from full simulation
alpha.summary <- all.measures %>% 
  filter(genotypes.sampled == 1) %>%
  group_by(df.sim.number) %>%
  summarise(alpha.mean = mean(total_complexity, na.rm = TRUE),
            alpha.min = 1) %>% # minimum value that food-web complexity can take
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
  as.matrix(.)

## calculate difference between maximum monoculture and 25-genotype polyculture for all different data sets of the simulation
agg.prop <- all.measures %>%
  filter(genotypes.sampled %in% c(1,25)) %>%
  group_by(df.sim.number, genotypes.sampled) %>%
  summarise(max.complex = max(total_complexity, na.rm = TRUE)) %>%
  spread(genotypes.sampled, max.complex) 

diff.25.max1 <- agg.prop$'25' - agg.prop$'1'
summary(diff.25.max1)
hist(diff.25.max1)

## calculate the mean food web complexity for the different levels of genetic variation
web.measures.summary <- all.measures %>% # all.measures
  group_by(df.sim.number, genotypes.sampled) %>%
  summarise(mean.complexity = mean(total_complexity, na.rm = TRUE))

mean.complexity.25 <- mean(filter(web.measures.summary, genotypes.sampled == 25)$mean.complexity)

web.measures.overall <- web.measures.summary %>%
  group_by(genotypes.sampled) %>%
  summarise(mean.complexity = mean(mean.complexity))

## identify starting values for Micheaelis-Menton function with a linear model. Subtracting 1 shifts the function to the right. I did this because the starting value is 1 genotype. Slope and intercept will correspond to the parameters a and d, respectively for the non-linear model
b <-2 # increased progressively from a value of 1 to maximize R2. 2 without switching gall_ptoid LD NA to zero.
full.init.lm <- lm(mean.complexity ~ I((genotypes.sampled - 1)/(genotypes.sampled - 1 + b)), web.measures.summary)
summary(full.init.lm) # new R2 = 0.9573

## non-linear model
full.nls <- nls(mean.complexity ~ a*(genotypes.sampled - 1)/(genotypes.sampled - 1 + b) + d, web.measures.summary,
                start = list(a = coef(full.init.lm)[2], 
                             b = b, 
                             d = coef(full.init.lm)[1]))
summary(full.nls)
full.nls.predict.df <- data.frame(mean.complexity = predict(full.nls, newdata = data.frame(genotypes.sampled = 1:25)),
                                  genotypes.sampled = 1:25)
max.full.nls <- max(predict(full.nls, newdata = data.frame(genotypes.sampled = 1:25))) # predict food-web complexity of 2.2
(max.full.nls - mean.complexity.25)/mean.complexity.25 # deviated by less than a 10th of 1% of observed value (over predicted)

pred.full.nls.100genos <- max(predict(full.nls, newdata = data.frame(genotypes.sampled = 1:100))) # hypothetical prediction of complexity of 100-genotype food-web

full.nls.lm <- lm(mean.complexity ~ I((genotypes.sampled - 1)/(genotypes.sampled - 1 + coef(full.nls)[2])), web.measures.summary)
summary(full.nls.lm) # R2 = 0.96
plot(full.nls.lm) # residuals look good
hist(residuals(full.nls.lm))
#plot(fitted(full.nls.lm) ~ web.measures.summary$genotypes.sampled)

## Testing non-asymptotic models for accumulation of food-web complexity

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

## non-linear model based on first 4 data points. Used values from full non-linear model as initial values in this model
full.nls.4pts <- nls(mean.complexity ~ a*(genotypes.sampled - 1)/(genotypes.sampled - 1 + b) + d, 
                     filter(web.measures.summary, genotypes.sampled < 5),
                     start = list(a = coef(full.nls)["a"], 
                                  b = coef(full.nls)["b"], 
                                  d = coef(full.nls)["d"]))
summary(full.nls.4pts)
max.4pts.nls <- max(predict(full.nls.4pts, newdata = data.frame(genotypes.sampled = 1:25)))
(max.4pts.nls - mean.complexity.25)/max.4pts.nls # error is still less than 10th of 1%. used to overestimates by about 2%

pred.4pts.nls.100genos <- max(predict(full.nls.4pts, newdata = data.frame(genotypes.sampled = 1:100)))

#full.nls.lm.4pts <- lm(mean.complexity ~ I((genotypes.sampled - 1)/(genotypes.sampled - 1 + coef(full.nls.4pts)[2])), filter(web.measures.summary, genotypes.sampled < 5))
#summary(full.nls.lm.4pts) # note that this R2 doesn't mean much, because it is being fit to a smaller dataset.

## find initial starting values for monoculture samples using a linear model
b <- 4 # progressively increased from 1 to maximize R2. Used to be 4
monos.init.lm <- lm(total_complexity.mean ~ I((plants.sampled - 1)/(plants.sampled - 1 + b)), monos.summary)
summary(monos.init.lm) # R2 = 0.8778

## non-linear model for monocultures. Used the coefficients from the full model as initial starting values so that the methods were comparable with the subsampling of 4 genotypes. Note however, that the coefficient estimates were robust to different starting value (e.g. based on prior iterations with linear model)
mono.nls <- nls(total_complexity.mean ~ a*(plants.sampled - 1)/(plants.sampled - 1 + b) + d, monos.summary,
                start = list(a = coef(full.nls)["a"], 
                             b = coef(full.nls)["b"], 
                             d = coef(full.nls)["d"]))
                        #list(a = coef(monos.init.lm)[2], 
                         #    b = b, 
                          #   d = coef(monos.init.lm)[1]))
summary(mono.nls)
monos.nls.predict <- data.frame(total_complexity = predict(mono.nls, newdata = data.frame(plants.sampled = 1:100)),
                                plants.sampled = 1:100)
max.mono.nls <- max(monos.nls.predict$total_complexity) # predict food-web complexity of 1.84 at 100 plants (comparable to 25 genotypes)

mono.nls.lm <- lm(total_complexity.mean ~ I((plants.sampled - 1)/(plants.sampled - 1 + coef(mono.nls)[2])), monos.summary)
summary(mono.nls.lm)

## log-log
mono.loglog.lm <- lm(log(total_complexity.mean) ~ log(plants.sampled), monos.summary)
summary(mono.loglog.lm)
exp(predict(mono.loglog.lm, newdata = data.frame(plants.sampled = 1:100))) # 2.45

## log-linear
mono.loglin.lm <- lm(total_complexity.mean ~ log(plants.sampled), monos.summary)
summary(mono.loglin.lm)
predict(mono.loglin.lm, newdata = data.frame(plants.sampled = 1:100)) # 2.17

## magnitude of increase in food-web complexity due to genotypes hosting distinct sets of trophic interactions
mean.complexity.25/max.mono.nls # 20% increase in food-web complexity.

## additive partition of food-web complexity
gamma.contrib <- mean.complexity.25 - 1 # subtracted 1 because this is the minimum value of food-web complexity for including in our simulation
alpha.contrib <- final.summary[ ,"alpha.mean_mean"] - 1
beta.area.contrib <- max.mono.nls - final.summary[ ,"alpha.mean_mean"]
beta.genetic.contrib <- mean.complexity.25 - max.mono.nls
beta.genetic.contrib/gamma.contrib
beta.area.contrib/gamma.contrib
alpha.contrib/gamma.contrib


# setup multipanel plot for combining ggplot and base graphics (ordination from "link_composition_plots")

## plot the relationship between genetic variation and total food web complexity (weighted linkage density)
total.labels <- paste(c(" ",1,5,10,15,20,25),
                      c("(1) "," (4)","(20)","(40)","(60)","(80)","(100)"), sep = "\n")
total <- ggplot(web.measures.summary,#all.measures, 
       aes(x = genotypes.sampled, 
           y = mean.complexity,# total_complexity, 
           group = genotypes.sampled)) + 
  #stat_summary(fun.y = mean, geom = "point", color = "grey",
   #            shape = 1, size = 1.5) + 
  geom_jitter(color = "grey", shape = 1, position = position_jitter(width = 0.25, height = NULL), size = 1) +
  stat_summary(fun.y = mean, geom = "point", color = "steelblue",
                           shape = 20, size = 3) +
  #geom_line(data = full.nls.predict.df, aes(x = genotypes.sampled, y = mean.complexity),
           # size = 1,
           # color = "steelblue", inherit.aes = FALSE) +
  #geom_boxplot(aes(group = genotypes.sampled), 
    #           notch = TRUE, outlier.shape = 1, 
     #          outlier.colour = "grey", fill = "steelblue") +
  #geom_jitter(shape = 1, 
   #           color = "grey", 
    #          position = position_jitter(width = 0.25, height = NULL), 
     #         size = 1) + 
  #stat_summary(fun.y = mean, geom = "point", color = "steelblue", size = 3) + 
  #geom_boxplot(data = monos.summary, 
   #          aes(x = plants.sampled/4, 
    #             y = total_complexity.mean,
     #            group = plants.sampled), 
      #       color = "black", 
       #      outlier.shape = NA, 
             #outlier.colour = "black", 
        #     inherit.aes = FALSE) +
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


# then plot ordination in "link_composition_plots.R" and save the figure as a pdf, portrait, 3.42" x 4"

## plot the relationship between genetic variation and willow-gall complexity (weighted linkage density)
# weird pattern emerging...
ggplot(all.measures, 
       aes(x = genotypes.sampled, y = link.density.plant_gall)) + 
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
       aes(x = genotypes.sampled, y = link.density.gall_ptoid)) + 
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

## nestedness vs. modularity

nested <- matrix(c(1,1,1,1,
                   1,1,1,0,
                   1,1,0,0,
                   1,0,0,0), ncol = 4, byrow = TRUE)
modular <- matrix(c(1,1,0,0,
                    1,1,0,0,
                    0,0,1,1,
                    0,0,1,1), ncol = 4, byrow = TRUE)
mean(rowSums(nested > 0))
mean(rowSums(modular > 0))
