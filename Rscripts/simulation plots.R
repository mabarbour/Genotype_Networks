
## load required libraries
#require(visreg)
require(ggplot2)
require(grid)
require(dplyr)
#require(piecewiseSEM)
#require(semPlot)

## upload results from simulation
all.measures <- read.csv("~/Documents/Genotype_Networks/data/food web complexity simulation 40 reps of 100 sims.csv")
dim(all.measures)[1] # 88605 unique simulations. 100,00 simulations originally run
table(all.measures$genotypes.sampled)

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

# setup multipanel plot for combining ggplot and base graphics (ordination from "link_composition_plots")

## plot the relationship between genetic variation and total food web complexity (weighted linkage density)
total <- ggplot(all.measures, 
       aes(x = genotypes.sampled, y = total_complexity)) + 
  geom_jitter(shape = 1, 
              color = "grey", 
              position = position_jitter(width = 0.25, height = NULL), 
              size = 1) + 
  stat_summary(fun.y = mean, geom = "point", color = "steelblue", size = 3) + 
  xlab("No. of willow genotypes") + 
  ylab(bquote('Food-web complexity ('*italic(LD[q])*')')) +
  scale_x_continuous(limits = c(0.5,25), 
                     breaks = c(1,5,10,15,20,25)) +
  scale_y_continuous(limits = c(1,2.5), #c(1,2.4)
                     breaks = seq(1, 2.5, by = 0.5)) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 9),#10
        axis.text.x = element_text(size = 9),#10
        axis.title.x = element_text(size = 11, vjust = 0.1),#vjust = 0.1, 
        axis.title.y = element_text(size = 11, vjust = 0.5),#vjust = 1, 
        panel.grid = element_blank()) 


# then plot ordination in "link_composition_plots.R" and save the figure as a pdf, portrait, 3.42" x 6"

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

## calculate the mean food web complexity for the different levels of genetic variation
web.measures.summary <- all.measures %>%
  group_by(genotypes.sampled) %>%
  summarise(mean.complexity = mean(total_complexity, na.rm = TRUE))
with(web.measures.summary, max(mean.complexity)/min(mean.complexity)) # 45% increase in food-web complexity
with(web.measures.summary, max(mean.complexity) - min(mean.complexity))
