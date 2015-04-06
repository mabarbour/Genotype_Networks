## load required libraries
require(visreg)
require(ggplot2)
require(grid)
require(dplyr)
require(gridBase)
#require(png)

## upload results from simulation
all.measures <- read.csv("~/Documents/Genotype_Networks/data/food web complexity simulation.csv")
dim(all.measures)[1] # 2221 unique simulations. 2425 simulations originally run

single.complexity.max <- max(filter(all.measures, Number.of.Genotypes == 1)$total_complexity, na.rm = TRUE)

## upload png of food web composition ordination
#ord.png <- readPNG("~/Documents/Genotype_Networks/figures/full.link.composition.png")
#ord.rast <- rasterGrob(ord.png, interpolate = TRUE)

# setup multipanel plot for combining ggplot and base graphics (ordination from "link_composition_plots")
par(mfrow = c(1,2))
plot.new()
vps <- baseViewports()
pushViewport(vps$figure)
vp1 <- plotViewport(margins = c(2.4, 0.75, 3.1, 0) - 1)#c(1.8, 1, 0, 1))

A.total <- data.frame(x = 0.5, y = 2.4, labels = "(A)")

## plot the relationship between genetic variation and total food web complexity (weighted linkage density)
total <- ggplot(all.measures, 
       aes(x = Number.of.Genotypes, y = total_complexity)) + 
  geom_hline(yintercept = single.complexity.max, linetype = "dashed") +
  geom_jitter(shape = 1, color = "grey", position = position_jitter(width = 0.15, height = NULL)) + 
  stat_summary(fun.y = mean, geom = "point", color = "steelblue", size = 5) + 
  xlab("Genotypic diversity (no. of genotypes)") + ylab("Food web complexity index") +
  scale_x_continuous(limits = c(0,25), breaks = 1:25) +
  scale_y_continuous(limits = c(1,2.4), breaks = seq(1,2.5, by = 0.25)) +
  geom_text(data = A.total, aes(x = x, y = y, label = labels), 
            inherit.aes = FALSE, size = 6) +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 16, vjust = 0.1),
        axis.title.y = element_text(size = 16, vjust = 1),
        panel.grid = element_blank()) 

print(total, vp = vp1)
# then plot ordination in "link_composition_plots.R" and save the figure as a pdf, landscape, 8.5" x 11"

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

## linear model with just the mean data. 
web.sum.lm <- lm(log(mean.complexity) ~ log(Number.of.Genotypes), data = web.measures.summary)
summary(web.sum.lm) # Need to understand interpretation here.
visreg(web.sum.lm) # not quite there

## linear model with all of the data. 
web.lm <- lm(log(total_complexity) ~ log(Number.of.Genotypes), data = all.measures)
summary(web.lm)
visreg(web.lm)
