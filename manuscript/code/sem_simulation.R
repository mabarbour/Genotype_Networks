#############################################################
#  Description: This script creates the output for the structural equation model and Figures S1 and S2 in supplementary material for the manuscript, "Genetic specificity of a plant-insect food web: Implications for linking genetic variation to food-web complexity"
#  Code author: Matthew A. Barbour
#  Email: barbour@zoology.ubc.ca
#############################################################

## Load required libraries ----
library(piecewiseSEM)
library(semPlot)
library(visreg)
library(dplyr)
library(ggplot2)

## Upload dataset ----
sem.df <- read.csv("manuscript/Dryad_data_copy/simulation_data_output/food web complexity simulation 50 reps of 200 sims 4 reps.csv") %>%
  filter(df.sim.number == 1) # note that results are qualitatively similar with all of the other data frames.

## Fig. S1 ----
# plot the relationship between genetic variation and total food web complexity (weighted linkage density)
fig.S1 <- ggplot(sem.df, 
                  aes(x = genotypes.sampled, 
                      y = total_complexity)) + 
  geom_jitter(color = "grey", shape = 1, 
              position = position_jitter(width = 0.25, 
                                         height = NULL), 
              size = 4) +
  stat_summary(fun.y = mean, geom = "point", color = "steelblue",
               shape = 20, size = 8) +
  xlab("No. of willow genotypes") + 
  ylab(bquote('Food-web complexity ('*italic(LD[q])*')')) +
  scale_x_continuous(limits = c(1,25), 
                     breaks = c(1,5,10,15,20,25)) +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),#9
        axis.title.x = element_text(size = 15, vjust = 0.1),
        axis.title.y = element_text(size = 15, vjust = 0.5),
        panel.grid = element_blank()) 
fig.S1

## Manage data for SEM analysis ----

# calculate gall richness, evenness, and total abundance
sem.df$gall.rich <- rowSums(select(sem.df, willow_vLG:willow_SG)>0)
sem.df$gall.even <- with(sem.df, vulnerability.plant_gall/gall.rich)
sem.df$gall.abund <- rowSums(select(sem.df, willow_vLG:willow_SG))

# Limit data frame to less than 5 genotypes to focus on the initial linear increase in food-web complexity
sem.df.less5 <- filter(sem.df, genotypes.sampled < 5) 

# examine correlationns among variables
corrless5 <- psych::corr.test(sem.df.less5[ ,c("genotypes.sampled","gall.even","gall.rich","gall.abund","vulnerability.gall_ptoid","generality.gall_ptoid","total_complexity")])
round(corrless5$r,2)

## Specify model ----

# create model pieces 
even.less5 <- lm(gall.even ~ genotypes.sampled, sem.df.less5)
visreg(even.less5)
summary(even.less5)

abund.less5 <- lm(gall.abund ~ genotypes.sampled, sem.df.less5)
visreg(abund.less5)
summary(abund.less5)

rich.less5 <- lm(gall.rich ~ genotypes.sampled, sem.df.less5)
visreg(rich.less5)
summary(rich.less5)

v.gp.less5 <- lm(vulnerability.gall_ptoid ~ gall.abund + gall.rich + gall.even, sem.df.less5) 
summary(v.gp.less5) 
car::vif(v.gp.less5) # variance inflation factor analysis okay

g.gp.less5 <- lm(generality.gall_ptoid ~ gall.rich + gall.even + gall.abund, sem.df.less5)
summary(g.gp.less5) 
car::vif(g.gp.less5)

c.less5 <- lm(total_complexity ~ generality.gall_ptoid + vulnerability.gall_ptoid + gall.rich + gall.even + gall.abund, sem.df.less5)
car::vif(c.less5) # okay
summary(c.less5) 

# piece models together into a list
list.less5 <- list(even.less5, rich.less5, abund.less5, v.gp.less5, g.gp.less5, c.less5)

## Goodness of fit for piecewiseSEM ----
# Note that when P > 0.05 suggest that the model provides an adequate fit to the data.
sem.fit(list.less5, sem.df.less5, 
        corr.errors = c("gall.even~~gall.rich",
                        "gall.even~~gall.abund",
                        "gall.abund~~gall.rich",
                        "vulnerability.gall_ptoid~~generality.gall_ptoid"))

## Fig. S2: SEM plot ----

# rewrite model for lavaan package
lav.less5 <- sem.lavaan(list.less5, 
                        sem.df.less5,c("gall.even~~gall.rich",
                                       "gall.even~~gall.abund",
                                       "gall.abund~~gall.rich",
                                       "vulnerability.gall_ptoid~~generality.gall_ptoid"))
lavaan::varTable(lav.less5)

# define the label that will go into the nodes
lbls<-c("Gall\nevenness","Gall\nrichness","Gall\nabundance",
        "Vulnerability\ngall-parasitoid",
        "Generality\ngall-parasitoid","Food-web\ncomplexity",
        "Genetic\nvariation")

# define the groups
grps<-list(Genotype=c("genotypes.sampled"),
           Galls=c("gall.even","gall.abund","gall.rich"),
           Gall.Ptoids = c("vulnerability.gall_ptoid","generality.gall_ptoid"),
           Food.web = c("total_complexity"))

# define the layout
ly<-matrix(c(0,0.5,
             0,-0.5,
             -0.2,0.15,
             
             0.55,0.4,
             0.5,-0.5,
             0.75, 0,
             -0.75,0),ncol=2,byrow=TRUE)

# generate plot
semPaths(lav.less5,what="std",
         #whatLabels = "omit",
         #details = TRUE,
         layout=ly,
         residuals=FALSE,
         nCharNodes=0, groups=grps,
         color=c("#009E73","#999999",#"#F0E442",#"#56B4E9",
                 "#999999","#E69F00"),
         nodeLabels=lbls, sizeMan=8, 
         posCol="blue",
         label.cex = 5,
         #edge.label.bg = FALSE,
         cut = 0.45,
         negCol = "red",
         edge.width = 2,
         asize = 1.75,
         #colFactor = 3,
         label.norm = "OOOOOOOOOOOOOOOOOOOOOOOOO",
         label.prop = 6,
         edge.label.cex=0.75,legend=FALSE)
