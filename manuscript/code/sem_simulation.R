
## load required libraries
require(piecewiseSEM)
require(semPlot)
require(visreg)
require(dplyr)
require(ggplot2)

## upload dataset
sem.df <- read.csv("manuscript/Dryad_data_copy/food web complexity simulation 50 reps of 200 sims 4 reps.csv") %>%
  filter(df.sim.number == 1) # note that results are qualitatively similar with all of the other data frames.


## raw data
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
  #scale_y_continuous(limits = c(1,2.25), #c(1,2.4)
   #                  breaks = seq(1, 2.25, by = 0.25)) +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),#9
        axis.title.x = element_text(size = 15, vjust = 0.1),
        axis.title.y = element_text(size = 15, vjust = 0.5),
        panel.grid = element_blank()) 
fig.S1

## calculate gall richness, evennes, and total abundance
sem.df$gall.rich <- rowSums(select(sem.df, willow_vLG:willow_SG)>0)
sem.df$gall.even <- with(sem.df, vulnerability.plant_gall/gall.rich)
sem.df$gall.abund <- rowSums(select(sem.df, willow_vLG:willow_SG))

## Limit data frame to less than 5 genotypes to focus on the initial linear increase in food-web complexity
sem.df.less5 <- filter(sem.df, genotypes.sampled < 5) 

corrless5 <- psych::corr.test(sem.df.less5[ ,c("genotypes.sampled","gall.even","gall.rich","gall.abund","vulnerability.gall_ptoid","generality.gall_ptoid","total_complexity")])
round(corrless5$r,2)

even.less5 <- lm(gall.even ~ genotypes.sampled, sem.df.less5)
#visreg(even.less5)
summary(even.less5)
#QuantPsyc::lm.beta(even.less5)

abund.less5 <- lm(gall.abund ~ genotypes.sampled, sem.df.less5)
#visreg(abund.less5)
summary(abund.less5)
#QuantPsyc::lm.beta(abund.less5)

rich.less5 <- lm(gall.rich ~ genotypes.sampled, sem.df.less5)
#visreg(rich.less5)
summary(rich.less5)
#plot(rich.less5)
#QuantPsyc::lm.beta(rich.less5)

#div.less5 <- lm(vulnerability.plant_gall ~ genotypes.sampled, sem.df.less5)
#summary(div.less5)
#QuantPsyc::lm.beta(div.less5)

v.gp.less5 <- lm(vulnerability.gall_ptoid ~ gall.abund + gall.rich + gall.even, 
                  sem.df.less5) # doesn't seem like I need
summary(v.gp.less5) #
#QuantPsyc::lm.beta(v.gp.less5)
car::vif(v.gp.less5) #okay
#visreg(v.gp.less5)
#plot(v.gp.less5) # residuals look okay

g.gp.less5 <- lm(generality.gall_ptoid ~ gall.rich + gall.even + gall.abund, sem.df.less5) # some low outliers
summary(g.gp.less5) # 
car::vif(g.gp.less5)
#QuantPsyc::lm.beta(g.gp.less5)
#visreg(g.gp.less5)
#plot(g.gp.less5) # residuals okay

c.less5 <- lm(total_complexity ~ generality.gall_ptoid + vulnerability.gall_ptoid + gall.rich + gall.even + gall.abund, sem.df.less5)
#QuantPsyc::lm.beta(c.less5)
car::vif(c.less5) # okay
summary(c.less5) # completely identified, can I use these coefficents?
#visreg(c.less5)
#plot(c.less5) # not good

list.less5 <- list(even.less5, rich.less5, abund.less5, v.gp.less5, g.gp.less5, c.less5)

# can't test goodness of fit because everything is identified
sem.fit(list.less5, sem.df.less5, 
        corr.errors = c("gall.even~~gall.rich",
                        "gall.even~~gall.abund",
                        "gall.abund~~gall.rich",
                        "vulnerability.gall_ptoid~~generality.gall_ptoid"))
coef.less5 <- sem.coefs(list.less5, sem.df.less5, 
          corr.errors = c("gall.even~~gall.rich",
                          "gall.even~~gall.abund",
                          "gall.abund~~gall.rich",
                          "vulnerability.gall_ptoid~~generality.gall_ptoid"), standardize = "scale")

# Effect of abundance. Mainly due to indirect effect via vulnerability
g.ab <- 0.685
ab.vul <- 0.623
vul.c <- 0.648
gen.c <- 0.260
ab.gen <- -0.099 # consider removing, because non-significant
ab.c <- 0.007 # similarly, consider removing, because non-significant
ab.vul.path <- g.ab*ab.vul*vul.c # abund -> vulnerability
ab.gen.path <- g.ab*ab.gen*gen.c # abund -> generality
ab.c.path <- g.ab*ab.c # abund -> complexity
ab.vul.path + ab.gen.path + ab.c.path # 0.26

# Total effect of evenness. Net negative effect.
g.ev <- -0.193
ev.vul <- -0.316
ev.gen <- 0.283
ev.c <- 0.583
ev.vul.path <- g.ev*ev.vul*vul.c # even -> vulnerability
ev.gen.path <- g.ev*ev.gen*gen.c # abund -> generality
ev.c.path <- g.ev*ev.c # abund -> complexity
ev.vul.path + ev.gen.path + ev.c.path # -0.09 

# Total effect of richness. Mainly due to direct effect on complexity
g.r <- 0.485
r.vul <- -0.169
r.gen <- 0.244
r.c <- 0.783
r.vul.path <- g.r*r.vul*vul.c # rich -> vulnerability
r.gen.path <- g.r*r.gen*gen.c # rich -> generality
r.c.path <- g.r*r.c # rich -> complexity
r.vul.path + r.gen.path + r.c.path # 0.36

lav.less5 <- sem.lavaan(list.less5, 
                        sem.df.less5,c("gall.even~~gall.rich",
                                       "gall.even~~gall.abund",
                                       "gall.abund~~gall.rich",
                                       "vulnerability.gall_ptoid~~generality.gall_ptoid"))
lavaan::varTable(lav.less5)


## Create plot of path analysis

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
#text(0.9,0.9,labels="Some text about\nthe model or\nthe weather in Indonesia")

#### Everything below this is old.
#### Other possible models

## more6 genotypes.Richness has maxed out, so I have removed it from all analyses.

sem.df.more6 <- filter(sem.df, genotypes.sampled > 4) 

corrmore6 <- psych::corr.test(sem.df.more6[ ,c("gall.even","gall.rich","gall.abund","willow_vLG","willow_rG","willow_aSG","willow_SG")])
round(corrmore6$r,2)

abund.more6 <- lm(gall.abund ~ genotypes.sampled, sem.df.more6)
visreg(abund.more6)
summary(abund.more6) # potentially too high of collinearity...
QuantPsyc::lm.beta(abund.more6)

even.more6 <- lm(gall.even ~ genotypes.sampled, sem.df.more6)
visreg(even.more6)
summary(even.more6)
plot(even.more6)
QuantPsyc::lm.beta(even.more6)

rich.more6 <- lm(gall.rich ~ genotypes.sampled, sem.df.more6)
visreg(rich.more6)
plot(rich.more6)
summary(rich.more6)
QuantPsyc::lm.beta(rich.more6)

v.gp.more6 <- lm(vulnerability.gall_ptoid ~ gall.even + gall.rich + gall.abund, 
                 sem.df.more6)
summary(v.gp.more6) #
QuantPsyc::lm.beta(v.gp.more6)
car::vif(v.gp.more6) #okay
visreg(v.gp.more6)
plot(v.gp.more6) # residuals look okay

g.gp.more6 <- lm(generality.gall_ptoid ~ gall.even + gall.rich + gall.abund, sem.df.more6) #
summary(g.gp.more6) # 
car::vif(g.gp.more6)
QuantPsyc::lm.beta(g.gp.more6)
visreg(g.gp.more6)
#plot(g.gp.more6) # residuals okay

c.more6 <- lm(total_complexity ~ generality.gall_ptoid + vulnerability.gall_ptoid + gall.even + gall.rich + gall.abund, sem.df.more6)
QuantPsyc::lm.beta(c.more6)
car::vif(c.more6) # okay
summary(c.more6) # completely identified, can I use these coefficents?
#visreg(c.more6)
plot(c.more6) # not good

list.more6 <- list(even.more6, v.gp.more6, g.gp.more6, c.more6)

# can't test goodness of fit because everything is identified
sem.fit(list.more6, sem.df.more6, corr.errors = c("gall.even~~gall.rich","vulnerability.gall_ptoid~~generality.gall_ptoid"))
sem.coefs(list.more6, sem.df.more6, corr.errors = c("gall.even~~gall.rich","vulnerability.gall_ptoid~~generality.gall_ptoid"), standardize = "scale")

lav.more6 <- sem.lavaan(list.more6, sem.df.more6,c("gall.even~~gall.rich","vulnerability.gall_ptoid~~generality.gall_ptoid"))
semPaths(lav.more6, 
         what = "std", 
         layout = "spring",
         nCharNodes = 0, 
         edge.label.cex = 1.5,
         label.cex = 2,
         label.prop = 3)

## less7_vLG genotypes. 
ggplot(sem.df, aes(x = genotypes.sampled, y = gall.even)) +
  geom_point() + geom_smooth()
ggplot(sem.df, aes(x = genotypes.sampled, y = gall.rich)) +
  geom_point() + geom_smooth()
ggplot(sem.df, aes(x = genotypes.sampled, y = gall.abund)) +
  geom_point() + geom_smooth()
ggplot(sem.df, aes(x = genotypes.sampled, y = willow_vLG/gall.abund)) +
  geom_point() + geom_smooth()

sem.df.less7_vLG <- filter(sem.df, genotypes.sampled < 5) 

corrless7_vLG <- psych::corr.test(sem.df.less7_vLG[ ,c("gall.even","gall.rich","gall.abund","willow_vLG","willow_rG","willow_aSG","willow_SG")])
round(corrless7_vLG$r,2)

abund.less7_vLG <- lm(gall.abund ~ genotypes.sampled, sem.df.less7_vLG)
visreg(abund.less7_vLG)
summary(abund.less7_vLG)

rich.less7_vLG <- lm(gall.rich ~ genotypes.sampled, sem.df.less7_vLG)
visreg(rich.less7_vLG)
summary(rich.less7_vLG)

#even.less7_vLG <- lm(gall.even ~ genotypes.sampled, sem.df.less7_vLG)
visreg(even.less7_vLG)
summary(even.less7_vLG)

v.gp.less7_vLG <- lm(vulnerability.gall_ptoid ~ gall.abund + gall.rich, 
                 sem.df.less7_vLG)
summary(v.gp.less7_vLG) #
QuantPsyc::lm.beta(v.gp.less7_vLG)
car::vif(v.gp.less7_vLG) #okay
visreg(v.gp.less7_vLG)
plot(v.gp.less7_vLG) # residuals look okay

g.gp.less7_vLG <- lm(generality.gall_ptoid ~ gall.abund + gall.rich, sem.df.less7_vLG) #
summary(g.gp.less7_vLG) # 
car::vif(g.gp.less7_vLG)
QuantPsyc::lm.beta(g.gp.less7_vLG)
visreg(g.gp.less7_vLG)
plot(g.gp.less7_vLG) # residuals okay

c.less7_vLG <- lm(total_complexity ~ generality.gall_ptoid + vulnerability.gall_ptoid + gall.abund + gall.rich, sem.df.less7_vLG)
QuantPsyc::lm.beta(c.less7_vLG)
car::vif(c.less7_vLG) # okay
summary(c.less7_vLG) # completely identified, can I use these coefficents?
#visreg(c.less7_vLG)
plot(c.less7_vLG) # not good

list.less7_vLG <- list(abund.less7_vLG, rich.less7_vLG, v.gp.less7_vLG, g.gp.less7_vLG, c.less7_vLG)

# can't test goodness of fit because everything is identified
# missing pathway of effect on total complexity, not made up for by direct effect of gall.abund
sem.fit(list.less7_vLG, sem.df.less7_vLG, 
        corr.errors = c("gall.abund~~gall.rich",
                        "vulnerability.gall_ptoid~~generality.gall_ptoid"))
sem.coefs(list.less7_vLG, sem.df.less7_vLG, 
          corr.errors = c("gall.abund~~gall.rich",
                          "vulnerability.gall_ptoid~~generality.gall_ptoid"), standardize = "scale")
0.59*0.72*0.57 # abund -> vulnerability
0.59*-0.25*0.35 # abund -> generality
0.59*-0.09 # abund -> complexity

0.49*0.16*0.35 # rich -> generality pathway
0.49*0.56 # rich direct pathway
0.49*-0.03*0.57 # rich -> vulnerability pathway

-0.19*-0.32*0.65 # even -> vulnerability
-0.19*0.28*0.26 # even -> generality
-0.19*0.58 # even direct

lav.less7_vLG <- sem.lavaan(list.less7_vLG, sem.df.less7_vLG,c("gall.even~~gall.rich","vulnerability.gall_ptoid~~generality.gall_ptoid"))
semPaths(lav.less7_vLG, 
         what = "std", 
         layout = "spring",
         nCharNodes = 0, 
         edge.label.cex = 1.5,
         label.cex = 2,
         label.prop = 3)

## 1 genotype 
sem.df.1 <- filter(sem.df, genotypes.sampled == 1)

corr1 <- psych::corr.test(sem.df.1[ ,c("gall.even","gall.rich","gall.abund","willow_vLG","willow_rG","willow_aSG","willow_SG")])
round(corr1$r,2)

v.gp.1 <- lm(vulnerability.gall_ptoid ~ gall.even + gall.abund + gall.rich, 
             sem.df.1)
summary(v.gp.1) # no richness effect
QuantPsyc::lm.beta(v.gp.1)
car::vif(v.gp.1) #okay
visreg(v.gp.1)
plot(v.gp.1) # residuals look okay

g.gp.1 <- lm(generality.gall_ptoid ~ gall.rich + gall.even + gall.abund, sem.df.1)
summary(g.gp.1) # no effects
car::vif(g.gp.1)
QuantPsyc::lm.beta(g.gp.1)
visreg(g.gp.1)
plot(g.gp.1) # residuals okay

c.1 <- lm(total_complexity ~ generality.gall_ptoid + vulnerability.gall_ptoid + gall.even + gall.rich, sem.df.1)
QuantPsyc::lm.beta(c.1)
car::vif(c.1) # okay
summary(c.1)
visreg(c.1)
plot(c.1) # okay

list.1 <- list(v.gp.1, g.gp.1, c.1)

# provides a good fit to the data
sem.fit(list.1, sem.df.1, corr.errors = c("gall.rich~~gall.even",
                                          "gall.rich~~gall.abund",
                                          "gall.even~~gall.abund",
                                          "vulnerability.gall_ptoid~~generality.gall_ptoid"))
sem.coefs(list.1, sem.df.1, corr.errors = c("gall.rich~~gall.even",
                                            "gall.rich~~gall.abund",
                                            "gall.even~~gall.abund",
                                            "vulnerability.gall_ptoid~~generality.gall_ptoid"), standardize = "scale")

lav.1 <- sem.lavaan(list.1, sem.df.1,c("gall.rich~~gall.even",
                                       "gall.rich~~gall.abund",
                                       "gall.even~~gall.abund",
                                       "vulnerability.gall_ptoid~~generality.gall_ptoid"))
semPaths(lav.1, 
         what = "std", 
         layout = "spring",
         nCharNodes = 0, 
         edge.label.cex = 1.5,
         label.cex = 2,
         label.prop = 3)

## 2 genotypes
sem.df.2 <- filter(sem.df, genotypes.sampled == 2)

corr2 <- psych::corr.test(sem.df.2[ ,c("gall.even","gall.rich","gall.abund","willow_vLG","willow_rG","willow_aSG","willow_SG")])
round(corr2$r,2)

v.gp.2 <- lm(vulnerability.gall_ptoid ~ gall.even + gall.abund + gall.rich, 
             sem.df.2)
summary(v.gp.2) # no richness effect
QuantPsyc::lm.beta(v.gp.2)
car::vif(v.gp.2) #okay
visreg(v.gp.2)
plot(v.gp.2) # residuals look okay

g.gp.2 <- lm(generality.gall_ptoid ~ gall.rich + gall.even + gall.abund, sem.df.2) # outlier at pt 37
summary(g.gp.2) # no effects
car::vif(g.gp.2)
QuantPsyc::lm.beta(g.gp.2)
visreg(g.gp.2)
plot(g.gp.2) # residuals okay

c.2 <- lm(total_complexity ~ generality.gall_ptoid + vulnerability.gall_ptoid + gall.even + gall.rich, sem.df.2)
QuantPsyc::lm.beta(c.2)
car::vif(c.2) # okay
summary(c.2)
#visreg(c.2)
plot(c.2) # weird

list.2 <- list(v.gp.2, g.gp.2, c.2)

# provides a good fit to the data
sem.fit(list.2, sem.df.2, corr.errors = c("gall.rich~~gall.even",
                                          "gall.rich~~gall.abund",
                                          "gall.even~~gall.abund",
                                          "vulnerability.gall_ptoid~~generality.gall_ptoid"))
sem.coefs(list.2, sem.df.2, corr.errors = c("gall.rich~~gall.even",
                                            "gall.rich~~gall.abund",
                                            "gall.even~~gall.abund",
                                            "vulnerability.gall_ptoid~~generality.gall_ptoid"), standardize = "scale")

lav.2 <- sem.lavaan(list.2, sem.df.2,c("gall.rich~~gall.even",
                                       "gall.rich~~gall.abund",
                                       "gall.even~~gall.abund",
                                       "vulnerability.gall_ptoid~~generality.gall_ptoid"))
semPaths(lav.2, 
         what = "std", 
         layout = "spring",
         nCharNodes = 0, 
         edge.label.cex = 1.5,
         label.cex = 2,
         label.prop = 3)

## 5 genotypes. 
sem.df.5 <- filter(sem.df, genotypes.sampled == 5)

corr5 <- psych::corr.test(sem.df.5[ ,c("gall.even","gall.rich","gall.abund","willow_vLG","willow_rG","willow_aSG","willow_SG")])
round(corr5$r,2)

v.gp.5 <- lm(vulnerability.gall_ptoid ~ gall.even + gall.abund + gall.rich, 
             sem.df.5)
summary(v.gp.5) #
QuantPsyc::lm.beta(v.gp.5)
car::vif(v.gp.5) #okay
visreg(v.gp.5)
plot(v.gp.5) # residuals look okay

g.gp.5 <- lm(generality.gall_ptoid ~ gall.even + gall.abund + gall.rich, sem.df.5) # outlier at pt 37
summary(g.gp.5) # no effects
car::vif(g.gp.5)
QuantPsyc::lm.beta(g.gp.5)
visreg(g.gp.5)
plot(g.gp.5) # residuals okay

c.5 <- lm(total_complexity ~ generality.gall_ptoid + vulnerability.gall_ptoid + gall.even + gall.rich, sem.df.5)
QuantPsyc::lm.beta(c.5)
car::vif(c.5) # okay
summary(c.5)
QuantPsyc::lm.beta(c.5)
#visreg(c.5)
plot(c.5) # weird

list.5 <- list(v.gp.5, g.gp.5, c.5)

# provides a good fit to the data
sem.fit(list.5, sem.df.5, corr.errors = c("gall.rich~~gall.even",
                                          "gall.rich~~gall.abund",
                                          "gall.even~~gall.abund",
                                          "vulnerability.gall_ptoid~~generality.gall_ptoid"))
sem.coefs(list.5, sem.df.5, corr.errors = c("gall.rich~~gall.even",
                                            "gall.rich~~gall.abund",
                                            "gall.even~~gall.abund",
                                            "vulnerability.gall_ptoid~~generality.gall_ptoid"), standardize = "scale")

lav.5 <- sem.lavaan(list.5, sem.df.5,c("gall.rich~~gall.even",
                                       "gall.rich~~gall.abund",
                                       "gall.even~~gall.abund",
                                       "vulnerability.gall_ptoid~~generality.gall_ptoid"))
semPaths(lav.5, 
         what = "std", 
         layout = "spring",
         nCharNodes = 0, 
         edge.label.cex = 1.5,
         label.cex = 2,
         label.prop = 3)

## 10 genotypes. Richness has maxed out, so I have removed it from all analyses.
sem.df.10 <- filter(sem.df, genotypes.sampled == 10)

corr10 <- psych::corr.test(sem.df.10[ ,c("gall.even","gall.abund","willow_vLG","willow_rG","willow_aSG","willow_SG")])
round(corr10$r,2)

v.gp.10 <- lm(vulnerability.gall_ptoid ~ gall.even + gall.abund, 
              sem.df.10)
summary(v.gp.10) #
QuantPsyc::lm.beta(v.gp.10)
car::vif(v.gp.10) #okay
visreg(v.gp.10)
plot(v.gp.10) # residuals look okay

g.gp.10 <- lm(generality.gall_ptoid ~ gall.even + gall.abund, sem.df.10) # outlier at pt 37
summary(g.gp.10) # no effects
car::vif(g.gp.10)
QuantPsyc::lm.beta(g.gp.10)
visreg(g.gp.10)
plot(g.gp.10) # residuals okay

c.10 <- lm(total_complexity ~ generality.gall_ptoid + vulnerability.gall_ptoid + gall.even, sem.df.10)
QuantPsyc::lm.beta(c.10)
car::vif(c.10) # okay
summary(c.10) # completely identified, can I use these coefficents?
QuantPsyc::lm.beta(c.10)
#visreg(c.10)
plot(c.10) # good

list.10 <- list(v.gp.10, g.gp.10, c.10)

# provides a good fit to the data
sem.fit(list.10, sem.df.10, corr.errors = c("gall.even~~gall.abund", "vulnerability.gall_ptoid~~generality.gall_ptoid"))
sem.coefs(list.10, sem.df.10, corr.errors = c("gall.even~~gall.abund",
                                              "vulnerability.gall_ptoid~~generality.gall_ptoid"), standardize = "scale")

lav.10 <- sem.lavaan(list.10, sem.df.10,c("gall.even~~gall.abund",
                                          "vulnerability.gall_ptoid~~generality.gall_ptoid"))
semPaths(lav.10, 
         what = "std", 
         layout = "spring",
         nCharNodes = 0, 
         edge.label.cex = 1.5,
         label.cex = 2,
         label.prop = 3)

## 15 genotypes. Richness has maxed out, so I have removed it from all analyses.
sem.df.15 <- filter(sem.df, genotypes.sampled == 15)

corr15 <- psych::corr.test(sem.df.15[ ,c("gall.even","gall.abund","willow_vLG","willow_rG","willow_aSG","willow_SG")])
round(corr15$r,2)

v.gp.15 <- lm(vulnerability.gall_ptoid ~ gall.even + gall.abund, 
              sem.df.15)
summary(v.gp.15) #
QuantPsyc::lm.beta(v.gp.15)
car::vif(v.gp.15) #okay
visreg(v.gp.15)
plot(v.gp.15) # residuals look okay

g.gp.15 <- lm(generality.gall_ptoid ~ gall.even + gall.abund, sem.df.15) # outlier at pt 37
summary(g.gp.15) # no effects
car::vif(g.gp.15)
QuantPsyc::lm.beta(g.gp.15)
visreg(g.gp.15)
plot(g.gp.15) # residuals okay

c.15 <- lm(total_complexity ~ generality.gall_ptoid + vulnerability.gall_ptoid + gall.even, sem.df.15)
QuantPsyc::lm.beta(c.15)
car::vif(c.15) # okay
summary(c.15) # completely identified, can I use these coefficents?
QuantPsyc::lm.beta(c.15)
#visreg(c.15)
plot(c.15) # good

list.15 <- list(v.gp.15, g.gp.15, c.15)

# provides a good fit to the data
sem.fit(list.15, sem.df.15, corr.errors = c("gall.even~~gall.abund", "vulnerability.gall_ptoid~~generality.gall_ptoid"))
sem.coefs(list.15, sem.df.15, corr.errors = c("gall.even~~gall.abund",
                                              "vulnerability.gall_ptoid~~generality.gall_ptoid"), standardize = "scale")

lav.15 <- sem.lavaan(list.15, sem.df.15,c("gall.even~~gall.abund",
                                          "vulnerability.gall_ptoid~~generality.gall_ptoid"))
semPaths(lav.15, 
         what = "std", 
         layout = "spring",
         nCharNodes = 0, 
         edge.label.cex = 1.5,
         label.cex = 2,
         label.prop = 3)

## 20 genotypes. Richness has maxed out, so I have removed it from all analyses.
sem.df.20 <- filter(sem.df, genotypes.sampled == 20)

corr20 <- psych::corr.test(sem.df.20[ ,c("gall.even","gall.abund","willow_vLG","willow_rG","willow_aSG","willow_SG")])
round(corr20$r,2)

v.gp.20 <- lm(vulnerability.gall_ptoid ~ gall.even + gall.abund, 
              sem.df.20)
summary(v.gp.20) #
QuantPsyc::lm.beta(v.gp.20)
car::vif(v.gp.20) #okay
visreg(v.gp.20)
plot(v.gp.20) # residuals look okay

g.gp.20 <- lm(generality.gall_ptoid ~ gall.even + gall.abund, sem.df.20) # outlier at pt 37
summary(g.gp.20) # no effects
car::vif(g.gp.20)
QuantPsyc::lm.beta(g.gp.20)
visreg(g.gp.20)
plot(g.gp.20) # residuals okay

c.20 <- lm(total_complexity ~ generality.gall_ptoid + vulnerability.gall_ptoid + gall.even, sem.df.20)
QuantPsyc::lm.beta(c.20)
car::vif(c.20) # okay
summary(c.20) # completely identified, can I use these coefficents?
QuantPsyc::lm.beta(c.20)
#visreg(c.20)
plot(c.20) # good

list.20 <- list(v.gp.20, g.gp.20, c.20)

# provides a good fit to the data
sem.fit(list.20, sem.df.20, corr.errors = c("gall.even~~gall.abund", "vulnerability.gall_ptoid~~generality.gall_ptoid"))
sem.coefs(list.20, sem.df.20, corr.errors = c("gall.even~~gall.abund",
                                              "vulnerability.gall_ptoid~~generality.gall_ptoid"), standardize = "scale")

lav.20 <- sem.lavaan(list.20, sem.df.20,c("gall.even~~gall.abund",
                                          "vulnerability.gall_ptoid~~generality.gall_ptoid"))
semPaths(lav.20, 
         what = "std", 
         layout = "spring",
         nCharNodes = 0, 
         edge.label.cex = 1.5,
         label.cex = 2,
         label.prop = 3)

## 20_vLG genotypes. Alternative model, replacing with just willow_vLG. Richness has maxed out, so I have removed it from all analyses.
sem.df.20_vLG <- filter(sem.df, genotypes.sampled == 20)

corr20_vLG <- psych::corr.test(sem.df.20_vLG[ ,c("gall.even","gall.abund","willow_vLG","willow_rG","willow_aSG","willow_SG")])
round(corr20_vLG$r,2)

v.gp.20_vLG <- lm(vulnerability.gall_ptoid ~ willow_vLG, 
                  sem.df.20_vLG)
summary(v.gp.20_vLG) #
QuantPsyc::lm.beta(v.gp.20_vLG)
car::vif(v.gp.20_vLG) #okay
visreg(v.gp.20_vLG)
plot(v.gp.20_vLG) # residuals look okay

g.gp.20_vLG <- lm(generality.gall_ptoid ~ willow_vLG, sem.df.20_vLG) # outlier at pt 37
summary(g.gp.20_vLG) # no effects
car::vif(g.gp.20_vLG)
QuantPsyc::lm.beta(g.gp.20_vLG)
visreg(g.gp.20_vLG)
plot(g.gp.20_vLG) # residuals okay

c.20_vLG <- lm(total_complexity ~ generality.gall_ptoid + vulnerability.gall_ptoid + willow_vLG, sem.df.20_vLG)
QuantPsyc::lm.beta(c.20_vLG)
car::vif(c.20_vLG) # okay
summary(c.20_vLG) # completely identified, can I use these coefficents?
QuantPsyc::lm.beta(c.20_vLG)
#visreg(c.20_vLG)
plot(c.20_vLG) # good

list.20_vLG <- list(v.gp.20_vLG, g.gp.20_vLG, c.20_vLG)

# can't test goodness of fit because everything is identified
sem.fit(list.20_vLG, sem.df.20_vLG, corr.errors = c("vulnerability.gall_ptoid~~generality.gall_ptoid"))
sem.coefs(list.20_vLG, sem.df.20_vLG, corr.errors = c("vulnerability.gall_ptoid~~generality.gall_ptoid"), standardize = "scale")

lav.20_vLG <- sem.lavaan(list.20_vLG, sem.df.20_vLG,c("vulnerability.gall_ptoid~~generality.gall_ptoid"))
semPaths(lav.20_vLG, 
         what = "std", 
         layout = "spring",
         nCharNodes = 0, 
         edge.label.cex = 1.5,
         label.cex = 2,
         label.prop = 3)

