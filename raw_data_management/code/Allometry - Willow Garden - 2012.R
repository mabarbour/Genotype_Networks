#### Metadata
# Data Collection Date: September 29, 2012 (end of growing season)

library(dplyr)

# setworking directions
#setwd('~/Documents/Genotype_Networks/data/')

# upload data
allom <- read.csv("raw_data_management/data/Allometry - Willow Garden - 2012.csv")
stem.diams <- read.csv("raw_data_management/data/survey_2_stem_diams.csv", skip=1)

# explore data
head(allom)
str(allom)
str(stem.diams)

# look at data distributions
hist(allom$Stem.Diamter)
hist(log(allom$Stem.Diamter)) # normalizes distribution

hist(allom$Leaf.Estimate)
hist(log(allom$Leaf.Estimate)) # normalizes distribution more

hist(allom$Shoot.Count)
hist(log(allom$Shoot.Count)) # distribution more normalized except for one outlier

### fit models. Shoot allometry appears to provide the tightest relationship. data point 18 (lowest value) appears to be an outlier for both shoot and leaf allometry. Removing this value doesn't alter the slope estimate much, but it does change the intercept. I will conduct analyses with both models to see if removing this data point affects the outcome of this relationship.

# shoot allometry. Equation = log(Shoot.Count) = 3.7245*log(Stem.Diamter) - 8.7794
plot(log(Shoot.Count)~log(Stem.Diamter), allom)
plot(Shoot.Count ~ Stem.Diamter, allom)

lm.shoot.notrans <- lm(Shoot.Count ~ Stem.Diamter, allom)
summary(lm.shoot.notrans)
abline(lm.shoot.notrans)
plot(lm.shoot.notrans)

plot(log(Shoot.Count) ~ log(Stem.Diamter), allom[-18,])
lm.shoot <- lm(log(Shoot.Count)~log(Stem.Diamter), allom[-18,])
summary(lm.shoot) # model explains 61.3% of the variance with the quadratic fit.
abline(lm.shoot)
plot(lm.shoot) # point 18 appears to be an outlier, but is also the smallest branch s

library(visreg)
lm.shoot.visreg <- visreg(lm.shoot, trans = exp)


# leaf allometry.
plot(log(Leaf.Estimate)~log(Stem.Diamter),allom)

lm.leaf <- lm(log(Leaf.Estimate)~log(Stem.Diamter),allom)
summary(lm.leaf) # slightly less variance explained compared to shoot model.
abline(lm.leaf)
plot(lm.leaf)

# allometric relationships. Focusing on shoots because this relationship was stronger 

# log(Shoot.Count) = 3.7245*log(Stem.Diamter) - 8.7794  # all data points.  61.3% of variance explained
# log(Shoot.Count) = log(Stem.Diameter)*3.0371 - 6.1602  # datapoint 18 removed. 54.6% of variance explained. Although lower variance explained, it appears to be a better fit of the data by inspecting the residuals.

stem.diams.shootEst <- transform(stem.diams, shootEst.all = exp(log(Stem.Diameter.mm)*3.7245 - 8.7794), shootEst.no18 = exp(log(Stem.Diameter.mm)*3.0371 - 6.1602))

stem.diams.shootEst_combine373 <- stem.diams.shootEst %>%
  group_by(Plant.Position) %>%
  filter(Plant.Position != 656) %>%
  summarise(Stem.Diameter.mm = sum(Stem.Diameter.mm), Galls.found = mean(Galls.found), shootEst.all = sum(shootEst.all), shootEst.no18 = sum(shootEst.no18))# %>% # combines 2 branches from plant #373 together also removes #656 from dataset, where the entire tree was sampled on both surveys



hist(stem.diams.shootEst_combine373$shootEst.all)
hist(stem.diams.shootEst_combine373$shootEst.no18) # much more normally distributed and not as large of a range. I think this is currently the better estimator.

write.csv(stem.diams.shootEst_combine373, file="survey_2_stem_diams_shootEsts.csv")
