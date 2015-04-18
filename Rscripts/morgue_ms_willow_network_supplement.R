```{r}
require(visreg)
abund <- which(full.df$vLG_abund > 0)
summary(full.df$vLG_abund[abund])
cut.vLG_abund <- cut(x = full.df$vLG_abund, breaks = c(1, 4, 22))
ptism.glm <- glm(vLG_parasitized/vLG_abund ~ vLG.height.mean+vLG_abund, data = full.df, 
                 weights = vLG_abund, family = "binomial")
summary(ptism.glm)
test <- visreg(ptism.glm, xvar = "vLG.height.mean", by = "vLG_abund", scale = "response")
```

```{r}
## SG parasitism
SG.ptized.glm <- glm(SG_parasitized/SG_abund ~ Genotype, data = full.df, weights = SG_abund, family = "binomial")
summary(SG.ptized.glm)
anova(SG.ptized.glm, test = "LR")

## vLG parasitism
vLG.ptized.df <- full.predictors %>%
  mutate(vLG_parasitized = vLG_Tory + vLG_Eulo + vLG_Platy + vLG_Mesopol + vLG_Mymarid)
table(vLG.ptized.df$Genotype) # data is unbalanced, but this doesn't qualitatively affect the outcome of the glm.
vLG.ptized.3plus.genos <- c("*","B","D","E","F","H","I","K","L","O","S","V","X","Y","Z")
colSums(vLG.ptized.df[ ,c("vLG_parasitized","vLG_abund","vLG_Tory","vLG_Eulo","vLG_Platy",
                          "vLG_Mesopol","vLG_Mymarid")])

# interesting, only marginally significant with quasibinomial
vLG.ptized.bin <- glm(cbind(vLG_parasitized, vLG_abund - vLG_parasitized) ~ Genotype,
                      data = filter(vLG.ptized.df, Genotype %in% vLG.ptized.3plus.genos), 
                      family = 'binomial')
summary(vLG.ptized.bin)
dfun(vLG.ptized.bin)
anova(vLG.ptized.bin, test = "LR")
require(dispmod)
anova(glm.binomial.disp(vLG.ptized.bin), test = "LR") # still significant after accounting for overdispersion
visreg(vLG.ptized.bin, scale = "response")

vLG.ptized.pred <- glm(cbind(vLG_parasitized, vLG_abund - vLG_parasitized) ~ vLG.height.mean,
                       data = vLG.ptized.df,
                       family = 'binomial')
summary(vLG.ptized.pred)
1-exp(coef(vLG.ptized.pred)[2]) # 25% reduction in the oddes of vLG being parasitized with every one unit increase in gall size.
dfun(vLG.ptized.pred) # overdispersion
anova(glm.binomial.disp(vLG.ptized.pred), test = "LR")
require(visreg)
visreg(vLG.ptized.pred, scale = "response")



# vLG_Platy test
vLG_Platy.bin <- glm(cbind(vLG_Platy, vLG_abund - vLG_Platy) ~ Genotype,
                     data = filter(full.df, Genotype %in% vLG.ptized.3plus.genos), 
                     family = 'quasibinomial')
plot(vLG_Platy.bin)
summary(vLG_Platy.bin)
anova(vLG_Platy.bin, test = "F")

library(dispmod)
vLG_Platy.bin.pred <- glm(cbind(vLG_Platy, vLG_abund - vLG_Platy) ~ vLG.height.mean,
                          data = full.predictors, 
                          family = 'binomial')
dfun(vLG_Platy.bin.pred)
summary(vLG_Platy.bin.pred)
glm.binomial.disp(vLG_Platy.bin.pred)
AIC(vLG_Platy.bin.pred)
plot(vLG_Platy.bin.pred)


# vLG_Mesopol test
vLG_Mesopol.bin <- glm(cbind(vLG_Mesopol, vLG_abund - vLG_Mesopol) ~ Genotype,
                       data = filter(full.df, Genotype %in% vLG.ptized.3plus.genos), 
                       family = 'binomial')
plot(vLG_Mesopol.bin)
summary(vLG_Mesopol.bin)
dfun(vLG_Mesopol.bin)
anova(vLG_Mesopol.bin, test = "LR")


# vLG_Tory test
vLG_Tory.bin <- glm(cbind(vLG_Tory, vLG_abund - vLG_Tory) ~ Genotype,
                    data = filter(full.df, Genotype %in% vLG.ptized.3plus.genos), 
                    family = 'binomial')
plot(vLG_Tory.bin)
summary(vLG_Tory.bin)
anova(vLG_Tory.bin, test = "LR")

vLG_Tory.bin.pred <- glm(cbind(vLG_Tory, vLG_abund - vLG_Tory) ~ vLG_abund + vLG.height.mean,
                         data = full.df, 
                         family = 'binomial')
summary(vLG_Tory.bin.pred)
AIC(vLG_Tory.bin.pred) # 136.1538 is the lowest score. Main effect model
visreg(vLG_Tory.bin.pred, scale = "response")

# vLG_Eulo test
vLG_Eulo.bin <- glm(cbind(vLG_Eulo, vLG_abund - vLG_Eulo) ~ Genotype,
                    data = filter(full.df, Genotype %in% vLG.ptized.3plus.genos), 
                    family = 'binomial')
#plot(vLG_Eulo.bin)
summary(vLG_Eulo.bin)
anova(vLG_Eulo.bin, test = "LR") # not significant

# not testing vLG_Mymarid because it is in such low abundance
predict.df <- data.frame(Genotype = vLG.ptized.3plus.genos)
vLG.genotype.predictions <- cbind.data.frame(predict.df, 
                                             vLG.Proportion.Ptized = predict(vLG.ptized.bin,
                                                                             newdata = predict.df,
                                                                             type = "response"))
summary(vLG.genotype.predictions)

ggplot(filter(vLG.ptized.df, Genotype %in% vLG.ptized.3plus.genos),
       aes(x = Genotype, y = vLG_parasitized/vLG_abund)) +
  geom_point(aes(size = vLG_abund), color = "grey", shape = 1,
             position = position_jitter(height = 0, width = 0.1)) + 
  scale_size(range = c(3,15)) + 
  geom_point(data = vLG.genotype.predictions, 
             aes(x = Genotype, y = vLG.Proportion.Ptized),
             color = "steelblue", size = 30, shape = "--") +
  theme_classic() +
  ylab("Proportion of galls parasitized")
range(filter(vLG.ptized.df, Genotype %in% vLG.ptized.3plus.genos)$vLG_abund)

## rG parasitism
rG.ptized.df <- full.df %>%
  filter(rG_abund > 0)
table(rG.ptized.df$Genotype)
rG.ptized.3plus.genos <- c("*","A","D", "F", "G", "H","I","K","O","Q", "R", "S","X", "Y", "Z") 
colSums(rG.ptized.df[ ,c("rG_parasitized","rG_abund","rG_Tory","rG_Eulo","rG_Platy",
                         "rG_Mesopol","rG_Lestodip")]) # too few rG_eulo, rG_Platy, rG_Mesopol, and rG_Lestodip to test for differences in parasitism among genotypes.

rG.ptized.bin <- glm(cbind(rG_parasitized, rG_abund - rG_parasitized) ~ Genotype,
                     data = filter(full.df, Genotype %in% rG.ptized.3plus.genos), 
                     family = 'binomial')
#plot(rG.ptized.bin) # residuals not great
summary(rG.ptized.bin)
anova(rG.ptized.bin, test = "LR") # marginally significant after retaining Genotypes with 3 plus replicates.

# rG_Tory test
rG_Tory.bin <- glm(cbind(rG_Tory, rG_abund - rG_Tory) ~ Genotype,
                   data = filter(full.df, Genotype %in% rG.ptized.3plus.genos), 
                   family = 'binomial')
plot(rG_Tory.bin)
summary(rG_Tory.bin)
anova(rG_Tory.bin, test = "LR") # not significant

rG_Eulo.bin <- glm(cbind(rG_Eulo, rG_abund - rG_Eulo) ~ Genotype,
                   data = filter(full.df, Genotype %in% rG.ptized.3plus.genos), 
                   family = 'binomial')
plot(rG_Eulo.bin)
summary(rG_Eulo.bin)
anova(rG_Eulo.bin, test = "LR") # not significant

predict.df <- data.frame(Genotype = rG.ptized.3plus.genos)
rG.genotype.predictions <- cbind.data.frame(predict.df, 
                                            rG.Proportion.Ptized = predict(rG.ptized.bin,
                                                                           newdata = predict.df,
                                                                           type = "response"))
summary(rG.genotype.predictions)

ggplot(filter(rG.ptized.df, Genotype %in% rG.ptized.3plus.genos),
       aes(x = Genotype, y = rG_parasitized/rG_abund)) +
  geom_point(aes(size = rG_abund), color = "grey", shape = 1,
             position = position_jitter(height = 0, width = 0.1)) + 
  scale_size(range = c(3,15)) + 
  geom_point(data = rG.genotype.predictions, 
             aes(x = Genotype, y = rG.Proportion.Ptized),
             color = "steelblue", size = 30, shape = "--") +
  theme_classic() +
  ylab("Proportion of galls parasitized")
range(filter(rG.ptized.df, Genotype %in% rG.ptized.3plus.genos)$rG_abund)

## aSG parasitism
aSG.ptized.df <- full.df %>%
  filter(aSG_abund > 0)
table(aSG.ptized.df$Genotype)
aSG.ptized.3plus.genos <- c("G","I","K","Q") 

aSG.ptized.bin <- glm(cbind(aSG_parasitized, aSG_abund - aSG_parasitized) ~ Genotype,
                      data = filter(full.df, Genotype %in% aSG.ptized.3plus.genos), 
                      family = 'binomial')
plot(aSG.ptized.bin) # residuals not great
summary(aSG.ptized.bin)
anova(aSG.ptized.bin, test = "LR") # marginally significant after retaining Genotypes with 3 plus replicates.
predict.df <- data.frame(Genotype = aSG.ptized.3plus.genos)
aSG.genotype.predictions <- cbind.data.frame(predict.df, 
                                             aSG.Proportion.Ptized = predict(aSG.ptized.bin,
                                                                             newdata = predict.df,
                                                                             type = "response"))
summary(aSG.genotype.predictions)

ggplot(filter(aSG.ptized.df, Genotype %in% aSG.ptized.3plus.genos),
       aes(x = Genotype, y = aSG_parasitized/aSG_abund)) +
  geom_point(aes(size = aSG_abund), color = "grey", shape = 1,
             position = position_jitter(height = 0, width = 0.1)) + 
  scale_size(range = c(3,15)) + 
  geom_point(data = aSG.genotype.predictions, 
             aes(x = Genotype, y = aSG.Proportion.Ptized),
             color = "steelblue", size = 30, shape = "--") +
  theme_classic() +
  ylab("Proportion of galls parasitized")
range(filter(aSG.ptized.df, Genotype %in% aSG.ptized.3plus.genos)$aSG_abund)

## SG parasitism
SG.ptized.df <- full.df %>%
  filter(SG_abund > 0)
table(SG.ptized.df$Genotype) # no Genotypes with 3+ replicates, so I'm unable to evaluate whether SG parasitism defers among willow genotypes. 

##
uni.check.bin <- glm(vLG_Platy/vLG_abund ~ vLG_abund*vLG.height.mean, full.predictors,
                     family = binomial,
                     weights = vLG_abund)
AIC(uni.check.bin)
summary(uni.check.bin)
visreg(uni.check.bin, xvar = "vLG_abund", by = "vLG.height.mean", scale = "response")


# NEED TO THINK ABOUT HOW I SCALE THIS UP TO DIFFERENT NUMBERS OF SHOOTS SAMPLED.

#### messing around with different ideas.
## Create mvabund model based on probability of observing an interaction
prob.df <- full.predictors %>%
  mutate(succ.vLG_Platy.succ = vLG_Platy,
         fail.vLG_Platy.fail = vLG_abund - vLG_Platy,
         succ.vLG_Tory.succ = vLG_Tory,
         fail.vLG_Tory.fail = vLG_abund - vLG_Tory,
         succ.vLG_Mesopol.succ = vLG_Mesopol,
         fail.vLG_Mesopol.fail = vLG_abund - vLG_Mesopol,
         succ.vLG_Eulo.succ = vLG_Eulo,
         fail.vLG_Eulo.fail = vLG_abund - vLG_Eulo,
         succ.rG_Platy.succ = rG_Platy,
         fail.rG_Platy.fail = rG_abund - rG_Platy,
         succ.rG_Tory.succ = rG_Tory,
         fail.rG_Tory.fail = rG_abund - rG_Tory,
         succ.rG_Mesopol.succ = rG_Mesopol,
         fail.rG_Mesopol.fail = rG_abund - rG_Mesopol,
         succ.rG_Eulo.succ = rG_Eulo,
         fail.rG_Eulo.fail = rG_abund - rG_Eulo,
         succ.SG_Platy.succ = SG_Platy,
         fail.SG_Platy.fail = SG_abund - SG_Platy,
         succ.aSG_Tory.succ = aSG_Tory,
         fail.aSG_Tory.fail = aSG_abund - aSG_Tory) %>%
  select(succ.vLG_Platy.succ:fail.aSG_Tory.fail) %>%
  mvabund()

vLG.interaction.df <- mvabund(select(full.predictors, vLG_Eulo, vLG_Platy, vLG_Mymarid, vLG_Tory, vLG_Mesopol))
vLG.interaction.df[vLG.interaction.df > 0] <- 1

test <- manyglm(vLG.interaction.df ~ vLG_density*vLG.height.mean, data = full.predictors, family = "binomial")
plot(test)
anova(test, p.uni = "adjusted")
length(prob.df[prob.df > 1] > 0)

test <- MASS::glm.nb(vLG_Platy ~ offset(log(vLG_abund)) + vLG.height.mean, data = full.predictors) #
summary(test)
coef(test)
AIC(test)
#exp(sum(-4.5 + .34*5 -0.25*8))
exp(sum(coef(test)[c(1:3)]*c(1,0.02,8))) # 8 mm gall at a density of 11 per 100 shoots. Would be 6.11 parasitoids per 100 shoots.
predict(test, type = "response")*10
visreg(test, scale = "response")
```

EVERYTHING BELOW THIS APPEARS TO BE OLD AND MAY NOT BE USEFUL.

```{r}
row.count <- dim(full.predictors)[1]
arrange(full.predictors, vLG.height.mean)$vLG.height.mean

## generate a data frame that varies on Iteomyia abundance but holds all other variables constant at their means.
predict.log.vLG_abund <- with(full.predictors,
                              data.frame(log.vLG_abund = seq(from = min(log.vLG_abund),
                                                             to = max(log.vLG_abund),
                                                             length.out = row.count),
                                         vLG.height.mean = rep(mean(vLG.height.mean), row.count),
                                         log.1.rG_abund = rep(mean(log.1.rG_abund), row.count),
                                         log.1.aSG_abund = rep(mean(log.1.aSG_abund), row.count)))

## predict link abundance based on the Iteomyia abundance variation.
log.vLG_abund.predict <- predict(net.mvabund.2, 
                                 newdata = predict.log.vLG_abund, 
                                 type = "response", se.fit = TRUE) 
colnames(log.vLG_abund.predict$se.fit) <- colnames(log.vLG_abund.predict$fit) # set column names for SE
log.vLG_abund.fit.df <- cbind.data.frame(log.vLG_abund.predict$fit, 
                                         predict.log.vLG_abund)
log.vLG_abund.SE.df <- cbind.data.frame(log.vLG_abund.predict$se.fit,
                                        predict.log.vLG_abund)
log.vLG_abund.predict.df <- rbind.data.frame(log.vLG_abund.fit.df,
                                             log.vLG_abund.SE.df) %>%
  mutate(type = c(rep("predict", row.count),
                  rep("SE", row.count))) %>%
  select(-log.1.aSG_abund, -log.1.rG_abund, -vLG.height.mean) %>%
  melt(id.vars = c("type","log.vLG_abund")) %>%
  spread(type, value)

ggplot(data = filter(log.vLG_abund.predict.df, 
                     variable %in% c("vLG_Platy","vLG_Mesopol","vLG_Tory")),
       aes(x = exp(log.vLG_abund), y = predict)) + 
  geom_line(aes(color = variable)) +
  geom_ribbon(aes(ymax = predict + SE, ymin = predict - SE, fill = variable), alpha = 0.25) +
  ylab("Link density (#/branch)") + xlab("Leaf gall midge density (#/branch)")


## generate a data frame that varies on Iteomyia diameter but holds all other variables constant at their means.
predict.vLG.height.mean <- with(full.predictors,
                                data.frame(vLG.height.mean = seq(from = min(vLG.height.mean),
                                                                 to = max(vLG.height.mean),
                                                                 length.out = row.count),
                                           log.vLG_abund = rep(mean(log.vLG_abund), row.count),
                                           log.1.rG_abund = rep(mean(log.1.rG_abund), row.count),
                                           log.1.aSG_abund = rep(mean(log.1.aSG_abund), row.count)))

## predict link abundance based on the Iteomyia diameter variation.
vLG.height.mean.predict <- predict(net.mvabund.2, 
                                   newdata = predict.vLG.height.mean, 
                                   type = "response", se.fit = TRUE) 
colnames(vLG.height.mean.predict$se.fit) <- colnames(vLG.height.mean.predict$fit) # set column names for SE
vLG.height.mean.fit.df <- cbind.data.frame(vLG.height.mean.predict$fit, 
                                           predict.vLG.height.mean)
vLG.height.mean.SE.df <- cbind.data.frame(vLG.height.mean.predict$se.fit,
                                          predict.vLG.height.mean)
vLG.height.mean.predict.df <- rbind.data.frame(vLG.height.mean.fit.df,
                                               vLG.height.mean.SE.df) %>%
  mutate(type = c(rep("predict", row.count),
                  rep("SE", row.count))) %>%
  select(-log.1.aSG_abund, -log.1.rG_abund, -log.vLG_abund) %>%
  melt(id.vars = c("type","vLG.height.mean")) %>%
  spread(type, value)

ggplot(data = filter(vLG.height.mean.predict.df, 
                     variable %in% c("vLG_Platy","vLG_Mesopol","vLG_Tory","rG_Tory")),
       aes(x = vLG.height.mean, y = predict)) + 
  geom_line(aes(color = variable)) +
  geom_ribbon(aes(ymax = predict + SE, ymin = predict - SE, fill = variable), alpha = 0.25) +
  ylab("Link density (#/branch)") + xlab("Leaf gall diameter (mm)")

## generate a data frame that varies on bud gall density but holds all other variables constant at their means.
predict.log.1.rG_abund <- with(full.predictors,
                               data.frame(log.1.rG_abund = seq(from = min(log.1.rG_abund),
                                                               to = max(log.1.rG_abund),
                                                               length.out = row.count),
                                          log.vLG_abund = rep(mean(log.vLG_abund), row.count),
                                          vLG.height.mean = rep(mean(vLG.height.mean), row.count),
                                          log.1.aSG_abund = rep(mean(log.1.aSG_abund), row.count)))

## predict link abundance based on the bud gall variation.
log.1.rG_abund.predict <- predict(net.mvabund.2, 
                                  newdata = predict.log.1.rG_abund, 
                                  type = "response", se.fit = TRUE) 
colnames(log.1.rG_abund.predict$se.fit) <- colnames(log.1.rG_abund.predict$fit) # set column names for SE
log.1.rG_abund.fit.df <- cbind.data.frame(log.1.rG_abund.predict$fit, 
                                          predict.log.1.rG_abund)
log.1.rG_abund.SE.df <- cbind.data.frame(log.1.rG_abund.predict$se.fit,
                                         predict.log.1.rG_abund)
log.1.rG_abund.predict.df <- rbind.data.frame(log.1.rG_abund.fit.df,
                                              log.1.rG_abund.SE.df) %>%
  mutate(type = c(rep("predict", row.count),
                  rep("SE", row.count))) %>%
  select(-log.1.aSG_abund, -vLG.height.mean, -log.vLG_abund) %>%
  melt(id.vars = c("type","log.1.rG_abund")) %>%
  spread(type, value)

ggplot(data = filter(log.1.rG_abund.predict.df, 
                     variable %in% c("rG_Tory","rG_Eulo")),
       aes(x = exp(log.1.rG_abund), y = predict)) + 
  geom_line(aes(color = variable)) +
  geom_ribbon(aes(ymax = predict + SE, ymin = predict - SE, fill = variable), alpha = 0.25) +
  ylab("Link density (#/branch)") + xlab("Bud gall density (#/branch)")

## NOTE THAT IT IS N'T WORKING TO HAVE THE RAW DATA BECAUSE THE SEQUENCE DATA OF THE PREDICTOR VARIABLE MAY HAVE SOME OVERLAPPING VALUES IN ITS ORIGINAL FORM.

## generate a data frame that varies on Iteomyia diameter but holds all other variables constant at their means.
#log.1.aSG_abund.seq <- arrange(full.predictors, log.1.aSG_abund)$log.1.aSG_abund
#log.1.aSG_abund.OG <- arrange(full.predictors, log.1.aSG_abund)[ ,interaxns_noPont]
predict.log.1.aSG_abund <- with(full.predictors,
                                data.frame(log.1.aSG_abund = seq(from = min(log.1.aSG_abund),
                                                                 to = max(log.1.aSG_abund),
                                                                 length.out = row.count),
                                           log.vLG_abund = rep(mean(log.vLG_abund), row.count),
                                           vLG.height.mean = rep(mean(vLG.height.mean), row.count),
                                           log.1.rG_abund = rep(mean(log.1.rG_abund), row.count)))

## predict link abundance based on the Iteomyia diameter variation.
log.1.aSG_abund.predict <- predict(net.mvabund.2, 
                                   newdata = predict.log.1.aSG_abund, 
                                   type = "response", se.fit = TRUE) 
colnames(log.1.aSG_abund.predict$se.fit) <- colnames(log.1.aSG_abund.predict$fit) # set column names for SE
#log.1.aSG_abund.OG.df <- cbind.data.frame(log.1.aSG_abund.OG,
#                                predict.log.1.aSG_abund)
log.1.aSG_abund.fit.df <- cbind.data.frame(log.1.aSG_abund.predict$fit, 
                                           predict.log.1.aSG_abund)
log.1.aSG_abund.SE.df <- cbind.data.frame(log.1.aSG_abund.predict$se.fit,
                                          predict.log.1.aSG_abund)
log.1.aSG_abund.predict.df <- rbind.data.frame(#log.1.aSG_abund.OG.df,
  log.1.aSG_abund.fit.df,
  log.1.aSG_abund.SE.df) %>%
  mutate(type = c(#rep("OG", row.count),
    rep("predict", row.count),
    rep("SE", row.count))) %>%
  select(-log.1.rG_abund, -vLG.height.mean, -log.vLG_abund) %>%
  melt(id.vars = c("type","log.1.aSG_abund")) %>%
  spread(type, value) # was getting an error about duplicate row identifiers, but it seems to be okay..

ggplot(data = filter(log.1.aSG_abund.predict.df, 
                     variable %in% c("aSG_Tory")),
       aes(x = exp(log.1.aSG_abund), y = predict)) + 
  geom_line(aes(color = variable)) +
  geom_ribbon(aes(ymax = predict + SE, ymin = predict - SE, fill = variable), alpha = 0.25) +
  ylab("Link density (#/branch)") + xlab("apical-Stem gall density (#/branch)")

# incorrectly plot the data with all of the points...not controlling for any other predictors though which are important...
test <- dplyr::select(full.predictors, vLG.height.mean, log.vLG_abund, log.1.rG_abund, vLG_Tory, vLG_Platy, vLG_Mesopol, rG_Tory, rG_Eulo)
test2 <- melt(as.data.frame(test), id.vars = c("vLG.height.mean", "log.vLG_abund", "log.1.rG_abund"))
test2 <- mutate(test2, value.bin = ifelse(value > 0, 1, 0))
test2 <- mutate(test2, log.vLG.cut = ifelse(log.vLG_abund < 1.343, 0, 1))

summary(full.predictors$log.vLG_abund)

library(MASS)
ggplot(filter(test2, variable %in% c("vLG_Tory","vLG_Platy","vLG_Mesopol")), aes(x = vLG.height.mean, y = value, color = factor(log.vLG.cut))) + geom_point() + geom_smooth(method = "glm.nb") + facet_wrap( ~ variable, ncol = 3)
ggplot(filter(test2, variable %in% c("vLG_Tory","vLG_Platy","vLG_Mesopol","rG_Tory")), aes(x = vLG.height.mean, y = value.bin)) + geom_point() + geom_smooth(method = "glm", family = "binomial") + facet_wrap( ~ variable, ncol = 1)

test3 <- glmer(vLG_Platy ~ vLG.height.mean + log(vLG_abund) + (1|Genotype), full.predictors, family = poisson)
overdisp_fun(test3) # no overdispersion
summary(test3)
plot(test3) # residuals don't look great though
visreg(test3, scale = "response", partial = TRUE) # why is the vLG_Platy point so high with the vLG_abund? It is impossible...but may suggest poor model fit, since these are the partial residual plots (ref: https://github.com/pbreheny/visreg/issues/1). Note however, that this issue is much improved when I used a genearlized linear mixed effects models
```

```{r}
full.df$vLG_link_rich <- rowSums(with(full.df, cbind(vLG_Platy, vLG_Tory, vLG_Mesopol, vLG_Mymarid)) > 0)
full.df$rG_link_rich <- rowSums(with(full.df, cbind(rG_Platy, rG_Tory, rG_Mesopol, rG_Lestodip)) > 0)
full.df$vLG_link_div <- diversity(with(full.df, cbind(vLG_Platy, vLG_Tory, vLG_Mesopol, vLG_Mymarid)))
full.df.sub <- full.df %>%
  filter(vLG_abund > 0) %>%
  mutate(vLG_abund.sub.log = log(vLG_abund))
ggplot(full.df.sub, aes(x = vLG_abund, y = vLG_link_rich)) + geom_point() + geom_smooth(method = "glm", family = "poisson") + scale_y_continuous(lim = c(0,3.5))

summary(glm(vLG_link_rich ~ log(vLG_abund), filter(full.df, vLG_abund > 0), family = poisson))
summary(glm(rG_link_rich ~ (rG_abund), filter(full.df, rG_abund > 0), family = poisson), scale = "response")
exp(0.11037)

rare.df <- full.df %>%
  mutate(vLG_link_abund = vLG_Platy + vLG_Tory + vLG_Mesopol + vLG_Mymarid) %>%
  filter(vLG_link_abund > 1) %>%
  mutate(vLG_link_rarerich = rarefy(cbind(vLG_Platy, vLG_Tory, vLG_Mesopol, vLG_Mymarid), 2)-1)
ggplot(rare.df, aes(x = vLG_abund, y = vLG_link_rarerich)) + geom_point() + geom_smooth(method = "lm")
summary(lm(vLG_link_rarerich ~ vLG_abund, rare.df))

anova(glm(vLG_link_rich ~ Genotype, full.df, family = poisson), test = "Chi")
vLG_link_rich.abund <- glm(vLG_link_rich ~ vLG_abund, filter(full.df, vLG_abund > 0), family = poisson)
summary(vLG_link_rich.abund)
dfun(vLG_link_rich.abund)
visreg(vLG_link_rich.abund, scale = "response")
plot(vLG_link_rich ~ vLG_abund, filter(full.df, vLG_abund > 0))

plus.links <- rowSums(full.df[ ,interaxns_noPont])
rarecurve(filter(full.df[,interaxns_noPont], plus.links > 0))
visreg(glm(link_richness ~ log(link_abund), filter(full.df, link_abund > 0), family = poisson), scale = "response")
```

```{r Interaction predictions}
# calculate the mean values, on the original scale, of all predictor variables
no.var.data <- with(full.predictors,
                    data.frame(vLG_abund = mean(vLG_abund),
                               rG_abund = mean(rG_abund),
                               aSG_abund = mean(aSG_abund),
                               vLG.height.mean = mean(vLG.height.mean)))

# http://www.ats.ucla.edu/stat/mult_pkg/faq/general/log_transformed_regression.htm Use this web page for interpreting log and non-transformed coefficients in regression.
library(visreg)
par(mfrow = c(1,2))
vLG_Platy.glm <- MASS::glm.nb(vLG_Platy ~ log(vLG_abund) + vLG.height.mean, full.predictors); summary(vLG_Platy.glm); jensen_magnitude(vLG_Platy.glm, no.var.data)
visreg(vLG_Platy.glm, scale = "response")
#plot(vLG_Platy.glm)

vLG_Mesopol.glm <- MASS::glm.nb(vLG_Mesopol ~ log(vLG_abund) + vLG.height.mean, full.predictors); summary(vLG_Mesopol.glm); jensen_magnitude(vLG_Mesopol.glm, no.var.data)
#plot(vLG_Mesopol.glm)
visreg(vLG_Mesopol.glm, scale = "response")

vLG_Tory.glm <- MASS::glm.nb(vLG_Tory ~ vLG_abund+vLG.height.mean, full.predictors); summary(vLG_Tory.glm); jensen_magnitude(vLG_Tory.glm, no.var.data) # gall height is marginally significant, but the AIC is slightly better than the model with only vLG_abund.
visreg(vLG_Tory.glm, scale = "response")

rG_Tory.glm <- glm(rG_Tory ~ rG_abund + vLG.height.mean, full.predictors, family = quasipoisson); summary(rG_Tory.glm); jensen_magnitude(rG_Tory.glm, no.var.data) # underdispersed, which is why I used quasipoisson.
#plot(rG_Tory.glm)
visreg(rG_Tory.glm, scale = "response")

rG_Platy.glm <- glm(rG_Platy ~ rG_abund, full.predictors, family = quasipoisson); summary(rG_Platy.glm); jensen_magnitude(rG_Platy.glm, no.var.data)  # underdispersed, which is why I used quasipoisson.

rG_Eulo.glm <- glm(rG_Eulo ~ rG_abund, full.predictors, family = quasipoisson); summary(rG_Eulo.glm); jensen_magnitude(rG_Eulo.glm, no.var.data) # underdispersed, which is why I used quasipoisson

aSG_Tory.glm <- glm(aSG_Tory ~ aSG_abund, full.predictors, family = quasipoisson); summary(aSG_Tory.glm); jensen_magnitude(aSG_Tory.glm, no.var.data) # underdispersed, so I used quasipoisson. 

rG_Lestodip.glm <- glm(rG_Lestodip ~ log(rG_abund + 1), full.predictors, family = quasipoisson)
summary(rG_Lestodip.glm); jensen_magnitude(rG_Lestodip.glm, no.var.data) # rG_abund is marginally significant.

## intercept only models.
vLG_Eulo.glm <- glm(vLG_Eulo ~ 1, full.predictors, family = poisson); summary(vLG_Eulo.glm); jensen_magnitude(vLG_Eulo.glm, no.var.data) # not a lot of support for including gall height as a predictor, althought it was marginally significant in multivariate model, so I dropped it from the model

SG_Platy.glm <- glm(SG_Platy ~ 1, full.predictors, family = poisson); summary(SG_Platy.glm); jensen_magnitude(SG_Platy.glm, no.var.data) # no evidence of vLG_abund being a clear predictor, so we dropped it from the model

rG_Mesopol.glm <- glm(rG_Mesopol ~ 1, full.predictors, family = poisson); summary(rG_Mesopol.glm); jensen_magnitude(rG_Mesopol.glm, no.var.data)

vLG_Mymarid.glm <- glm(vLG_Mymarid ~ 1, full.predictors, family = poisson); summary(vLG_Mymarid.glm); jensen_magnitude(vLG_Mymarid.glm, no.var.data) # vLG_abund is a significant predictor, but the visreg predictions look so wonky that I don't really trust it, so I just modelled it with an intercept
```

```{r visreg plots}
# Platygaster interactions
layout(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = FALSE))
visreg(vLG_Platy.glm, scale = "response", ask = F)
plot(x = range(full.predictors$vLG_abund), y = c(0, 1), type = "n", ylab = "SG_Platy", xlab = "vLG_abund"); abline(h = exp(coef(SG_Platy.glm)), col = "steelblue", lwd = 3) # intercept only models
visreg(rG_Platy.glm, scale = "response", ask = F)

# Mesopolobus interactions
layout(matrix(c(1, 2, 3, 0), nrow = 2, ncol = 2, byrow = FALSE))
visreg(vLG_Mesopol.glm, scale = "response", ask = F)
plot(x = range(full.predictors$rG_abund), y = c(0, 1), type = "n", ylab = "rG_Mesopol", xlab = "rG_abund"); abline(h = exp(coef(rG_Mesopol.glm)), col = "steelblue", lwd = 3) # intercept only models

# Torymus interactions
layout(matrix(c(1, 2, 3, 4, 5, 0), nrow = 2, ncol = 3, byrow = FALSE))
visreg(vLG_Tory.glm, scale = "response", ask = F)
visreg(rG_Tory.glm, scale = "response", ask = F)
visreg(aSG_Tory.glm, scale = "response", ask = F)

# Eulophid interactions and other less abundant ones
layout(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = FALSE))
plot(x = range(full.predictors$vLG_abund), y = c(0, 1), type = "n", ylab = "vLG_Eulo", xlab = "vLG_abund"); abline(h = exp(coef(vLG_Eulo.glm)), col = "steelblue", lwd = 3)
visreg(rG_Eulo.glm, scale = "response", ask = F)
plot(x = range(full.predictors$vLG_abund), y = c(0, 1), type = "n", ylab = "vLG_Mymarid", xlab = "vLG_abund"); abline(h = exp(coef(vLG_Mymarid.glm)), col = "steelblue", lwd = 3)
visreg(rG_Lestodip.glm, scale = "response", ask = F)

par(mfrow = c(1,1))
```

```{r table of jensen effects}
jensen.df <- rbind_all(list(jensen_magnitude(vLG_Platy.glm, no.var.data),
                            jensen_magnitude(vLG_Mesopol.glm, no.var.data),
                            jensen_magnitude(vLG_Tory.glm, no.var.data),
                            jensen_magnitude(vLG_Eulo.glm, no.var.data),
                            jensen_magnitude(vLG_Mymarid.glm, no.var.data),
                            jensen_magnitude(rG_Tory.glm, no.var.data),
                            jensen_magnitude(rG_Mesopol.glm, no.var.data),
                            jensen_magnitude(rG_Platy.glm, no.var.data),
                            jensen_magnitude(rG_Mesopol.glm, no.var.data),
                            jensen_magnitude(rG_Lestodip.glm, no.var.data),
                            jensen_magnitude(aSG_Tory.glm, no.var.data),
                            jensen_magnitude(SG_Platy.glm, no.var.data)))
jensen.df <- mutate(jensen.df, Interactions = c("vLG_Platy","vLG_Mesopol","vLG_Tory","vLG_Eulo","vLG_Mymarid","rG_Tory","rG_Mesopol","rG_Platy","rG_Mesopol","rG_Lestodip","aSG_Tory","SG_Platy"))
plot(log(Magnitude.Difference) ~ log(mean.no.var), filter(jensen.df, Magnitude.Difference != "1"))
visreg(lm(log(Magnitude.Difference) ~ log(mean.no.var), filter(jensen.df, Magnitude.Difference != "1")), scale = "response", xtrans = log, xlab = "log(mean.no.var)")

interactions.split <- colsplit(jensen.df$Interactions, "_", names = c("gall","parasitoid"))
jensen.df <- cbind(jensen.df, interactions.split)

var.web <- cast(jensen.df, gall ~ parasitoid, mean, value = "mean.with.var")
rownames(var.web) <- var.web$gall
var.web <- as.matrix.data.frame(var.web[ ,-1])

no.var.web <- cast(jensen.df, gall ~ parasitoid, mean, value = "mean.no.var")
rownames(no.var.web) <- no.var.web$gall
no.var.web <- as.matrix.data.frame(no.var.web[ ,-1])

library(bipartite)
plotweb(var.web, method = "normal")
plotweb(no.var.web, method = "normal")
```

```{r remove weak interactions}
colSums(full.predictors[ ,interaxns_noPont]) # remove rG_Lestodip, rG_Mesopol, vLG_Mymarid, and SG_Platy. SG_Platy had 9 interactions, but all of these interactions came from a single, multichambered gall.

net.trait.noweak <- mvabund(full.predictors[ ,c("aSG_Tory", "rG_Eulo", "rG_Tory", "vLG_Eulo", "vLG_Mesopol","vLG_Platy", "vLG_Tory")])

net.mvabund.2.noweak <- manyglm(net.trait.noweak ~ log.vLG_abund + log.vLG.height.mean + 
                                  log.1.rG_abund + log.1.aSG_abund,
                                data = full.predictors,
                                family = "negative.binomial")
```

```{r}
## Identify the effect of phenotypic variation (gall density and gall size) on food web structure

## create data with the same mean, but NO VARIANCE in predictor variables. For some reason, the predict.manyglm needs a dataset that is the same dimension as the original one.
row.count <- dim(full.predictors)[1]
no.var.data <- with(full.predictors, 
                    data.frame(log.vLG_abund = rep(log(mean(vLG_abund)), row.count),
                               log.1.rG_abund = rep(log(mean(rG_abund)+1), row.count),
                               log.1.aSG_abund = rep(log(mean(aSG_abund)+1), row.count),
                               log.vLG.height.mean = rep(log(mean(vLG.height.mean)), row.count)))
# increase mean by one standard deviation. With this type of code, I should be able to make predictions for how evolution will affect the food web.
no.var.data.vLG.size.Increase <- with(full.predictors, 
                                      data.frame(log.vLG_abund = rep(log(mean(vLG_abund)), row.count),
                                                 log.1.rG_abund = rep(log(mean(rG_abund)+1), row.count),
                                                 log.1.aSG_abund = rep(log(mean(aSG_abund)+1), row.count),
                                                 log.vLG.height.mean = rep(log(mean(vLG.height.mean)+sd(vLG.height.mean)), row.count)))

jensen.net.mvabund.2 <- jensen_magnitude_manyglm(net.mvabund.2, original.data = full.predictors,
                                                 comparison.data = no.var.data)
jensen.net.mvabund.2.noweak <- jensen_magnitude_manyglm(net.mvabund.2.noweak, original.data = full.predictors, comparison.data = no.var.data)
jensen.net.mvabund.2.meanIncrease <- jensen_magnitude_manyglm(net.mvabund.2, original.data = no.var.data, comparison.data = no.var.data.vLG.size.Increase)

# used a wilcoxon signed rank test because the differences were not normally distributed
with(jensen.net.mvabund.2, hist(var.predict - no.var.predict, breaks = seq(0,2,0.01)))
with(jensen.net.mvabund.2, wilcox.test(var.predict, no.var.predict, 
                                       paired = TRUE, conf.int = TRUE))
magnitude.lm <- lm(log(Percent.Difference) ~ log(no.var.predict), jensen.net.mvabund.2)
summary(magnitude.lm); visreg(magnitude.lm, scale = "response", xtrans = log)

# I get the same qualitative results even after removing the weakest interactions from the dataset
with(jensen.net.mvabund.2.noweak, wilcox.test(var.predict, no.var.predict, 
                                              paired = TRUE, conf.int = TRUE)) 
magnitude.lm.noweak <- lm(log(Percent.Difference) ~ log(no.var.predict), 
                          jensen.net.mvabund.2.noweak)
summary(magnitude.lm.noweak); visreg(magnitude.lm.noweak, scale = "response", xtrans = log)
```

```{r effects of variation on food web structure}
library(bipartite)
jensen.to.use <- jensen.net.mvabund.2.noweak
interactions.split <- colsplit(jensen.to.use$response.variables, 
                               "_", names = c("gall","parasitoid"))
jensen.to.use <- cbind(jensen.to.use, interactions.split)

var.web <- cast(jensen.to.use, gall ~ parasitoid, mean, value = "var.predict")
rownames(var.web) <- var.web$gall
var.web <- as.matrix.data.frame(var.web[ ,-1])

no.var.web <- cast(jensen.to.use, gall ~ parasitoid, mean, value = "no.var.predict")
rownames(no.var.web) <- no.var.web$gall
no.var.web <- as.matrix.data.frame(no.var.web[ ,-1])

par(mfrow = c(2,1))
plotweb(var.web, method = "normal")
plotweb(no.var.web, method = "normal")

cbind(networklevel(var.web, 
                   index = c("linkage density", "interaction evenness")),
      networklevel(no.var.web,
                   index = c("linkage density", "interaction evenness")))
```

```{r}
vLG.ptized.glm <- glm(cbind(vLG_parasitized, vLG_abund - vLG_parasitized) ~ vLG.height.mean, full.df, family = binomial)
summary(vLG.ptized.glm)
visreg(vLG.ptized.glm, scale = "response")
vLG.ptized.no.var.data <- with(full.predictors, data.frame(vLG.height.mean = mean(vLG.height.mean)))
jensen_magnitude(vLG.ptized.glm, no.var.data = vLG.ptized.no.var.data)
#plot(vLG.ptized.glm)

vLG.Platy.glm <- glm(cbind(vLG_Platy, vLG_abund - vLG_Platy) ~ vLG.height.mean, full.df, family = binomial)
summary(vLG.Platy.glm)
visreg(vLG.Platy.glm, scale = "response")
vLG.Platy.no.var.data <- with(full.predictors, data.frame(vLG.height.mean = mean(vLG.height.mean)))
jensen_magnitude(vLG.Platy.glm, no.var.data = vLG.Platy.no.var.data)

vLG.Tory.glm <- glm(cbind(vLG_Tory, vLG_abund - vLG_Tory) ~ vLG.height.mean, full.df, family = binomial)
summary(vLG.Tory.glm)
visreg(vLG.Tory.glm, scale = "response")
vLG.Tory.no.var.data <- with(full.predictors, data.frame(vLG.height.mean = mean(vLG.height.mean)))
jensen_magnitude(vLG.Tory.glm, no.var.data = vLG.Tory.no.var.data)

#plot(vLG_Platy ~ Genotype, full.df)
#plot(vLG_Tory ~ Genotype, full.df)

vLG.ptized.glmer <- glmer(cbind(vLG_parasitized, vLG_abund - vLG_parasitized) ~ log(vLG.height.mean) + (1|Genotype), full.df, family = binomial)
summary(vLG.ptized.glmer)
visreg(vLG.ptized.glmer, scale = "response")
mean(predict(vLG.ptized.glmer, type = "response"))
#mean(predict(vLG.ptized.glmer, newdata = vLG.ptized.no.var.data, type = "response"))

vLG_Platy.glmer <- glmer(vLG_Platy ~ log(vLG_abund) + log(vLG.height.mean) + (1|Genotype), full.df, family = poisson)
summary(vLG_Platy.glmer)
overdisp_fun(vLG_Platy.glmer) # no overdispersion when I model genotype as a randome effect.
visreg(vLG_Platy.glmer, scale = "response")
```

```{r}
## Examine whether the proportion of galls parasitized varies among willow genotypes

## vLG parasitism
vLG.ptized.df <- full.predictors %>%
  mutate(vLG_parasitized = vLG_Tory + vLG_Eulo + vLG_Platy + vLG_Mesopol + vLG_Mymarid)
table(vLG.ptized.df$Genotype) # data is unbalanced, but this doesn't qualitatively affect the outcome of the glm.
vLG.ptized.3plus.genos <- c("*","B","D","E","F","H","I","K","L","O","S","V","X","Y","Z")
colSums(vLG.ptized.df[ ,c("vLG_parasitized","vLG_abund","vLG_Tory","vLG_Eulo","vLG_Platy",
                          "vLG_Mesopol","vLG_Mymarid")])

# interesting, only marginally significant with quasibinomial
vLG.ptized.bin <- glm(cbind(vLG_parasitized, vLG_abund - vLG_parasitized) ~ Genotype,
                      data = filter(vLG.ptized.df, Genotype %in% vLG.ptized.3plus.genos), 
                      family = 'binomial')
summary(vLG.ptized.bin)
dfun(vLG.ptized.bin)
anova(vLG.ptized.bin, test = "LR")
require(dispmod)
anova(glm.binomial.disp(vLG.ptized.bin), test = "LR") # still significant after accounting for overdispersion
visreg(vLG.ptized.bin, scale = "response")

vLG.ptized.pred <- glm(cbind(vLG_parasitized, vLG_abund - vLG_parasitized) ~ vLG.height.mean,
                       data = vLG.ptized.df,
                       family = 'binomial')
summary(vLG.ptized.pred)
1-exp(coef(vLG.ptized.pred)[2]) # 25% reduction in the oddes of vLG being parasitized with every one unit increase in gall size.
dfun(vLG.ptized.pred) # overdispersion
anova(glm.binomial.disp(vLG.ptized.pred), test = "LR")
require(visreg)
visreg(vLG.ptized.pred, scale = "response")


vLG.ptized.bin.glmer <- glmer(cbind(vLG_parasitized, vLG_abund - vLG_parasitized) ~ (1|Genotype),
                              data = vLG.ptized.df, # not removing Genotypes with less than 2 replicates
                              family = 'binomial')
overdisp_fun(vLG.ptized.bin.glmer)
summary(vLG.ptized.bin.glmer)

test.2 <- rpt.binomGLMM.add(y = with(vLG.ptized.df, cbind(vLG_parasitized, vLG_abund - vLG_parasitized)), groups = vLG.ptized.df$Genotype) 
test.2
library(emdbook)
library(MASS)
vLG.ptized.df <- mutate(vLG.ptized.df, 
                        plant.position = factor(plant.position))
library(lme4)
sc.abund <- scale(vLG.ptized.df$vLG_abund)
sc.height <- scale(vLG.ptized.df$vLG.height.mean)
test.4 <- glmer(cbind(vLG_parasitized, vLG_abund - vLG_parasitized) ~  vLG.height.mean + (1|Genotype),
                data = vLG.ptized.df,
                family = 'binomial')
summary(test.4)
AIC(test.4)
confint(test.4, method = "profile") # variance > 0!
overdisp_fun(test.4) # no clear evidence of overdispersion
ggQQ_ranef(ranef(test.4)$Genotype$'(Intercept)') # and the residuals look great!
visreg(test.4, scale = "response")
R.link  <- .6411 / (.6411 +pi^2 /3) # heritability is reasonable as well
R.link.low  <- .1153673 / (.1153673 +pi^2 /3) 
R.link.high  <- 1.44 / (1.44 +pi^2 /3) 
profile(test.4, which = 1)
summary(test.4)
test.3 <- rpt.binomGLMM.multi(y = with(vLG.ptized.df, cbind(vLG_parasitized, vLG_abund - vLG_parasitized)), groups = vLG.ptized.df$Genotype) 

library(RLRsim)

#plot(vLG.ptized.bin)
summary(vLG.ptized.bin)
dfun(vLG.ptized.bin)
anova(vLG.ptized.bin, test = "LR") # still significant after retaining Genotypes with 3 plus replicates.

vLG_ptized.bin.pred <- gam(cbind(vLG_parasitized, vLG_abund - vLG_parasitized) ~ s(vLG.height.mean),
                           data = full.df, 
                           family = 'binomial')
summary(vLG_ptized.bin.pred)
AIC(vLG_ptized.bin.pred)
visreg(vLG_ptized.bin.pred, scale = "response")
#jensen.single(vLG_ptized.bin.pred)

# vLG_Platy test
vLG_Platy.bin <- glm(cbind(vLG_Platy, vLG_abund - vLG_Platy) ~ Genotype,
                     data = filter(full.df, Genotype %in% vLG.ptized.3plus.genos), 
                     family = 'quasibinomial')
plot(vLG_Platy.bin)
summary(vLG_Platy.bin)
anova(vLG_Platy.bin, test = "F")

library(dispmod)
vLG_Platy.bin.pred <- glm(cbind(vLG_Platy, vLG_abund - vLG_Platy) ~ vLG.height.mean,
                          data = full.predictors, 
                          family = 'binomial')
dfun(vLG_Platy.bin.pred)
summary(vLG_Platy.bin.pred)
glm.binomial.disp(vLG_Platy.bin.pred)
AIC(vLG_Platy.bin.pred)
plot(vLG_Platy.bin.pred)
visreg(vLG_Platy.bin.pred, scale = "response")
jensen.single(vLG_Platy.bin.pred)
glm.binomial.disp(vLG_Platy.bin.pred)
qAIC(vLG_Platy.bin.pred, dispersion = dfun(vLG_Platy.bin.pred))

dfun(vLG_Platy.bin.pred)

vLG_Platy.bin.pred.glmer <- glmer(cbind(vLG_Platy, vLG_abund - vLG_Platy) ~ vLG.height.mean + (1|Genotype),
                                  data = full.df, 
                                  family = 'binomial')
summary(vLG_Platy.bin.pred.glmer)
overdisp_fun(vLG_Platy.bin.pred.glmer)
AIC(vLG_Platy.bin.pred.glmer)
plot(vLG_Platy.bin.pred.glmer)
ggQQ_ranef(ranef(vLG_Platy.bin.pred.glmer)$Genotype$'(Intercept)') # not too bad.

# AIC gam interaction 207.6398, AIC = 204 with main effect model in gam. 213 with main model

## Identify the effect of phenotypic variation (gall density and gall size) on food web structure

jensen.single(final.vLG_Platy)
jensen.single(vLG_Mesopol.bin.pred)
jensen.single(vLG_Tory.bin.pred)

# vLG_Mesopol test
vLG_Mesopol.bin <- glm(cbind(vLG_Mesopol, vLG_abund - vLG_Mesopol) ~ Genotype,
                       data = filter(full.df, Genotype %in% vLG.ptized.3plus.genos), 
                       family = 'binomial')
plot(vLG_Mesopol.bin)
summary(vLG_Mesopol.bin)
dfun(vLG_Mesopol.bin)
anova(vLG_Mesopol.bin, test = "LR")

vLG_Mesopol.bin.pred <- gam(cbind(vLG_Mesopol, vLG_abund - vLG_Mesopol) ~ vLG.height.mean,
                            data = full.df, 
                            family = 'quasibinomial')
summary(vLG_Mesopol.bin.pred)
visreg(vLG_Mesopol.bin.pred, xvar = "vLG.height.mean", scale = "response")
dfun(vLG_Mesopol.bin.pred)

jensen.single(vLG_Mesopol.bin.pred)

# vLG_Tory test
vLG_Tory.bin <- glm(cbind(vLG_Tory, vLG_abund - vLG_Tory) ~ Genotype,
                    data = filter(full.df, Genotype %in% vLG.ptized.3plus.genos), 
                    family = 'binomial')
plot(vLG_Tory.bin)
summary(vLG_Tory.bin)
anova(vLG_Tory.bin, test = "LR")

vLG_Tory.bin.pred <- glm(cbind(vLG_Tory, vLG_abund - vLG_Tory) ~ vLG_abund + vLG.height.mean,
                         data = full.df, 
                         family = 'binomial')
summary(vLG_Tory.bin.pred)
AIC(vLG_Tory.bin.pred) # 136.1538 is the lowest score. Main effect model
visreg(vLG_Tory.bin.pred, scale = "response")

# vLG_Eulo test
vLG_Eulo.bin <- glm(cbind(vLG_Eulo, vLG_abund - vLG_Eulo) ~ Genotype,
                    data = filter(full.df, Genotype %in% vLG.ptized.3plus.genos), 
                    family = 'binomial')
#plot(vLG_Eulo.bin)
summary(vLG_Eulo.bin)
anova(vLG_Eulo.bin, test = "LR") # not significant

# not testing vLG_Mymarid because it is in such low abundance
predict.df <- data.frame(Genotype = vLG.ptized.3plus.genos)
vLG.genotype.predictions <- cbind.data.frame(predict.df, 
                                             vLG.Proportion.Ptized = predict(vLG.ptized.bin,
                                                                             newdata = predict.df,
                                                                             type = "response"))
summary(vLG.genotype.predictions)

ggplot(filter(vLG.ptized.df, Genotype %in% vLG.ptized.3plus.genos),
       aes(x = Genotype, y = vLG_parasitized/vLG_abund)) +
  geom_point(aes(size = vLG_abund), color = "grey", shape = 1,
             position = position_jitter(height = 0, width = 0.1)) + 
  scale_size(range = c(3,15)) + 
  geom_point(data = vLG.genotype.predictions, 
             aes(x = Genotype, y = vLG.Proportion.Ptized),
             color = "steelblue", size = 30, shape = "--") +
  theme_classic() +
  ylab("Proportion of galls parasitized")
range(filter(vLG.ptized.df, Genotype %in% vLG.ptized.3plus.genos)$vLG_abund)

## rG parasitism
rG.ptized.df <- full.df %>%
  filter(rG_abund > 0)
table(rG.ptized.df$Genotype)
rG.ptized.3plus.genos <- c("*","A","D", "F", "G", "H","I","K","O","Q", "R", "S","X", "Y", "Z") 
colSums(rG.ptized.df[ ,c("rG_parasitized","rG_abund","rG_Tory","rG_Eulo","rG_Platy",
                         "rG_Mesopol","rG_Lestodip")]) # too few rG_eulo, rG_Platy, rG_Mesopol, and rG_Lestodip to test for differences in parasitism among genotypes.

rG.ptized.bin <- glm(cbind(rG_parasitized, rG_abund - rG_parasitized) ~ Genotype,
                     data = filter(full.df, Genotype %in% rG.ptized.3plus.genos), 
                     family = 'binomial')
#plot(rG.ptized.bin) # residuals not great
summary(rG.ptized.bin)
anova(rG.ptized.bin, test = "LR") # marginally significant after retaining Genotypes with 3 plus replicates.

# rG_Tory test
rG_Tory.bin <- glm(cbind(rG_Tory, rG_abund - rG_Tory) ~ Genotype,
                   data = filter(full.df, Genotype %in% rG.ptized.3plus.genos), 
                   family = 'binomial')
plot(rG_Tory.bin)
summary(rG_Tory.bin)
anova(rG_Tory.bin, test = "LR") # not significant

rG_Eulo.bin <- glm(cbind(rG_Eulo, rG_abund - rG_Eulo) ~ Genotype,
                   data = filter(full.df, Genotype %in% rG.ptized.3plus.genos), 
                   family = 'binomial')
plot(rG_Eulo.bin)
summary(rG_Eulo.bin)
anova(rG_Eulo.bin, test = "LR") # not significant

predict.df <- data.frame(Genotype = rG.ptized.3plus.genos)
rG.genotype.predictions <- cbind.data.frame(predict.df, 
                                            rG.Proportion.Ptized = predict(rG.ptized.bin,
                                                                           newdata = predict.df,
                                                                           type = "response"))
summary(rG.genotype.predictions)

ggplot(filter(rG.ptized.df, Genotype %in% rG.ptized.3plus.genos),
       aes(x = Genotype, y = rG_parasitized/rG_abund)) +
  geom_point(aes(size = rG_abund), color = "grey", shape = 1,
             position = position_jitter(height = 0, width = 0.1)) + 
  scale_size(range = c(3,15)) + 
  geom_point(data = rG.genotype.predictions, 
             aes(x = Genotype, y = rG.Proportion.Ptized),
             color = "steelblue", size = 30, shape = "--") +
  theme_classic() +
  ylab("Proportion of galls parasitized")
range(filter(rG.ptized.df, Genotype %in% rG.ptized.3plus.genos)$rG_abund)

## aSG parasitism
aSG.ptized.df <- full.df %>%
  filter(aSG_abund > 0)
table(aSG.ptized.df$Genotype)
aSG.ptized.3plus.genos <- c("G","I","K","Q") 

aSG.ptized.bin <- glm(cbind(aSG_parasitized, aSG_abund - aSG_parasitized) ~ Genotype,
                      data = filter(full.df, Genotype %in% aSG.ptized.3plus.genos), 
                      family = 'binomial')
plot(aSG.ptized.bin) # residuals not great
summary(aSG.ptized.bin)
anova(aSG.ptized.bin, test = "LR") # marginally significant after retaining Genotypes with 3 plus replicates.
predict.df <- data.frame(Genotype = aSG.ptized.3plus.genos)
aSG.genotype.predictions <- cbind.data.frame(predict.df, 
                                             aSG.Proportion.Ptized = predict(aSG.ptized.bin,
                                                                             newdata = predict.df,
                                                                             type = "response"))
summary(aSG.genotype.predictions)

ggplot(filter(aSG.ptized.df, Genotype %in% aSG.ptized.3plus.genos),
       aes(x = Genotype, y = aSG_parasitized/aSG_abund)) +
  geom_point(aes(size = aSG_abund), color = "grey", shape = 1,
             position = position_jitter(height = 0, width = 0.1)) + 
  scale_size(range = c(3,15)) + 
  geom_point(data = aSG.genotype.predictions, 
             aes(x = Genotype, y = aSG.Proportion.Ptized),
             color = "steelblue", size = 30, shape = "--") +
  theme_classic() +
  ylab("Proportion of galls parasitized")
range(filter(aSG.ptized.df, Genotype %in% aSG.ptized.3plus.genos)$aSG_abund)

## SG parasitism
SG.ptized.df <- full.df %>%
  filter(SG_abund > 0)
table(SG.ptized.df$Genotype) # no Genotypes with 3+ replicates, so I'm unable to evaluate whether SG parasitism defers among willow genotypes. 

##
uni.check.bin <- glm(vLG_Platy/vLG_abund ~ vLG_abund*vLG.height.mean, full.predictors,
                     family = binomial,
                     weights = vLG_abund)
AIC(uni.check.bin)
summary(uni.check.bin)
visreg(uni.check.bin, xvar = "vLG_abund", by = "vLG.height.mean", scale = "response")


# NEED TO THINK ABOUT HOW I SCALE THIS UP TO DIFFERENT NUMBERS OF SHOOTS SAMPLED.

#### messing around with different ideas.
## Create mvabund model based on probability of observing an interaction
prob.df <- full.predictors %>%
  mutate(succ.vLG_Platy.succ = vLG_Platy,
         fail.vLG_Platy.fail = vLG_abund - vLG_Platy,
         succ.vLG_Tory.succ = vLG_Tory,
         fail.vLG_Tory.fail = vLG_abund - vLG_Tory,
         succ.vLG_Mesopol.succ = vLG_Mesopol,
         fail.vLG_Mesopol.fail = vLG_abund - vLG_Mesopol,
         succ.vLG_Eulo.succ = vLG_Eulo,
         fail.vLG_Eulo.fail = vLG_abund - vLG_Eulo,
         succ.rG_Platy.succ = rG_Platy,
         fail.rG_Platy.fail = rG_abund - rG_Platy,
         succ.rG_Tory.succ = rG_Tory,
         fail.rG_Tory.fail = rG_abund - rG_Tory,
         succ.rG_Mesopol.succ = rG_Mesopol,
         fail.rG_Mesopol.fail = rG_abund - rG_Mesopol,
         succ.rG_Eulo.succ = rG_Eulo,
         fail.rG_Eulo.fail = rG_abund - rG_Eulo,
         succ.SG_Platy.succ = SG_Platy,
         fail.SG_Platy.fail = SG_abund - SG_Platy,
         succ.aSG_Tory.succ = aSG_Tory,
         fail.aSG_Tory.fail = aSG_abund - aSG_Tory) %>%
  select(succ.vLG_Platy.succ:fail.aSG_Tory.fail) %>%
  mvabund()

vLG.interaction.df <- mvabund(select(full.predictors, vLG_Eulo, vLG_Platy, vLG_Mymarid, vLG_Tory, vLG_Mesopol))
vLG.interaction.df[vLG.interaction.df > 0] <- 1

test <- manyglm(vLG.interaction.df ~ vLG_density*vLG.height.mean, data = full.predictors, family = "binomial")
plot(test)
anova(test, p.uni = "adjusted")
length(prob.df[prob.df > 1] > 0)

test <- MASS::glm.nb(vLG_Platy ~ offset(log(vLG_abund)) + vLG.height.mean, data = full.predictors) #
summary(test)
coef(test)
AIC(test)
#exp(sum(-4.5 + .34*5 -0.25*8))
exp(sum(coef(test)[c(1:3)]*c(1,0.02,8))) # 8 mm gall at a density of 11 per 100 shoots. Would be 6.11 parasitoids per 100 shoots.
predict(test, type = "response")*10
visreg(test, scale = "response")
```

```{r, functions, echo=FALSE}
kable.adonis <- function(adonis.object){
  require(knitr)
  kable(adonis.object$aov.tab[ ,-5], # remove R2 calculation
        col.names = c("df","SS","MS","F","P"), 
        digits = c(0,2,2,2,3))
}

heritability.adonis <- function(adonis.object, data){
  gen.row <- which(rownames(adonis.object$aov.tab) == "Genotype")
  res.row <- which(rownames(adonis.object$aov.tab) == "Residuals")
  gen.MS <- adonis.object$aov.tab$MeanSqs[gen.row]
  res.MS <- adonis.object$aov.tab$MeanSqs[res.row]
  mean.sample.size <- mean(table(data$Genotype))
  gen.variance <- (gen.MS - res.MS)/mean.sample.size
  heritability <- gen.variance/(gen.variance + res.MS) # note that residual MS is equal to the residual variance
  print(heritability)
}

heritability <- function(merMod.object){
  require(lme4)
  require(RLRsim)
  
  var.comps <- as.data.frame(VarCorr(merMod.object))
  genetic.var <- var.comps$vcov[1]
  residual.var <- var.comps$vcov[2]
  heritability <- genetic.var/(genetic.var + residual.var)
  restricted.likelihood.ratio.test <- exactRLRT(merMod.object)
  data.frame(#Model = formula(merMod.object),
    Heritability = round(heritability, 2),
    RLRT_statistic = round(restricted.likelihood.ratio.test$statistic, 2),
    P = round(restricted.likelihood.ratio.test$p.value, 3))
}



jensen.single <- function(model.object, predict.type = "response"){
  model <- model.object
  
  var.predict <- predict(model.object, type = predict.type)
  
  predictor.names <- names(model.object$model)[-1] # omits response variable
  data <- data.frame(model.object$model[ ,predictor.names])
  
  no.var.data <- list()
  for(i in 1:length(data)){
    no.var.data[[i]] <- mean(data[ ,i])
  }
  no.var.data <- as.data.frame(no.var.data)
  names(no.var.data) <- predictor.names
  
  no.var.predict <- predict(model.object, newdata = no.var.data, type = predict.type)
  
  data.frame(mean.with.var = mean(var.predict, na.rm = TRUE),
             mean.no.var = mean(no.var.predict, na.rm = TRUE))
}

jensen.single <- function(model.object, no.var.data, predict.type = "response"){
  model <- model.object
  
  var.predict <- predict(model.object, type = predict.type)
  
  predictor.names <- names(model.object$model)[-1] # omits response variable
  data <- data.frame(model.object$model[ ,predictor.names])
  
  #no.var.data <- list()
  #for(i in 1:length(data)){
  # no.var.data[[i]] <- mean(data[ ,i])
  #}
  #no.var.data <- as.data.frame(no.var.data)
  #names(no.var.data) <- predictor.names
  
  no.var.predict <- predict(model.object, newdata = no.var.data, type = predict.type)
  
  data.frame(mean.with.var = mean(var.predict, na.rm = TRUE),
             mean.no.var = mean(no.var.predict, na.rm = TRUE))
}

```

```{r}
library(psych)
gall.df.cors <- full.df %>%
  mutate(vLG_density = log(vLG_abund+1) - log(shootEst.no18),
         vLG_density2 = log(vLG_abund/shootEst.no18+0.01),
         vLG_density3 = sqrt(vLG_abund/shootEst.no18),
         rG_density = log(vLG_abund+1) - log(shootEst.no18),
         rG_density2 = log(rG_abund/shootEst.no18+0.01),
         rG_density3 = sqrt(rG_abund/shootEst.no18)) %>%
  select(vLG.height.mean, vLG_density:rG_density3)
trait.df.cors <- full.df %>%
  select(Total_Area:flavanonOLES.PC1)

gall.trait.cors <- corr.test(x = trait.df.cors, y = gall.df.cors)
round(gall.trait.cors$r,2)
# further examine all correlations greater than 0.20
# vLG.height.mean: salicortin__A270nm, luteolin.glucoside.5.glu.__A320nm, quercetyin.der1, luteolin.der1, luteolin der2, termulacin, Carbon, sal_tannin.PC1, flavonoOLES.PC1
# vLG_density: Total_Area, Height, D_mean_smoothed, ampelopsin.der__A320nm, myricitrin__A320nm, eriodictyol.7.glucoside__A270nm, luteolin.der1__A320nm, C_N_imputed (sal_tannin.PC1? a bit low). Note though the Nitrogen and C_N_ratio are not good predictors for the "un_imputed" data. 
# rG_density: total_area, height, D_mean_smoothed, ampelopsin.der__A320nm, eriodictyol.7.glucoside__A270nm, luteolin.der1__A320nm, C_N_imputed

library(car)
# need to work on code below.
scatterplotMatrix(cbind.data.frame(gall.df.cors$vLG.height.mean, 
                                   trait.df.cors[ ,"salicortin__A270nm", 
                                                 "luteolin.glucoside.5.glu.__A320nm",
                                                 "quercetin.der1__A320nm", 
                                                 "luteolin.der1__A320nm", 
                                                 "luteolin der2__A320nm", 
                                                 "termulacin__A320nm",
                                                 "Carbon", "sal_tannin.PC1", 
                                                 "flavonoOLES.PC1"]))
```

```{r}
interaction_abund <- rowSums(full.df[ ,interaxns_noPont])
interaction_richness <- rowSums(full.df[ ,interaxns_noPont] > 0)
ggplot(full.df, aes(x = Genotype, y = interaction_abund)) +
  geom_boxplot()

ggplot(full.df, aes(x = Genotype, y = interaction_abund/shootEst.no18)) +
  geom_boxplot() # G, H, M, N, O, P, R, U all have 75% of their observations at zero for interaction density

ggplot(full.df, aes(x = Genotype, y = interaction_richness)) +
  geom_boxplot() # G, H, M, N, O, P, R, U all have 75% of their observations at zero for interaction richness

non.zero.genos <- c("*","A","B","D","E","F","I","J","K","L","Q","S","T","V","W","X","Y","Z")

df.non.zero.genos <- filter(full.df, Genotype %in% non.zero.genos)
hist(log(rowSums(df.non.zero.genos[ ,interaxns_noPont]) +1))
interaction_abund.non.zero.genos <- rowSums(df.non.zero.genos[ ,interaxns_noPont])
library(lme4)
interaction.abund.lmer <- lmer(log(interaction_abund+1) ~ offset(log(shootEst.no18)) + (1|Genotype), full.df)
heritability(interaction.abund.lmer) # ~0.28 with non.zero.genos. 0.33 with full dataset, but the residuals look pretty good actually. Actually, setting an offset and log(x+1) transforming the response variable seems to do a pretty good job of making the random effect normally distributed.
ggQQ_ranef(ranef(interaction.abund.lmer)$Genotype$'(Intercept)')
hist(ranef(interaction.abund.lmer)$Genotype$'(Intercept)')
```

Table 1: Quantitative full data set with euclidean distance
```{r, results='asis', echo=FALSE}
full.adonis <- adonis(log(full.df.quant+1) ~ log(shootEst.no18) + Genotype, 
                      data = full.df, method = "euclidean")
heritability.adonis(full.adonis, full.df)
kable.adonis(full.adonis)
```

Table S2: Qualitative full data set with euclidean distance
```{r, results='asis', echo=FALSE}
full.adonis.qual <- adonis(full.df.qual ~ log(shootEst.no18) + Genotype, 
                           data = full.df, method = "euclidean")
heritability.adonis(full.adonis.qual, full.df)
kable.adonis(full.adonis.qual)
```

Table S3: Qualitative full data set with Raup-crick null model
```{r, results='asis', echo=FALSE}
rc.dist.full <- raupcrick(full.df.qual) # calculate raupcrick dissimilarity for teasing apart differences in alpha vs. beta-diversity
rc.dist.full.adonis <- adonis(rc.dist.full ~ log(shootEst.no18) + Genotype, 
                              data = full.df)
heritability.adonis(rc.dist.full.adonis, full.df)
kable.adonis(rc.dist.full.adonis)
```

Table S3: Qualitative sub data set with Raup-crick null model
```{r, results='asis', echo=FALSE}
rc.dist.sub <- raupcrick(sub.df.qual) # calculate raupcrick dissimilarity for teasing apart differences in alpha vs. beta-diversity
rc.dist.sub.adonis <- adonis(rc.dist.sub ~ log(shootEst.no18) + Genotype, 
                             data = sub.df)
heritability.adonis(rc.dist.sub.adonis, sub.df)
kable.adonis(rc.dist.sub.adonis)
```

Table S2: Quantitative data subset, euclidean
```{r, results='asis', echo=FALSE}
sub.adonis.euclid <- adonis(log(sub.df.quant+1) ~ log(shootEst.no18) + Genotype, 
                            data = sub.df, method = "euclidean")
kable.adonis(sub.adonis.euclid)
```

Table S3: Quantitative data subset, with bray-curtis dissimilarity
```{r, results='asis', echo=FALSE}
sub.adonis.bray <- adonis(sub.df.quant ~ log(shootEst.no18) + Genotype, 
                          data = sub.df, method = "bray")
kable.adonis(sub.adonis.bray)
```

Table S3: Qualitative data subset, with Jaccard dissimilarity
```{r, results='asis', echo=FALSE}
sub.adonis.jac <- adonis(sub.df.qual ~ log(shootEst.no18) + Genotype, 
                         data = sub.df, method = "jaccard")
heritability.adonis(adonis.object = sub.adonis.jac, data = sub.df)
kable.adonis(adonis.object = sub.adonis.jac)
```

Table S4: Qualitative data subset, with euclidean dissimilarity
```{r, results='asis', echo=FALSE}
sub.adonis.euc <- adonis(sub.df.qual ~ log(shootEst.no18) + Genotype, 
                         data = sub.df, method = "euclidean")
kable.adonis(adonis.object = sub.adonis.euc)
```

Heritability of focal gall-parasitoid guild interactions
```{r gall-ptoid heritability, echo=FALSE}
## Variation among genotypes
gall.ptoid.variation <- cbind(Genotype = full.df$Genotype, full.df[ ,c("vLG_egg","vLG_ecto","rG_ecto")]/full.df$shootEst.no18) %>%
  group_by(Genotype) %>%
  select(vLG_egg, vLG_ecto, rG_ecto) %>%
  summarise_each(funs(mean))

max(filter(gall.ptoid.variation, vLG_egg > 0)$vLG_egg)/min(filter(gall.ptoid.variation, vLG_egg > 0)$vLG_egg) # 35-fold variation in vLG_egg parasitism among genotypes that hosted this interaction

max(filter(gall.ptoid.variation, vLG_ecto > 0)$vLG_ecto)/min(filter(gall.ptoid.variation, vLG_ecto > 0)$vLG_ecto) # 20-fold variation in vLG_ecto parasitism among genotypes that hosted this interaction

max(filter(gall.ptoid.variation, rG_ecto > 0)$rG_ecto)/min(filter(gall.ptoid.variation, rG_ecto > 0)$rG_ecto) # 9-fold variation in vLG_ecto parasitism among genotypes that hosted this interaction

## Iteomyia (vLG) - Egg parasitoid 
vLG_egg.lmer.full <- lmer(log(vLG_egg + 1) ~ log(shootEst.no18) + (1|Genotype), 
                          data = full.df)
heritability(vLG_egg.lmer.full)
ggQQ_ranef(ranef(vLG_egg.lmer.full)$Genotype$"(Intercept)") # non-normal random effect

vLG_egg.lmer.sub <- lmer(log(vLG_egg + 1) ~ log(shootEst.no18) + (1|Genotype), 
                         data = sub.df)
heritability(vLG_egg.lmer.sub)
ggQQ_ranef(ranef(vLG_egg.lmer.sub)$Genotype$"(Intercept)") # normal random effect

## Iteomyia (vLG) - Larval parasitoid (ecto)
vLG_ecto.lmer.full <- lmer(log(vLG_ecto + 1) ~ log(shootEst.no18) + (1|Genotype), 
                           data = full.df)
heritability(vLG_ecto.lmer.full)
ggQQ_ranef(ranef(vLG_ecto.lmer.full)$Genotype$"(Intercept)") # close to normal random effect

vLG_ecto.lmer.sub <- lmer(log(vLG_ecto + 1) ~ log(shootEst.no18) + (1|Genotype), 
                          data = sub.df)
heritability(vLG_ecto.lmer.sub) # marginally significant
ggQQ_ranef(ranef(vLG_ecto.lmer.sub)$Genotype$"(Intercept)") # normal random effect

## Rabdophaga salicisbrassicoides (rG) - Larval parasitoid (ecto)
rG_ecto.lmer.full <- lmer(log(rG_ecto + 1) ~ log(shootEst.no18) + (1|Genotype), 
                          data = full.df)
heritability(rG_ecto.lmer.full)
ggQQ_ranef(ranef(rG_ecto.lmer.full)$Genotype$"(Intercept)") # non-normal random effect

rG_ecto.lmer.sub <- lmer(log(rG_ecto + 1) ~ log(shootEst.no18) + (1|Genotype), 
                         data = sub.df)
heritability(rG_ecto.lmer.sub)
ggQQ_ranef(ranef(rG_ecto.lmer.sub)$Genotype$"(Intercept)") # close to normal random effect
```

Gall density and gall size variation among willow genotypes
```{r}
gall.density.variation <- cbind(Genotype = full.df$Genotype, full.df[ ,c("vLG_abund","rG_abund")]/full.df$shootEst.no18) %>%
  group_by(Genotype) %>%
  select(vLG_abund, rG_abund) %>%
  summarise_each(funs(mean))

max(filter(gall.density.variation, vLG_abund > 0)$vLG_abund)/min(filter(gall.density.variation, vLG_abund > 0)$vLG_abund) # 67-fold variation in Iteomyia among willow genotypes

max(filter(gall.density.variation, rG_abund > 0)$rG_abund)/min(filter(gall.density.variation, rG_abund > 0)$rG_abund) # 62-fold variation in Rabdophaga-bud among willow genotypes.

vLG.lmer.full <- lmer(log(vLG_abund + 1) ~ log(shootEst.no18) + (1|Genotype),
                      data = full.df)
heritability(merMod.object = vLG.lmer.full)

rG.lmer.full <- lmer(log(rG_abund + 1) ~ log(shootEst.no18) + (1|Genotype),
                     data = full.df)
heritability(merMod.object = rG.lmer.full)

## vLG size. To test this, I took advantage of the full gall data set and used a nested random effect model.
vLG.size.df <- filter(gall.size.df, gall.sp == "vLG")
table(vLG.size.df$Genotype) # large heterogeneity in sample sizes
table(vLG.size.df$plant.position) # large heterogeneity in sample sizes

vLG.size.lmer <- lmer(gall.height ~ 1 + (1 | Genotype) + (1 |plant.position:Genotype),
                      data = vLG.size.df)
vLG.size.lmer.gen <- update(vLG.size.lmer, .~. - (1 | plant.position:Genotype)) # only Genotype as random effect
vLG.size.lmer.pp <- update(vLG.size.lmer, .~. - (1 | Genotype)) # only plant.position as random effect

vLG.size.lmer.sum <- summary(vLG.size.lmer)
summary(vLG.size.lmer.gen)
summary(vLG.size.lmer.pp)

gen.var <- vLG.size.lmer.sum$varcor$Genotype[1] # Variance due to Genotype
pp.var <- vLG.size.lmer.sum$varcor$'plant.position:Genotype'[1] # variance due to plant position
res.var <- 4.1133 # residual variance - don't know how to extract this from summary object.

H2.vLG.size <- gen.var/(gen.var + pp.var + res.var) # 0.15
exactRLRT(m = vLG.size.lmer.gen, mA = vLG.size.lmer, m0 = vLG.size.lmer.pp) # this model tests whether 

## rG size. To test this, I took advantage of the full gall data set and used a nested random effect model.
rG.size.df <- filter(gall.size.df, gall.sp == "rG")
table(rG.size.df$Genotype) # large heterogeneity in sample sizes
table(rG.size.df$plant.position) # large heterogeneity in sample sizes

rG.size.lmer <- lmer(gall.height ~ 1 + (1 | Genotype) + (1 |plant.position:Genotype),
                     data = rG.size.df)
rG.size.lmer.gen <- update(rG.size.lmer, .~. - (1 | plant.position:Genotype)) # only Genotype as random effect
rG.size.lmer.pp <- update(rG.size.lmer, .~. - (1 | Genotype)) # only plant.position as random effect

rG.size.lmer.sum <- summary(rG.size.lmer)
summary(rG.size.lmer.gen)
summary(rG.size.lmer.pp)

gen.var <- rG.size.lmer.sum$varcor$Genotype[1] # Variance due to Genotype
pp.var <- rG.size.lmer.sum$varcor$'plant.position:Genotype'[1] # variance due to plant position
res.var <- 1.80518 # residual variance - don't know how to extract this from summary object.

H2.rG.size <- gen.var/(gen.var + pp.var + res.var) # 0.0418
exactRLRT(m = rG.size.lmer.gen, mA = rG.size.lmer, m0 = rG.size.lmer.pp) # this model tests whether 
```


RDA analysis of community response to Iteomyia density and size as well as R. salicisbrassicoides density
```{r RDA, echo=FALSE}
full.rda.df <- full.df %>%
  filter(vLG.height.mean > 0) %>% 
  mutate(vLG_density = vLG_abund/shootEst.no18,
         rG_density = rG_abund/shootEst.no18)
full.rda.df.interaxn <- full.rda.df[ ,interaxns_noPont]
full.rda.df.interaxn.qual <- ifelse(full.rda.df[ ,interaxns_noPont] > 0, 1, 0)

full.rda.df.interaxn.log1 <- log(full.rda.df.interaxn+1)
full.rda.df.interaxn.log1.offset <- full.rda.df.interaxn.log1 - log(full.rda.df$shootEst.no18) # this is the same as applying an "offset" to all of the data. In other words, to interpret it in terms of density

rda.full <- rda(full.rda.df.interaxn.log1.offset ~  
                  vLG_density*vLG.height.mean + rG_density, data = full.rda.df)
RsquareAdj(rda.full)
anova(rda.full)
anova(rda.full, by = "margin")

rda.traits.df <- na.omit(data.frame(full.rda.df.interaxn, 
                                    shootEst.no18 = full.rda.df$shootEst.no18,
                                    C_N_imputed = full.rda.df$C_N_imputed, 
                                    Height = full.rda.df$Height, 
                                    sal_tannin.PC1 = full.rda.df$sal_tannin.PC1, 
                                    flavonOLES.PC1 = full.rda.df$flavonOLES.PC1))
rda.full.traits <- rda(log(rda.traits.df[ ,interaxns_noPont]+1) ~ 
                         Condition(log(shootEst.no18)) + 
                         sal_tannin.PC1, data = rda.traits.df)
RsquareAdj(rda.full.traits)
anova(rda.full.traits)
anova(rda.full.traits, by = "margin")
vif.cca(rda.full.traits)
plot(rda.full.traits)

rda.full.qual <- rda(full.rda.df.interaxn.qual ~  
                       vLG_density + vLG.height.mean + rG_density, data = full.rda.df)
RsquareAdj(rda.full.qual)
anova(rda.full.qual)
anova(rda.full.qual, by = "margin")

pos.interactions <- which(rowSums(full.rda.df.interaxn.qual) > 0)
rda.full.qual.jac <- capscale(full.rda.df.interaxn.qual[pos.interactions, ] ~  
                                vLG_density+vLG.height.mean + rG_density, 
                              data = full.rda.df[pos.interactions, ], distance = "jaccard")
RsquareAdj(rda.full.qual.jac)
plot(rda.full.qual.jac)
anova(rda.full.qual.jac)
anova(rda.full.qual.jac, by = "margin")
vif.cca(rda.full.qual.jac)

## vLG-egg parasitoid response
vLG_egg.lm <- lm(log(vLG_egg + 1) ~ offset(log(shootEst.no18)) + 
                   vLG_density*vLG.height.mean,
                 data = full.rda.df)
summary(vLG_egg.lm)
plot(vLG_egg.lm)
library(visreg)
visreg(vLG_egg.lm, xvar = "vLG_density", by = "vLG.height.mean")

## vLG-ecto parasitoid response
vLG_ecto.lm <- lm(log(vLG_ecto+1) ~ offset(log(shootEst.no18)) + vLG_density + vLG.height.mean,
                  data = full.rda.df)
summary(vLG_ecto.lm)
plot(vLG_ecto.lm)

## rG-ecto parasitoid response
rG_ecto.lm <- lm(log(rG_ecto+1) ~ offset(log(shootEst.no18)) + rG_density,
                 data = full.rda.df)
summary(rG_ecto.lm)

```

```{r functional responses}
library(ggplot2)
full.rda.df.small.galls <- filter(full.rda.df, vLG.height.mean < 8)
ggplot(full.rda.df.small.galls, aes(x = vLG_density, y = vLG_egg/vLG_abund)) +
  geom_point(aes(size = vLG_abund)) +
  stat_smooth(method = "glm", family = binomial, aes(weight = full.rda.df.small.galls$vLG_abund))

vLG_egg.bin <- glm(cbind(vLG_egg, vLG_abund - vLG_egg) ~ vLG_density*vLG.height.mean, 
                   data = full.rda.df, family = binomial)
summary(vLG_egg.bin)
anova(update(vLG_egg.bin, .~. - vLG_density*vLG.height.mean), vLG_egg.bin, test = "Chi")

plot(vLG_egg.bin)
visreg(vLG_egg.bin, xvar = "vLG_density", by = "vLG.height.mean", scale = "response")
visreg2d(vLG_egg.bin, x = "vLG.height.mean", y = "vLG_density", scale = "response")

vLG_ecto.bin <- glm(cbind(vLG_ecto, vLG_abund - vLG_ecto) ~ vLG_density + vLG.height.mean,
                    data = full.rda.df, family = binomial)
summary(vLG_ecto.bin)
anova(update(vLG_ecto.bin, .~. - vLG_density - vLG.height.mean), vLG_ecto.bin, test = "Chi")
plot(vLG_ecto.bin)
visreg(vLG_ecto.bin, scale = "response")
visreg2d(vLG_ecto.bin, x = "vLG.height.mean", y = "vLG_density", scale = "response")

rG_ecto.bin <- glm(cbind(rG_ecto, rG_abund - rG_ecto) ~ I(rG_abund/shootEst.no18),
                   data = full.df, family = binomial)
summary(rG_ecto.bin)
anova(update(rG_ecto.bin, .~. - I(rG_abund/shootEst.no18)), rG_ecto.bin, test = "Chi")
plot(rG_ecto.bin)
```

```{r}
# try mantel test
data.for.mantel <- na.omit(select(full.df, plant.position, shootEst.no18,
                                  aSG_Tory:rG_Platy, rG_Tory, SG_Platy, vLG_Eulo:vLG_Tory, # all interactions
                                  Total_Area:HCH.tremulacin.__A220nm, water_content,
                                  specific_leaf_area, C_N_imputed))
traits.stand <- decostand(select(data.for.mantel, Total_Area:C_N_imputed), method = "standardize")
traits.dist <- vegdist(traits.stand, method = "euclidean")
network.stand <- (select(data.for.mantel, aSG_Tory:vLG_Tory)/data.for.mantel$shootEst.no18)
network.dist <- vegdist(network.stand, method = "euclidean")

network.for.dist <- mantel(xdis = traits.dist, ydis = network.dist)

```

Gall density mechanisms
```{r}
traits <- full.df %>%
  select(plant.position,
         Total_Area, Height, Density, # architecture traits
         Trichome.No., specific_leaf_area, water_content, C_N_imputed, # other leaf quality traits
         sal_tannin.PC1, cinn.PC1, cinn.PC2, # salicylates/tannins and phenolic acids
         flavonOLES.PC1, flavonOLES.PC2, flavanonOLES.PC1) %>% # flavonoids
  mutate(log_size = log(Total_Area),
         log_trichomes = log(Trichome.No.+1))
traits <- na.omit(traits)
height_resid <- residuals(lm(Height ~ log_size, traits))
density_resid <- residuals(lm(Density ~ log_size, traits))
sla_resid <- residuals(lm(specific_leaf_area ~ water_content, traits))

traits.df <- cbind.data.frame(select(traits, plant.position, water_content:log_trichomes), 
                              height_resid, density_resid, sla_resid)

trait.names <- colnames(traits.df)[-1]

gall.density.size.df <- full.df %>% 
  select(vLG_abund, rG_abund, shootEst.no18, plant.position, vLG.height.mean) %>%
  mutate(vLG_density = vLG_abund/shootEst.no18,
         rG_density = rG_abund/shootEst.no18) %>%
  select(plant.position:rG_density)

gall.density.size.traits.df <- left_join(gall.density.size.df, traits.df)
vLG.df <- na.omit(select(gall.density.size.traits.df, vLG_density, water_content:sla_resid))
vLG.size.df <- na.omit(select(gall.density.size.traits.df, vLG.height.mean, water_content:sla_resid))
rG.df <- na.omit(select(gall.density.size.traits.df, rG_density, water_content:sla_resid))
#focal.df <- na.omit(gall.density.size.traits.df)
#focal.df.traits <- focal.df[ ,trait.names]

# Iteomyia salicisverruca (vLG) abundance model
vLG.dens.glm <- glm(vLG_abund > 0 ~ offset(log(shootEst.no18)) + C_N_imputed + Total_Area + sal_tannin.PC1, full.df, family = binomial)
summary(vLG.dens.glm)
plot(vLG.dens.glm)

vLG.dens.lm <- lm(log(vLG_abund) ~ offset(log(shootEst.no18)) + C_N_imputed + Total_Area + sal_tannin.PC1, filter(full.df, vLG_abund > 0))
summary(vLG.dens.lm)
plot(vLG.dens.lm)
visreg(vLG.dens.lm)

plot(log(vLG_abund/shootEst.no18) ~ sal_tannin.PC1, filter(full.df, vLG_abund > 0))
library(MuMIn)

dredge(vLG.full)
vLG.null <- lm(vLG.df$vLG_density ~ 1, vLG.df[ ,trait.names])
vLG.full <- lm(vLG.df$vLG_density ~ ., vLG.df[ ,trait.names])
summary(vLG.full) # adj. R2 = 0.03 DOESN'T PASS THE TEST!!!
vLG.step <- step(vLG.null, scope = formula(vLG.full), direction = "forward")
vLG.update <- update(vLG.null, .~. + sal_tannin.PC1 + log_size + height_resid) 
summary(vLG.update)
#plot(vLG.update)

vLG.size.null <- lm(vLG.size.df$vLG.height.mean ~ 1, vLG.size.df[ ,trait.names])
vLG.size.full <- lm(vLG.size.df$vLG.height.mean ~ ., vLG.size.df[ ,trait.names])
summary(vLG.size.full) # adj. R2 = 0.03 DOESN'T PASS THE TEST!!!
vLG.size.step <- step(vLG.size.null, scope = formula(vLG.size.full), direction = "forward")
vLG.size.update <- update(vLG.size.null, .~. + flavonOLES.PC1) 
summary(vLG.size.update)

rG.null <- lm(rG.df$rG_density ~ 1, rG.df[ ,trait.names])
rG.full <- lm(rG.df$rG_density ~ ., rG.df[ ,trait.names])
summary(rG.full) # adj. R2 = 0.12 DOESN'T PASS THE TEST!!!
rG.step <- step(rG.null, scope = formula(rG.full), direction = "forward")
rG.update <- update(rG.null, .~. + log_size  + height_resid + 
                      density_resid) 
summary(rG.update)

```

Analysis of dissimilarity in interaction composition
```{r}
links <- full.df[ ,interaxns_noPont] 
link_abund <- rowSums(links)

links.dist <- cbind.data.frame(Genotype = full.df$Genotype, links) %>%
  filter(link_abund > 0, Genotype %in% c("*","A","B","D","E","F","I","K","L","Q","S","T","V","W","X","Y","Z"))

links.geno.dist <- meandist(vegdist(links.dist[ ,-1], "bray"), links.dist$Genotype)
summary(links.geno.dist) # on average, gall communities were 57% dissimilar from each other
plot(links.geno.dist, cluster = "ward")
```

```{r genotype level}
geno.level.interactions.df <- full.df %>%
  select(Genotype, aSG_Tory:rG_Platy, rG_Tory, SG_Platy, vLG_Eulo:vLG_Tory, aSG_abund:vLG_abund, vLG.height.mean) %>%
  group_by(Genotype) %>%
  summarise_each(funs(mean.na.rm = mean(., na.rm = TRUE))) %>%
  #na.omit() %>%
  mutate(log.vLG_abund = log(vLG_abund),
         log.vLG.height.mean = log(vLG.height.mean),
         log.1.rG_abund = log(rG_abund+1),
         log.1.aSG_abund = log(aSG_abund+1))

library(psych)
corr.test(geno.level.interactions.df[ ,c("vLG_abund","rG_abund","aSG_abund","vLG.height.mean")])

car::scatterplotMatrix(select(geno.level.interactions.df, vLG_abund, rG_abund, aSG_abund, vLG.height.mean))
corr.test(geno.level.interactions.df[ ,c("log.vLG_abund","log.1.rG_abund","log.1.aSG_abund","log.vLG.height.mean")])

# correlations are rather weak
corr.test(select(geno.level.interactions.df, aSG_Tory:vLG_Tory))
corr.test(select(geno.level.interactions.df, vLG_Platy, vLG_Mesopol, vLG_Tory, rG_Tory))
car::scatterplotMatrix(select(geno.level.interactions.df, vLG_Platy, vLG_Mesopol, vLG_Tory, rG_Tory))

plot(rG_Eulo ~ vLG_Mymarid, geno.level.interactions.df)

plot(log.1.aSG_abund ~ log.vLG_abund, geno.level.interactions.df)

geno.level.mvabund <- mvabund(geno.level.interactions.df[ ,interaxns_noPont])

geno.level.manyglm <- manylm(geno.level.mvabund ~ vLG_abund*vLG.height.mean + rG_abund + aSG_abund, data = geno.level.interactions.df)

plot(geno.level.manyglm, which = 1:3)
anova.manylm(geno.level.manyglm, p.uni = "unadjusted")

hist(geno.level.interactions.df$vLG_Platy)
plot(vLG_Mesopol ~ vLG.height.mean, geno.level.interactions.df)
vLG_Platy.lm <- lm(vLG_Mesopol ~ vLG_abund+vLG.height.mean, geno.level.interactions.df)
summary(vLG_Platy.lm)
plot(vLG_Platy.lm)
AIC(vLG_Platy.lm)
library(visreg)
visreg(vLG_Platy.lm, scale = "response")
```

net.trait <- mvabund(full.predictors[ ,interaxns_noPont])

interactions <- full.predictors[ ,interaxns_noPont]
interactions[interactions > 0] <- 1

net.trait.sub <- mvabund(interactions)
# even though manyglm doesn't appear to test interaction terms appropriately (sequential vs. dropping different terms), there is not interactive effect of vLG height and vLG abundance in univariate negative binomial models for any vLG_parasitoid interaction, so I believe this final model is appropriate for the data.
net.mvabund.pres <- manyglm(net.trait.sub ~ log.vLG_abund + vLG.height.mean + log.1.rG_abund + log.1.aSG_abund, data = full.predictors, family = "binomial")
plot(net.mvabund.pres, which = 1:3)
anova.qual <- anova(net.mvabund.pres, p.uni = "unadjusted")
anova.qual

summary(glm(cbind(rG_Tory, rG_Platy + rG_Eulo + rG_Lestodip + rG_Mesopol) ~ vLG.height.mean, data = full.predictors, family = "binomial"))
visreg(glm(cbind(vLG_Tory, vLG_Platy + vLG_Eulo + vLG_Mymarid + vLG_Mesopol) ~ vLG.height.mean + vLG_abund, data = full.predictors, family = "binomial"), scale = "response")



```{r}
link.predict <- predict(manyglm.full, 
                        newdata = data.frame(Genotype = levels(full.df$Genotype)), 
                        type = "response")
link.predict.df <- cbind.data.frame(link.predict, Genotype = levels(full.df$Genotype))
library(vegan)
mean(vegdist(select(link.predict.df, -Genotype), method = "horn"))
sd(vegdist(select(link.predict.df, -Genotype), method = "horn"))
hist(vegdist(select(link.predict.df, -Genotype), method = "horn"))

link.predict.melt <- melt(link.predict.df)

ggplot(link.predict.melt, aes(x = Genotype, y = value, fill = variable)) + geom_bar(stat = "identity")
```

```{r}
full.mean <- cbind.data.frame(Genotype = full.df$Genotype, full.mvabund) %>%
  group_by(Genotype) %>%
  summarise_each(funs(mean.narm = mean(., na.rm = TRUE)))

# calculate average dissimilarity in link community composition among willow genotypes using average link abundance data.
require(vegan)
full.mean.df.noU <- full.mean %>% filter(Genotype != "U") # genotype U never had any interactions
full.mean.bray <- vegdist(select(full.mean.df.noU, -Genotype), method = "bray")
mean(full.mean.bray) # on average, link composition was 77% dissimilar from each other.
sd(full.mean.bray)
range(full.mean.bray)

plot(hclust(full.mean.bray, method = "average"), labels = full.mean.df.noU$Genotype)
#NMDS <- metaMDS(comm = select(full.mean.df.noU, -Genotype), distance = "bray")
#plot(NMDS, type = "t")

require(tidyr)
beta_partition_lists <- as.data.frame(full.mean.df.noU) %>%
  gather(key = Genotype, value = link_abund, aSG_Tory:vLG_Tory) %>% # why don't variable and value columns change?
  #mutate(variable = as.character(variable)) %>%
  separate(variable, into = c("Gall","Parasitoid"), sep = "_") %>%
  filter(value > 0) %>%
  cast(Gall ~ Parasitoid | Genotype, sum) # note that sum function doesn't do anything. I'm just reshaping the dataframe into multiple lists.
beta_partition_lists <- beta_partition_lists[-21]

for(i in 1:length(beta_partition_lists)) {
  rownames(beta_partition_lists[[i]]) <- as.character(beta_partition_lists[[i]]$Gall) # change row names to gall specie names
  beta_partition_lists[[i]] <- as.data.frame(select(beta_partition_lists[[i]], -Gall)) # remove gall.sp as a column of data. Turning this into a matrix is important for easier analysis with betalink package
  beta_partition_lists[[i]] <- as.matrix(beta_partition_lists[[i]])  # need to be converted into matrices for "as.table" function to work for quantitative data.
}

### Beta-diversity of interaction network analysis
source('~/Documents/betalink/R/betalink.R')
source('~/Documents/betalink/R/vec2data.frame.R')
source('~/Documents/betalink/R/measures.r')
source('~/Documents/betalink/R/betalink.b.r')
source('~/Documents/betalink/R/betalink.dist.r')
source('~/Documents/betalink/R/betalink.q.R')

# Betalink Quantitative
betalink_quantitative <- betalink.dist(beta_partition_lists, bf = "bray", triangular = T) #
betalink_qualitative <- betalink.dist(beta_partition_lists, bf = B01, triangular = T)

mean(betalink_quantitative$WN) # 77% dissimilarity with bray
mean(betalink_qualitative$WN) # 63% dissimilarity with whittaker

mean(betalink_quantitative$OS) # 33% with bray
mean(betalink_qualitative$OS) # 6% with whittaker

mean(betalink_quantitative$ST) # 44% with bray
mean(betalink_qualitative$ST) # 57% with whittaker

mean(betalink_quantitative$contrib) # species turnover contributes to 53% of the dissimilarity in interaction networks, according to bray-curtis
mean(betalink_qualitative$contrib) # species turnover contributes to 86% of the dissimilarity in interaction networks, according to whittaker index on qualitative data.
```

# subset of full data set. Only contains Genotype with 3+ replicates in the subset of data that only contains willows with at least 1 gall-parasitoid interaction found on them.
abund_interaxns_noPont <- rowSums(full.df[ ,interaxns_noPont]) # used to filter willow replicates with at least 1 gall-parasitoid interaction found on them.
table(filter(full.df, abund_interaxns_noPont > 0)$Genotype) # number of replicates for all genotypes in reduced dataset.
geno_3plus <- c("*","B","D","I","K","L","Q","S","V","W","X","Y","Z") # Genotypes with 3+ replicates
sub.df <- filter(full.df, abund_interaxns_noPont > 0, Genotype %in% geno_3plus)

# quantitative interaction data
full.df.quant <- full.df[ ,interaxns_noPont]
sub.df.quant <- sub.df[ ,interaxns_noPont]

# qualitative interaction data
full.df.qual <- ifelse(full.df[ ,interaxns_noPont] > 0, 1, 0)
sub.df.qual <- ifelse(sub.df[ ,interaxns_noPont] > 0, 1, 0)

library(lme4)
library(RLRsim)
plot(vLG_abund ~ Genotype, full.df)
vLG.lmer <- lmer(log(vLG_abund+1) ~ (1|Genotype), full.df)
heritability(vLG.lmer)

with(full.df, table(vLG_abund, Genotype))
with(full.df, table(rG_abund, Genotype))
hist(full.df$rG_abund)
hist(log(full.df$rG_abund+1))

library(glmmADMB)
hist(subset(full.df, vLG_abund > 0)$vLG_abund)
fit1 <- glmmadmb(vLG_abund ~ (1|Genotype), subset(full.df, vLG_abund > 0), family="truncpoiss" )
summary(fit1)
t <- ranef(fit1)$Genotype
ggQQ_ranef(t)

fit2 <- glmmadmb(vLG_abund >0 ~ (1|Genotype), subset(full.df, vLG_abund > 0), family="binomial" )
summary(fit1)

vLG.glmer <- glmer(vLG_abund>0 ~ (1|Genotype), full.df, family = binomial)
summary(vLG.glmer)
plot(vLG.glmer)
overdisp_fun(vLG.glmer) # individual plant effect cleans up overdispersion.
ggQQ_ranef(ranef(vLG.glmer)$Genotype$'(Intercept)') # looks pretty good
confint.merMod(vLG.glmer, method = "profile")
1.479/(1.479+(pi^2)/3) # 0.29 H2

plot(vLG_Platy ~ Genotype, full.df)
vLG.Platy.lmer <- lmer(log(vLG_Platy+1) ~ (1|Genotype), full.df)
summary(vLG.Platy.lmer)
heritability(vLG.Platy.lmer)

with(full.df, table(vLG_Platy, Genotype))
with(full.df, table(vLG_Mesopol, Genotype))
with(full.df, table(vLG_Tory, Genotype))
vLG.Platy.glmer <- glmer(vLG_Tory > 0 ~ (1|Genotype), full.df, family = binomial)
summary(vLG.Platy.glmer)
overdisp_fun(vLG.Platy.glmer) 
plot(vLG.Platy.glmer)
ggQQ_ranef(ranef(vLG.Platy.glmer)$Genotype$'(Intercept)') # looks pretty good
ggQQ_ranef(ranef(vLG.Platy.glmer)$plant.position$'(Intercept)') # pretty skewed near middle...
confint.merMod(vLG.Platy.glmer, method = "profile")
2.962/(2.962+1+log(1/exp(-1.8822)+1))

library(MCMCglmm)
rpt.poisGLMM.add(y = full.df$vLG_Platy, groups = full.df$Genotype)
vLG_Platy.MCMC <- MCMCglmm(fixed = vLG_Platy > 0 ~ 1, random = ~ Genotype, family = "categorical", data = full.df)
summary(vLG_Platy.MCMC)
allChains <- as.mcmc(cbind(vLG_Platy.MCMC$Sol, vLG_Platy.MCMC$VCV))
xyplot(allChains)

rpt.binomGLMM.add(y = full.df$vLG_Platy > 0, groups = full.df$Genotype)

library(mvabund)

coef.df.3 <- mutate(melt(coef(net.mvabund.3)),
                    predictor_response = paste(X1, X2, sep = "_")) %>%
  select(predictor_response, value)
coef.df.4 <- mutate(melt(coef(net.mvabund.4)),
                    predictor_response = paste(X1, X2, sep = "_")) %>%
  select(predictor_response, value)

coef.df.all <- join_all(list(coef.df.2, coef.df.3, coef.df.4), by = "predictor_response")