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