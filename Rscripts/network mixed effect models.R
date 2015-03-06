## network mixed-effect models

## source in required data
source('~/Documents/Genotype_Networks/Rscripts/network_management_gall_level.R')

## source in required libraries
library(lme4)
library(visreg)
library(gamm4)

gall_mech.df

## testing for overdispersion
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

## vLG mechanisms df
vLG_mech.df <- filter(gall_mech.df, vLG_total > 0 & gall.height > 0) %>%
  mutate(sc.gall.height = scale(gall.height, center = TRUE, scale = TRUE),
         sc.vLG_density = scale(vLG_density, center = TRUE, scale = TRUE),
         sc.vLG_abund = scale(vLG_abund, center = TRUE, scale = TRUE),
         vLG_egg = vLG_Platy + vLG_Mymarid,
         vLG_ecto = vLG_Eulo + vLG_Mesopol + vLG_Tory,
         vLG_parasitized = vLG_egg + vLG_ecto,
         gall.id.nest = factor(gall.id.nest))
str(vLG_mech.df)

## vLG_parasitized
library(MCMCglmm)
t <- rpt.binomGLMM.add(y = with(vLG_mech.df, cbind(vLG_parasitized, vLG_total-vLG_parasitized)), groups = vLG_mech.df$Genotype.x)
t
library(embd)
q <- rpt.binomGLMM.multi(y = with(vLG_mech.df, cbind(vLG_parasitized, vLG_total-vLG_parasitized)), groups = vLG_mech.df$Genotype.x)
test <- glmmPQL(cbind(vLG_parasitized, vLG_total-vLG_parasitized) ~ 1,
                random = ~ 1|Genotype.x/plant.position,
                data = vLG_mech.df,
                family = binomial)
summary(test)
0.695707^2/(0.695707^2 + 0.3712638^2 + 1.108155^2)

vLG_ptized.glmer <- glmer(cbind(vLG_parasitized, vLG_total-vLG_parasitized) ~
                            sc.vLG_abund + sc.gall.height + (1|Genotype.x/plant.position),
                          data = vLG_mech.df,
                          family = binomial)
summary(vLG_ptized.glmer)
confint.merMod(vLG_ptized.glmer, method = "profile")
overdisp_fun(vLG_ptized.glmer)
ggQQ_ranef(ranef(vLG_ptized.glmer)$Genotype.x$'(Intercept)') # looks great
ggQQ_ranef(ranef(vLG_ptized.glmer)$plant.position$'(Intercept)') # looks great
profile(vLG_ptized.glmer, which = 1)
AIC(vLG_ptized.glmer)

test <- ranef(vLG_ptized.glmer, condVar = TRUE, whichel = "Genotype.x")
dotplot(test)
test2 <- ranef(vLG_ptized.glmer, condVar = TRUE, whichel = "plant.position:Genotype.x")
dotplot(test2)
qqmath(test2)
qqmath(test)
library(arm)
t <- REsim(vLG_ptized.glmer, whichel = "Genotype.x", nsims = 1000)
plotREsim(t, scale = 1, var = "mean", sd = "sd")

vLG_ptized.glm <- glm(cbind(vLG_parasitized, vLG_total - vLG_parasitized) ~
                        vLG_abund*gall.height,
                      data = vLG_mech.df,
                      family = binomial)
summary(vLG_ptized.glm)
plot(vLG_ptized.glm)
plot(residuals(vLG_ptized.glm) ~ Genotype.x, vLG_mech.df)

vLG_ptized.glmer <- glmer(cbind(vLG_parasitized, vLG_total - vLG_parasitized) ~
                           sc.gall.height + (1|Genotype.x/plant.position),
                          data = vLG_mech.df,
                          family = binomial)
summary(vLG_ptized.glmer)
AIC(vLG_ptized.glmer)
overdisp_fun(vLG_ptized.glmer)
visreg(vLG_ptized.glmer, scale = "response")


library(mgcv)
vLG_ptized.gamm <- gamm4(cbind(vLG_parasitized, vLG_total - vLG_parasitized) ~
                           s(gall.height),
                         random = ~ (1|Genotype.x/plant.position),
                          data = vLG_mech.df,
                          family = binomial)
summary(vLG_ptized.gamm$gam)
summary(vLG_ptized.gamm$mer)
plot(vLG_ptized.gamm$mer)
AIC(vLG_ptized.gamm$mer)
vis.gam(vLG_ptized.gamm$gam, type = "response")

## vLG_Platy
vLG_Platy.glmer <- glmer(cbind(vLG_Platy, vLG_total - vLG_Platy) ~
                            sc.gall.height + (1|Genotype.x/plant.position),
                          data = vLG_mech.df,
                          family = binomial)
summary(vLG_Platy.glmer)
AIC(vLG_Platy.glmer)
overdisp_fun(vLG_Platy.glmer)
visreg(vLG_Platy.glmer, scale = "response")

vLG_Platy.gamm <- gamm4(cbind(vLG_Platy, vLG_total - vLG_Platy) ~
                          s(gall.height) + vLG_abund,
                        random = ~ (1|Genotype.x/plant.position),
                        data = vLG_mech.df,
                        family = binomial)
summary(vLG_Platy.gamm$gam)
summary(vLG_Platy.gamm$mer)
plot(vLG_Platy.gamm$mer)
AIC(vLG_Platy.gamm$mer)
vis.gam(vLG_Platy.gamm$gam, type = "response")

## vLG_egg
vLG_egg.glmer <- glmer(cbind(vLG_egg, vLG_total - vLG_egg) ~ sc.vLG_density*sc.gall.height + 
                         (1 | plant.position:Genotype.x) + (1 | Genotype.x), data = vLG_mech.df,
                       family = binomial(link = "logit"))
summary(vLG_egg.glmer)
AIC(vLG_egg.glmer)
overdisp_fun(vLG_egg.glmer)

vLG_egg.glmer <- glmer(cbind(vLG_egg, vLG_total - vLG_egg) ~ sc.vLG_density + sc.gall.height + sc.vLG_density:sc.gall.height +
                         (1 | Genotype.x/plant.position), 
                       data = vLG_mech.df,
                       family = binomial(link = "logit"),
                       glmerControl(optimizer="bobyqa"))
summary(vLG_egg.glmer)
AIC(vLG_egg.glmer)
vLG_egg.glmer.noInteraction <- update(vLG_egg.glmer, .~. -sc.vLG_density:sc.gall.height)
library(pbkrtest)
PBmodcomp(largeModel = vLG_egg.glmer, smallModel = vLG_egg.glmer.noInteraction)

vLG_egg.glmer.vis <- glmer(cbind(vLG_egg, vLG_total - vLG_egg) ~ vLG_density*gall.height + 
                         (1 | plant.position), data = vLG_mech.df,
                       family = binomial(link = "logit"))
summary(vLG_egg.glmer.vis) # same overall picture as standardized model that includes Genotype as a random effect.
AIC(vLG_egg.glmer.vis)
range(predict(vLG_egg.glmer.vis, type = "response", newdata = filter(vLG_mech.df, gall.height < 8)))
0.886/0.085

## vLG_ecto
vLG_ecto.glmer <- glmer(cbind(vLG_ecto, vLG_total - vLG_ecto) ~ sc.vLG_density + sc.gall.height +
                         (1 | Genotype.x/plant.position), data = vLG_mech.df, 
                        family = binomial(link = "logit"))
summary(vLG_ecto.glmer)
AIC(vLG_ecto.glmer)

vLG_ecto.glmer.vis <- glmer(cbind(vLG_ecto, vLG_total - vLG_ecto) ~ vLG_density + gall.height + 
                             (1 | plant.position), data = vLG_mech.df,
                           family = binomial(link = "logit"))
summary(vLG_ecto.glmer.vis) # same overall picture as standardized model that includes Genotype as a random effect.
AIC(vLG_ecto.glmer.vis)
range(predict(vLG_ecto.glmer.vis, type = "response", newdata = filter(vLG_mech.df, gall.height < 8)))

## rG mechanisms df
rG_mech.df <- filter(gall_mech.df, rG_total > 0) %>%
  mutate(sc.gall.height = scale(gall.height, center = TRUE, scale = TRUE),
         sc.rG_density = scale(rG_density, center = TRUE, scale = TRUE),
         rG_ecto = rG_Tory + rG_Mesopol + rG_Eulo)

## rG_ecto
rG_ecto.glmer <- glmer(cbind(rG_ecto, rG_total - rG_ecto) ~ 1 + 
                          (1 | Genotype.x/plant.position), data = rG_mech.df, family = binomial(link = "logit"))
summary(rG_ecto.glmer)
AIC(rG_ecto.glmer)

## plots for manuscript. Consider theta = 45
grey.scale <- colorRampPalette(colors = c("white","grey","black"))

#pdf(file="~/Documents/Genotype_Networks/iteomyia egg parasitoid mixed effect.png",width=11, height=8.5)
visreg2d(vLG_egg.glmer.vis, y = "vLG_density", x = "gall.height", 
         scale = "response", color.palette = grey.scale, 
         zlim = c(0,1), main = "Probability of Iteomyia - Egg parasitoid interaction",
         xlab = "Gall diameter (mm)", ylab = "Gall density (count per 100 shoots)")
#dev.off()
visreg2d(vLG_ecto.glmer.vis, y = "vLG_density", x = "gall.height", 
         scale = "response", color.palette = grey.scale, 
         zlim = c(0,1), main = "Probability of Iteomyia - Larval parasitoid interaction",
         xlab = "Gall diameter (mm)", ylab = "Gall density (count per 100 shoots)")


