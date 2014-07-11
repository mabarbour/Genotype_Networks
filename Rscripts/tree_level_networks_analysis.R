######## Tree level Network, Community, and Species level analysis

#### upload necessary packages

# data manipulation. order of libraries is important. always load dplyr last.
library(reshape2)
library(reshape)
library(dplyr)

# plotting
library(ggplot2)

# analysis
library(mvabund)
library(vegan)

# upload molten gall network data
gall_net_melt <- read.csv("~/Documents/Genotype_Networks/data/gall_network_data.csv")
gall_net_melt <- tbl_df(gall_net_melt)
gall_net_melt <- filter(gall_net_melt, plant.position != 9999) # removed unknown plant.position

### NEED TO ADD PLANTS WITH NO GALLS

##### Does the proportion of parasitized vLG galls vary among genotypes?

### contingency table analyses for each gall species separately. 

# Iteomyia salicisverruca
vLG_geno_parasitism <- gall_net_melt %>%
  dcast(Genotype ~ gall.sp + gall_contents, sum) %>%
  select(Genotype, starts_with("vLG")) %>%
  mutate(vLG_ptoid_attack = vLG_Eulo.fem + vLG_Eulo.mal + vLG_Mesopol + vLG_Mymarid + vLG_Platy + vLG_Ptero.2 + vLG_Tory.fem + vLG_Tory.mal, vLG_prop_ptoid_attack = vLG_ptoid_attack/(vLG_vLG.pupa + vLG_ptoid_attack)) %>% # omitting exit hole because of potential overlap with other ectoparasitoids. Need to probably resolve this by not using "unk" gall.id
  select(Genotype, vLG_ptoid_attack, vLG_vLG.pupa, vLG_prop_ptoid_attack)

vLG_for_chi <- vLG_geno_parasitism %>%
  filter(vLG_ptoid_attack > 0 & vLG_vLG.pupa > 0) %>% # remove all observations with zeros in either column
  select(vLG_ptoid_attack, vLG_vLG.pupa)

chisq.test(vLG_for_chi)
chisq.test(vLG_for_chi, simulate.p.value = TRUE, B = 10000) # same results as above

# Rabdophaga salicisbrassicoides
rG_geno_parasitism <- gall_net_melt %>%
  dcast(Genotype ~ gall.sp + gall_contents, sum) %>%
  select(Genotype, starts_with("rG")) %>%
  mutate(rG_ptoid_attack = rG_diff.or.larv + rG_Eulo.fem + rG_Eulo.mal + rG_Lestodip + rG_Mesopol + rG_Platy + rG_Tory.fem + rG_Tory.mal + rG_unk.ptoid, rG_prop_ptoid_attack = rG_ptoid_attack/(rG_rG.larv + rG_ptoid_attack)) %>% # omitting exit hole because of potential overlap with other ectoparasitoids. Need to probably resolve this by not using "unk" gall.id
  select(Genotype, rG_ptoid_attack, rG_rG.larv, rG_prop_ptoid_attack)

rG_for_chi <- rG_geno_parasitism %>%
  filter(rG_ptoid_attack > 0 & rG_rG.larv > 0) %>% # remove all observations with zeros in either column
  select(rG_ptoid_attack, rG_rG.larv)

chisq.test(rG_for_chi) # qualitatively the same even if I only remove observations with > 0 in either ptoid attack or rG.larva
chisq.test(rG_for_chi, simulate.p.value = TRUE, B = 10000)

# Rabdophaga salicisbattatus
SG_geno_parasitism <- gall_net_melt %>%
  dcast(Genotype ~ gall.sp + gall_contents, sum) %>%
  select(Genotype, starts_with("SG")) %>%
  mutate(SG_prop_ptoid_attack = SG_Platy/(SG_Platy + SG_SG.larv)) %>% # omitting exit hole because of potential overlap with other ectoparasitoids. Need to probably resolve this by not using "unk" gall.id. SG_Platy was the only species reared from SG
  select(Genotype, SG_Platy, SG_SG.larv, SG_prop_ptoid_attack)

SG_for_chi <- SG_geno_parasitism %>%
  filter(SG_Platy > 0 & SG_SG.larv > 0) %>% # unable to run the anlaysis unless I permit there to be zeros in at least one column
  select(SG_Platy, SG_SG.larv)

chisq.test(SG_for_chi) # I don't know how much I trust this...appears to be driven by one data point.
chisq.test(SG_for_chi, simulate.p.value = TRUE, B = 10000)

# Cecidomyiidae sp. A (aSG)
aSG_geno_parasitism <- gall_net_melt %>%
  dcast(Genotype ~ gall.sp + gall_contents, sum) %>%
  select(Genotype, starts_with("aSG")) %>%
  mutate(aSG_ptoid_attack = aSG_Tory.fem + aSG_Tory.mal, aSG_prop_ptoid_attack = aSG_ptoid_attack/(aSG_ptoid_attack + aSG_aSG.larv)) %>% 
  select(Genotype, aSG_ptoid_attack, aSG_aSG.larv, aSG_prop_ptoid_attack)

aSG_for_chi <- aSG_geno_parasitism %>%
  filter(aSG_ptoid_attack > 0 & aSG_aSG.larv > 0) %>% # 
  select(aSG_ptoid_attack, aSG_aSG.larv)

chisq.test(aSG_for_chi) # qualitatively the same results no matter how I decide to retain the data.
chisq.test(aSG_for_chi, simulate.p.value = TRUE, B = 10000)

# Pontania californica
rsLG_geno_parasitism <- gall_net_melt %>%
  dcast(Genotype ~ gall.sp + gall_contents, sum) %>%
  select(Genotype, starts_with("rsLG")) %>%
  mutate(rsLG_surv = rsLG_Pont.ad + rsLG_Pont.prep, rsLG_ptoid_attack = rsLG_Eury.fem + rsLG_Eury.mal + rsLG_Lathro.fem + rsLG_Lathro.mal, rsLG_prop_ptoid_attack = rsLG_ptoid_attack/(rsLG_surv + rsLG_ptoid_attack)) %>%
  select(rsLG_ptoid_attack, rsLG_surv, rsLG_prop_ptoid_attack)

rsLG_for_chi <- rsLG_geno_parasitism %>%
  filter(rsLG_ptoid_attack > 0 & rsLG_surv > 0) %>% # 
  select(rsLG_ptoid_attack, rsLG_surv)

chisq.test(rsLG_for_chi) # qualitatively the same results no matter how I decide to retain the data.
chisq.test(rsLG_for_chi, simulate.p.value = TRUE, B = 10000)

##### How dissimilar are gall-parasitoid interaction networks among genotypes?

### Genotype-level networks
levels(gall_net_melt$gall_contents)

genotype_net_filter <- gall_net_melt %>%
  filter(gall_contents %in% c("Eulo.fem", "Eulo.mal", "Eury.fem", "Eury.mal", "Lathro.fem", "Lathro.mal", "Lestodip", "Mesopol", "Mymarid", "Platy", "Ptero.2", "Tory.fem", "Tory.mal"))

genotype_net_add <- mutate(genotype_net_filter, gall_contents_collapse = revalue(gall_contents, c("Eulo.fem" = "Eulo", "Eulo.mal" = "Eulo", "Eury.fem" = "Eury", "Eury.mal" = "Eury", "Lathro.fem" = "Lathro", "Lathro.mal" = "Lathro", "Ptero.2" = "Mesopol", "Tory.fem" = "Tory", "Tory.mal" = "Tory") ))

total_net <- cast(genotype_net_add, gall.sp ~ gall_contents_collapse, sum)
total_net_gephi <- melt(total_net)
total_net_gephi <- select(total_net_gephi, Source = gall.sp, Target = gall_contents_collapse, Weight = value)
write.csv(total_net_gephi, "~/Documents/Genotype_Networks/data/total_gall_parasitoid_network.csv")

genotype_net <- cast(genotype_net_add, gall.sp ~ gall_contents_collapse | Genotype, sum)

genotype_net <- genotype_net[-c(2,5,7,8,15,16,20,21)] # betalink doesn't like it if you only have 1 column. Need to figure this out...

# this loop prepares each network for analysis in betalink
for(i in 1:length(genotype_net)) {
  rownames(genotype_net[[i]]) <- as.character(genotype_net[[i]]$gall.sp) # change row names to gall specie names
  genotype_net[[i]] <- as.matrix(select(genotype_net[[i]], -gall.sp)) # remove gall.sp as a column of data. Turning this into a matrix is important for easier analysis with betalink package
}

genotype_net_qualitative <- list()
for(i in 1:length(genotype_net)) {
  genotype_net_qualitative[[i]] <- genotype_net[[i]]
  genotype_net_qualitative[[i]][genotype_net_qualitative[[i]] > 0] <- 1
}
names(genotype_net_qualitative) <- names(genotype_net)

library(betalink)
betalink_genotype_net_qualitative <- betalink.dist(genotype_net_qualitative, bf=B01) # tried using Sorenson's index (B11), but it was giving very wonky results. Also giving wonky results when I try to use Paul Rabie's qualitative code...

mean(betalink_genotype_net_qualitative$WN) # average dissimilarity is 49% among genotypes
mean(betalink_genotype_net_qualitative$OS)
mean(betalink_genotype_net_qualitative$S)
mean(betalink_genotype_net_qualitative$ST)
mean(betalink_genotype_net_qualitative$contrib) # 71% of the dissimilarity is due to species turnover.

betalink.plot(as.matrix(genotype_net_qualitative[["S"]]), as.matrix(genotype_net_qualitative[["I"]]))

betalink_genotype_net_quantitative <- betalink.dist(genotype_net, bf = "bray")

w1 = genotype_net[["S"]]
w2 = genotype_net[["I"]]

mean(betalink_genotype_net_quantitative$WN) # 68% dissimilarity in interaction networks
mean(betalink_genotype_net_quantitative$OS)
mean(betalink_genotype_net_quantitative$S_all.taxa)
mean(betalink_genotype_net_quantitative$S_insects)
mean(betalink_genotype_net_quantitative$S_plants)
mean(betalink_genotype_net_quantitative$ST)
mean(betalink_genotype_net_quantitative$contrib) # species turnover contributes less to dissimilarity in interaction networks than qualitative data (40%)







# vLG proportion analysis
vLG_parasitism <- gall_net_melt %>%
  dcast(Gender + Genotype + plant.position ~ gall.sp + gall_contents, sum) %>%
  mutate(vLG_ptoid_attack = vLG_Eulo.fem + vLG_Eulo.mal + vLG_exit.hole + vLG_Mesopol + vLG_Mymarid + vLG_Platy + vLG_Ptero.2 + vLG_Tory.fem + vLG_Tory.mal) %>%
  select(Gender, Genotype, plant.position, vLG_vLG.pupa, vLG_ptoid_attack) %>%
  filter(vLG_vLG.pupa > 0 | vLG_ptoid_attack > 0)

table(vLG_parasitism$Genotype)

vLG_parasitism_sub2 = filter(vLG_parasitism, Genotype %in% c("*","B","D","E","F","H","I","J","K","L","O","S","V","X","Y","Z") )

vLG_parasitism_glm <- glm(cbind(vLG_ptoid_attack, vLG_vLG.pupa) ~ Genotype, data = vLG_parasitism, family = "binomial")
summary(vLG_parasitism_glm)
anova(vLG_parasitism_glm, test = "Chi")
plot(vLG_parasitism_glm) # need to see what is going on here


# create a tree level data
tree_gall_net <- dcast(gall_net_melt, Gender + Genotype + plant.position ~ gall.sp + gall_contents, sum) 
colSums(select(tree_gall_net, -c(Gender, Genotype, plant.position))) # take a look the abundance of these different "interactions"

focal = select(tree_gall_net, -c(Gender, Genotype, plant.position))
adonis(focal ~ Genotype, data = tree_gall_net, method = "bray")

test = capscale(focal ~ Genotype, data = tree_gall_net, distance = "bray")
plot(test, display = "sp")

focal_mvabund <- mvabund(focal)
test2 <- manylm(focal_mvabund ~ Genotype, data = tree_gall_net) 
anova(test2, test = "F", p.uni = "unadjusted")

plot(vLG_Platy ~ Genotype, tree_gall_net)
plot(vLG_Tory.fem ~ Genotype, tree_gall_net)

#### EVERYTHING BELOW THIS IS OLD BUT MAY BE USEFUL

#library(reshape2)

# plotting
#library(ggplot2)
#library(car) # for scatterplot matrices

# network analysis
#library(bipartite)

### source in required functions
#source('~/Documents/Genotype_Networks/data/fxns_networks.R')

### Upload necessary data

# upload plant position info
plant.position.info <- read.csv("~/Documents/Genotype_Willow_Community/datasets_&_Rscripts/Willow Garden Positions.csv")
plant.position_df <- tbl_df(plant.position.info) %>%
  mutate(plant.position = Plant.Position) %>%
  select(plant.position, Genotype, Gender)

# upload complete network data
gall_network_data <- read.csv("~/Documents/Genotype_Networks/data/gall_network_data.csv")
gall_network_df <- tbl_df(gall_network_data)

# Shoot count estimates as well as plant positions where zero galls were collected
shoot.countEst <- read.csv("~/Documents/Genotype_Networks/data/survey_2_stem_diams_shootEsts.csv")
shoot.count_df <- tbl_df(shoot.countEst)
shoot.count_df <- shoot.count_df %>%
  mutate(plant.position = Plant.Position) %>%
  select(plant.position, Galls.found, shootEst.all, shootEst.no18)

### Create dataset for gall-parasitoid interactions that includes intraspecific variation within the parasitoids
tree_gallnet_intraptoid <- gall_network_df %>%
  select(plant.position, gall.sp, Eury.fem:Eulo.fem, Lestodip:Lathro.mal, Eulo.mal, Mymarid, Genotype, Gender) %>%
  melt(id.vars = c("plant.position","Genotype","Gender","gall.sp"), variable_name = "gall_contents") %>%
  dcast(Gender + Genotype + plant.position ~ gall.sp + gall_contents, sum) %>%
  select(Gender, Genotype, plant.position, aSG_Tory.mal, aSG_Tory.fem, rG_Platy, rG_Tory.mal, rG_Tory.fem, rG_Mesopol, rG_Eulo.fem, rG_Lestodip, rG_Eulo.mal, rsLG_Eury.fem, rsLG_Eury.mal, rsLG_Lathro.fem, rsLG_Lathro.mal, SG_Platy, vLG_Platy, vLG_Tory.mal, vLG_Tory.fem, vLG_Mesopol, vLG_Ptero.2, vLG_Eulo.fem, vLG_Eulo.mal, vLG_Mymarid) # only retain gall-ptoid interactions that had at least 1 observation 

#tree_gallnet_intraptoid <- inner_join(tree_gallnet_intraptoid, shoot.count_df) # retains all data from shoot.count_df which contains relevant information for when no galls were collected from plants. NEED TO DECIDE THE TYPE OF JOIN REQUIRED.
not_in_shootcount <- anti_join(tree_gallnet_intraptoid, shoot.count_df) 
not_in_gallnet <- anti_join(shoot.count_df, tree_gallnet_intraptoid)

with(tree_gallnet_intraptoid, cbind.data.frame(plant.position, Galls.found, vLG_Platy, shootEst.all)) # 0 galls found and NA make sense. But other categories do not, so I need to figure out what is going on here.

library("mvabund")
test = filter(tree_gallnet_intraptoid, shootEst.all > 0)
test.offset = matrix(rep(test$shootEst.all, 22), ncol = 22)
gallnet_mvabund = mvabund(select(test, aSG_Tory.mal:vLG_Mymarid))
gallnet_manylm = manylm(gallnet_mvabund ~ Genotype, data = test, offset = test.offset)
gallnet_anova.manylm = anova.manylm(gallnet_manylm, p.uni = "adjusted")
summary.manylm(gallnet_manylm, p.uni = "adjusted")

#library("MASS") # note that this masks "select" from dplyr
t1 = glm.nb(vLG_Platy ~ Genotype + offset(log(shootEst.all)), test, control = glm.control(maxit = 100))
t2 = glm.nb(vLG_Platy ~ offset(log(shootEst.all)), test, control = glm.control(maxit = 100))
anova(t2,t1)
plot(vLG_Platy ~ Genotype, test) # need to figure out how to plot with offset.

### Note that including the offset changes the interpretation to "density". Including it as a covariate allows it to change at a slower or faster rate. I think density is the original context in which I was doing this. However, the AIC is MUCH lower if I include density as a covariate. had to up the iterations to get rid of the "alternation" error.


### Create gall-parasitoid network at the Genotype level. 
geno_gallnet <- gall_network_df %>%
  select(plant.position, gall.sp, Eury.fem:Eulo.fem, Lestodip:Lathro.mal, Eulo.mal, Mymarid, Genotype, Gender) %>%
  melt(id.vars = c("plant.position","Genotype","Gender","gall.sp"), variable_name = "gall_contents") %>%
  dcast(Gender + Genotype ~ gall.sp + gall_contents, sum) %>%
  select(Gender, Genotype, aSG_Tory.mal, aSG_Tory.fem, rG_Platy, rG_Tory.mal, rG_Tory.fem, rG_Mesopol, rG_Eulo.fem, rG_Lestodip, rG_Eulo.mal, rsLG_Eury.fem, rsLG_Eury.mal, rsLG_Lathro.fem, rsLG_Lathro.mal, SG_Platy, vLG_Platy, vLG_Tory.mal, vLG_Tory.fem, vLG_Mesopol, vLG_Ptero.2, vLG_Eulo.fem, vLG_Eulo.mal, vLG_Mymarid) %>% # only retain gall-ptoid interactions that had at least 1 observation
  mutate(aSG_Tory = aSG_Tory.mal + aSG_Tory.fem, rG_Tory = rG_Tory.mal + rG_Tory.fem, rG_Eulo = rG_Eulo.fem + rG_Eulo.mal, rsLG_Eury = rsLG_Eury.fem + rsLG_Eury.mal, rsLG_Lathro = rsLG_Lathro.fem + rsLG_Lathro.mal, vLG_Tory = vLG_Tory.fem + vLG_Tory.mal, vLG_Mesopol = vLG_Mesopol + vLG_Ptero.2, vLG_Eulo = vLG_Eulo.fem + vLG_Eulo.mal) %>%
  select(Gender, Genotype, aSG_Tory, rG_Platy, rG_Tory, rG_Mesopol, rG_Eulo, rG_Lestodip, rsLG_Eury, rsLG_Lathro, SG_Platy, vLG_Platy, vLG_Tory, vLG_Mesopol, vLG_Eulo, vLG_Mymarid) 

geno_shootcount = left_join(shoot.count_df, plant.position_df) %>%
  group_by(Genotype) %>%
  summarise(geno_shootEst.all = sum(shootEst.all, na.rm=TRUE), geno_n = n()) # 10 reps of genotype Z? NEED TO MAKE SURE THIS SHOOT COUNT CORRESPONDS EXACTLY WITH THE SAMPLING FOR THE GALLING INSECTS.

geno_gallnet_shootcount = left_join(geno_gallnet, geno_shootcount)
geno_gallnet_density = round((select(geno_gallnet_shootcount, aSG_Tory:vLG_Mymarid)/geno_gallnet_shootcount$geno_shootEst.all)*1000) # number of interactions per 1000 shoots sampled, which was close to the minimum number sampled for each genotype. I rounded the numbers to the nearest whole digit for modularity analysis. 
rownames(geno_gallnet_density) <- geno_gallnet_shootcount$Genotype

library("bipartite")
visweb(geno_gallnet_density)
visweb(geno_gallnet_density, type = "diagonal")

## modularity and nestedness analysis. still need to make sure my null models are appropriate and that I have accounted for shoot densities appropriately.

# create a list with 100 of the same gall networks data frames to calculate and identify the maximum modularity value.
geno_gallnet_samelist <- list()
for(i in 1:100){
  geno_gallnet_samelist[[i]] <- geno_gallnet_density
}

### Still need to try this code after I account for shoot estimates. It will probably take a very long time, which is why I have coded it out. Identify iteration generating the maximum modularity value "Q". Note that any adjustments here should also occur for the null model below!!!
vector.Q_geno_link <- sapply(geno_gallnet_samelist, computeModules)
max.Q_geno_link <- max(sapply(vector.Q_geno_link, function(x) max(x@likelihood)))

# Nestedness
nested_data <- data[order(rowSums(data), decreasing = TRUE), order(colSums(data), decreasing = TRUE)] # creates a nested matrix based on the visweb(type="nested") method
(wt_NODF <- nested(nested_data, method = "weighted NODF"))
# czvalues(geno_link_modul, weighted = TRUE) # defaults to hight trophic level. Note this function remains experimental and still unclear of its relevance to quantitative bipartite networks...
# plotModuleWeb(geno_link_modul)

### Null model analysis. Taken from Dormann and Strauss 2013; Methods in Ecology and Evolution supplement. Note that I did receive a warning indicating that "Some empty columns or rows were deleted"... This take a while to run so don't do it every time!
# r2dtable = Patefield algorithm (preserves marginal totals of network)
# shuffle.web = preserves connectance of observed web
# swap.web = preserves marginal total (r2dtable) and connectance of observed web (shuffle.web)
# note the qualitative connectance, and not weighted connectance is preserved by these null models.
# According to Dormann et al. 2009, using the different null models may give us insight to different processes influencing the structure of the network. 
null.geno_link <- nullmodel(data, N = 100, method = "swap.web")
modules.null <- sapply(null.geno_link, computeModules)
like.nulls <- sapply(modules.null, function(x) x@likelihood)
(z.geno_link <- (max.Q_geno_link - mean(like.nulls))/sd(like.nulls)) 
# Patefield algortihm z-score
# swap.web algortihm z-score: 2.64
# shuffle.web algorithm z-score: 

# Null model for nestedness
wt_NODF_fxn <- function(x) nested(x, method = "weighted NODF")
nested.null <- sapply(null.geno_link, wt_NODF_fxn)
(z.geno_link.wtNODF = (wt_NODF - mean(nested.null))/sd(nested.null))

networklevel(null.geno_link[[2]], index = c("connectance","linkage density","weighted connectance"))



########## EVERYTHING BELOW IS OLD BUT MAYBE USEFUL
### Create data frame of interest
net_focal <- subset(complete_network_data, select = c("Gender","Genotype","plant.position.new","X","gall.sp","gall.id","point.count","g.height","g.width","g.3meas","vLG.pupa","rG.surv","Pont.surv","SG.larv","aSG.larv","Eury.tot","Tory.tot","Meso.tot","Eulo.tot","Lathro.tot","Mymarid","Lestodip","Platy","g.total"))

### Visualize relationship among gall measurement variables. I first explored PCAs with all three gall measurements, however, g.3meas never appeared to contribute unique information. And since I only have g.3meas for a subset of the galls, I have decided to omit it from future analysis and only conduct PCAs on g.height and g.width.

# vLG
vLG.size.data <- na.omit(subset(net_focal, gall.sp == "vLG", select = g.height:g.width))
scatterplotMatrix(vLG.size.data)
vLG.size.pca <- rda(vLG.size.data)
summary(vLG.size.pca)
plot(vLG.size.pca, display="sp")

# rG
rG.size.data <- na.omit(subset(net_focal, gall.sp == "rG", select = g.height:g.width))
scatterplotMatrix(rG.size.data)
rG.size.pca <- rda(rG.size.data)
summary(rG.size.pca)
plot(rG.size.pca, display="sp")

# rsLG
rsLG.size.data <- na.omit(subset(net_focal, gall.sp == "rsLG", select = g.height:g.width))
scatterplotMatrix(rsLG.size.data)
rsLG.size.pca <- rda(rsLG.size.data)
summary(rsLG.size.pca)
plot(rsLG.size.pca, display="sp")

# SG
SG.size.data <- na.omit(subset(net_focal, gall.sp == "SG", select = g.height:g.width))
scatterplotMatrix(SG.size.data)
SG.size.pca <- rda(SG.size.data)
summary(SG.size.pca)
plot(SG.size.pca, display="sp")

# aSG
aSG.size.data <- na.omit(subset(net_focal, gall.sp == "aSG", select = g.height:g.width))
scatterplotMatrix(aSG.size.data) # not really correlated with each other...
aSG.size.pca <- rda(aSG.size.data)
summary(aSG.size.pca)
plot(aSG.size.pca, display="sp")

###### Finalize data set and melt it

# merge gall size PCs back into big dataframe
size.PCs <- rbind.data.frame(scores(vLG.size.pca, display="sites"), scores(rG.size.pca, display="sites"), scores(rsLG.size.pca, display="sites"), scores(SG.size.pca, display="sites"), scores(aSG.size.pca,display="sites"))
size.PCs <- transform(size.PCs, id = rownames(size.PCs))
net_focal$id = factor(rownames(net_focal))

net_focal <- join(net_focal, size.PCs, by="id")
net_focal <- subset(net_focal, select = c(Gender:g.width,vLG.pupa:g.total,PC1:PC2))

# examine whether gall size has an effect on parasitism rates.
## THIS JUST LOOKS AT SURVIVAL AND DOESN'T FOCUS ON PARASITISM AND THEREFORE MAY BE CONFOUNDED A BIT THE "NOTHINGS" INCLUDED IN G.TOTAL. I THINK I NEED TO HAVE SUMMED UP THE PTOID ATTACKS (TRY INCLUDING AND EXCLUDING EXIT HOLES) OVER THE G.TOTAL (AND PTOID ATTACKS + VLG SURVIVAL) TO SEE IF THIS GIVES ANY INSIGHT TO THE CUTOFF VALUE FOR GALL SIZE, OR IF THE MODEL IS A STRONGER PREDICTOR OF PARASITISM.
vLG.bin <- glm(cbind(vLG.pupa, g.total-vLG.pupa) ~ PC1, subset(net_focal, gall.sp == "vLG"), family="binomial")
vLG.lm <- lm(vLG.pupa/g.total ~ PC1, subset(net_focal, gall.sp == "vLG"))
summary(vLG.lm)
plot(vLG.lm)
summary(vLG.bin)
anova(vLG.bin, test = "Chisq")
plot(vLG.bin)
library(MASS)
library(pscl)
library(rms) # look into the potential of this package.
#vLG.lrm <- lrm(cbind(vLG.pupa, g.total-vLG.pupa) ~ g.height, subset(net_focal, gall.sp == "vLG"))
pR2(vLG.bin)
confint(vLG.bin, level=0.95)
dose.p(vLG.bin, p=0.5)
yhat <- fitted(vLG.bin)

plot(vLG.pupa/(g.total) ~ g.height, subset(net_focal, gall.sp == "vLG"))
lines(yhat[order(g.height)] ~ g.height[order(g.height)], subset(net_focal, gall.sp == "vLG"))



gall_net_melt <- melt(net_focal, id.vars = c("Genotype","Gender","plant.position.new","X","gall.sp","gall.id","point.count","g.total","g.height","g.width","PC1","PC2"), variable_name = "gall_contents") # preserve zeros for network analysis, because they identify herbivores that were present, just not connected to any parasitoids.

t <- dcast(gall_net_melt, plant.position.new ~ gall.sp + gall_contents, sum)
plot(rowSums(t[ ,c(60:62,64,66)]) ~ rowSums(t[ ,54:66])) # now I need to account for density...

summary(lm(rowSums(t[ ,c(60:62,64,66)]) ~ rowSums(t[ ,54:66])))


### Regional-scale network (all tree samples summed together). Consider using this regional network to identify the most important species (czvalues function)
total_net <- dcast(gall_net_melt, gall.sp ~ gall_contents, sum)
rownames(total_net) <- total_net[ ,1]
total_net_ptoids <- subset(total_net, select=Eury.tot:Platy) # remove gall.sp
gall_abunds <- rowSums(total_net[,-1])
names(gall_abunds) <- rownames(total_net)

visweb(total_net_ptoids, type="diagonal") # Compartments sort according to phylogeny
visweb(total_net_ptoids)
plotweb(total_net_ptoids, low.abun = gall_abunds - rowSums(total_net_ptoids), abuns.type = "additional")
networklevel(total_net_ptoids)

### Genotype level network with gall-parasitoid links. Need to account for number of shoots sampled on each genotype.

# create sum of estimated shoots sampled for each genotype
shoot.countEst.plantinfo <- join(shoot.countEst, plant.position.info[ ,c("Row..","Plant.Position","Genotype","Gender")])
# note that for plant position 656 I don't have a shoot count estimate for it (DID I REALLY NEVER MEASURE ITS BASAL DIAMETER???)
plot(shootEst.no18 ~ Genotype, shoot.countEst.plantinfo)
summary(lm(shootEst.no18 ~ Genotype, shoot.countEst.plantinfo)) # marginally significant variation in the number of shoots sampled per genotype
shootEst.geno <- aggregate(cbind(shootEst.no18, shootEst.all) ~ Genotype, FUN = sum, shoot.countEst.plantinfo) # calculate the estimated number of shoots sampled per genotype by summing up across all plants.
shootEst.geno.trans <- transform(shootEst.geno, per1000.no18 = shootEst.no18/1000, per1000.all = shootEst.all/1000)

# create genotype by gall-ptoid bipartite network (not controlling for number of shoots sampled)
geno_link_net <- dcast(gall_net_melt, Genotype ~ gall.sp + gall_contents, sum)
geno_link_net_sub <- geno_link_net[ ,-1]
rownames(geno_link_net_sub) <- geno_link_net$Genotype
geno_link_net_sub <- geno_link_net_sub[ ,-which(colSums(geno_link_net_sub) == 0)] # remove gall-parasitoid interactions that never occurred
geno_link_net_sub_noSurv <- geno_link_net_sub[ ,-c(1,3,9,12,14)] # remove columns that indicate survival

networklevel(geno_link_net_sub_noSurv)
visweb(geno_link_net_sub_noSurv, type="diagonal")
visweb(geno_link_net_sub_noSurv, type="nested")
plotweb(geno_link_net_sub_noSurv)

# Control for number of shoots sampled for genotype by gall-ptoid bipartite network. RIGHT NOW I HAVEN'T REMOVED PLANT POSITION 656!!! FOR WHICH I DON'T HAVE A SHOOT ESTIMATE FOR
geno_link_net_shoots <- join(geno_link_net, shootEst.geno.trans)
geno_link_net_shoots_sub <- geno_link_net_shoots[ ,-1]
rownames(geno_link_net_shoots_sub) <- geno_link_net_shoots$Genotype
geno_link_net_shoots_sub <- geno_link_net_shoots_sub[ ,-which(colSums(geno_link_net_shoots_sub) == 0)] # remove gall-parasitoid interactions that never occurred
geno_link_net_shoots_sub_noSurv <- geno_link_net_shoots_sub[ ,-c(1,3,9,12,14)] # remove columns that indicate survival
no18.link.dens.data <- subset(geno_link_net_shoots_sub_noSurv, select = aSG_Tory.tot:vLG_Platy)/geno_link_net_shoots_sub_noSurv$per1000.no18 # frequencies are on a per 1000 shoot basis based on the shoot estimates from omitting no18
all.link.dens.data <- subset(geno_link_net_shoots_sub_noSurv, select = aSG_Tory.tot:vLG_Platy)/geno_link_net_shoots_sub_noSurv$per1000.all # frequencies are on a per 1000 shoot basis based on the shoot estimates from including all inputs for shoot data.

visweb(no18.link.dens.data, type = "diagonal")
visweb(no18.link.dens.data)
plotweb(no18.link.dens.data)

# the functions ceiling, floor, round(x,0) differentially create integers. I will create datasets using all of the different combinations and make sure my results are robust to them. However, if they all show qualitatively the same results, I will use the round(x,0) and the shoot count estimate with all of the shoots for the data.


# create a list with 100 of the same geno.link data frame to calculate and identify the maximum modularity value.
data = round(all.link.dens.data,0)

geno_link_samelist <- list()
for(i in 1:100){
  geno_link_samelist[[i]] <- data
}

### Still need to try this code after I account for shoot estimates. It will probably take a very long time, which is why I have coded it out. Identify iteration generating the maximum modularity value "Q". Note that any adjustments here should also occur for the null model below!!!
vector.Q_geno_link <- sapply(geno_link_samelist, computeModules)
max.Q_geno_link <- max(sapply(vector.Q_geno_link, function(x) max(x@likelihood)))

# Nestedness
nested_data <- data[order(rowSums(data), decreasing = TRUE), order(colSums(data), decreasing = TRUE)] # creates a nested matrix based on the visweb(type="nested") method
(wt_NODF <- nested(nested_data, method = "weighted NODF"))
# czvalues(geno_link_modul, weighted = TRUE) # defaults to hight trophic level. Note this function remains experimental and still unclear of its relevance to quantitative bipartite networks...
# plotModuleWeb(geno_link_modul)

### Null model analysis. Taken from Dormann and Strauss 2013; Methods in Ecology and Evolution supplement. Note that I did receive a warning indicating that "Some empty columns or rows were deleted"... This take a while to run so don't do it every time!
# r2dtable = Patefield algorithm (preserves marginal totals of network)
# shuffle.web = preserves connectance of observed web
# swap.web = preserves marginal total (r2dtable) and connectance of observed web (shuffle.web)
# note the qualitative connectance, and not weighted connectance is preserved by these null models.
# According to Dormann et al. 2009, using the different null models may give us insight to different processes influencing the structure of the network. 
null.geno_link <- nullmodel(data, N = 100, method = "swap.web")
modules.null <- sapply(null.geno_link, computeModules)
like.nulls <- sapply(modules.null, function(x) x@likelihood)
(z.geno_link <- (max.Q_geno_link - mean(like.nulls))/sd(like.nulls)) 
# Patefield algortihm z-score
# swap.web algortihm z-score: 2.64
# shuffle.web algorithm z-score: 

# Null model for nestedness
wt_NODF_fxn <- function(x) nested(x, method = "weighted NODF")
nested.null <- sapply(null.geno_link, wt_NODF_fxn)
(z.geno_link.wtNODF = (wt_NODF - mean(nested.null))/sd(nested.null))

networklevel(null.geno_link[[2]], index = c("connectance","linkage density","weighted connectance"))

############ CODE FOR "Gall-parasitoid network for each genotype (sum all samples within genotype)" BELOW NEEDS TO BE REVIEWED AND MODIFIED

#### Gall-parasitoid network for each genotype (sum all samples within genotype). Right now, this does not account for differences in number of shoots sampled between genotypes (which is likely small)
geno_net <- cast(gall_net_melt, gall.sp ~ gall_contents | Genotype, sum)
#geno_net <- geno_net[c(1:3,5:27)]# remove genotype C as it was never part of the dataset

geno_net_new <- list()
for(i in 1:length(geno_net)){
  tmp <- geno_net[[i]][,-c(1:6)] # remove gall.sp and surviving gall larva data.
  geno_net_new[[i]] <- as.data.frame(empty.col(tmp))
  rownames(geno_net_new[[i]]) <- geno_net[[i]][,1] 
} # empty.col functions removes parasitoids with no interactions but preserves basal species that were observed in the network.

geno_net_metrics <- list()
for(i in 1:length(geno_net_new)){
  geno_net_metrics[[i]] <- network_metrics(as.matrix(geno_net_new[[i]]))
} 
geno_net_data <- ldply(geno_net_metrics)
geno_net_data <- transform(geno_net_data, Geno.ID = c("*","A","B",LETTERS[4:26]), pred.preyRich = predRich/preyRich) 

# Plot relationships among network metrics.
scatterplotMatrix(geno_net_data[,c(1:7,9)])


##### Create networks for each plant position

# create data frames
pp_net <- cast(gall_net_melt, gall.sp ~ gall_contents | plant.position.new, sum) # network for each plant position
pp_gall.ptoid.comm <- dcast(gall_net_melt, plant.position.new ~ gall_contents, sum) # community composition of surviving insects
colnames(pp_gall.ptoid.comm)[1] <- "pp" # for joining with other datasets later

pp_gall.height <- dcast(gall_net_melt, plant.position.new ~ gall.sp, mean, value.var = "g.height")
colnames(pp_gall.height) <- c("pp", paste("height",colnames(pp_gall.height)[-1],sep="."))

pp_gall.width <- dcast(gall_net_melt, plant.position.new ~ gall.sp, mean, value.var = "g.width")
colnames(pp_gall.width) <- c("pp", paste("width",colnames(pp_gall.width)[-1],sep="."))

pp_gall.sizePC1 <- dcast(gall_net_melt, plant.position.new ~ gall.sp, mean, value.var = "PC1")
colnames(pp_gall.sizePC1) <- c("pp", paste("PC1",colnames(pp_gall.sizePC1)[-1],sep="."))

pp_gall.sizePC2 <- dcast(gall_net_melt, plant.position.new ~ gall.sp, mean, value.var = "PC2")
colnames(pp_gall.sizePC2) <- c("pp", paste("PC2",colnames(pp_gall.sizePC2)[-1],sep="."))

# run loop to remove unwanted columns and change rownames of network level data
pp_net_new <- list()
for(i in 1:length(pp_net)){
  tmp <- pp_net[[i]][,-c(1:6)] # # remove gall.sp and surviving gall larva data.
  pp_net_new[[i]] <- as.data.frame(empty.col(tmp))
  rownames(pp_net_new[[i]]) <- pp_net[[i]][,1] 
} # empty.col functions removes parasitoids with no interactions but preserves basal species that were observed in the network.

# calculate network metrics at each plant position
pp_net_metrics <- list()
for(i in 1:length(pp_net_new)){
  pp_net_metrics[[i]] <- network_metrics(as.matrix(pp_net_new[[i]]))
} 
pp_net_data <- ldply(pp_net_metrics) # turn it into a dataframe

# run for loop to extract plant position ID from all networks
pp <- vector()
for(i in 1:length(pp_net)){
  pp[i] <- names(pp_net[i])
}
pp_net_data_trans <- transform(pp_net_data, pp = pp, pred.preyRich = predRich/preyRich) # add plant position ID and predator:prey richness ratios

### merge Gender and Genotype info into pp_net_data as well as shoot count estimates. Add some new variables too.
pp.info <- data.frame(pp = paste("pp",plant.position.info$Plant.Position,sep=""),Geno.ID = plant.position.info$Genotype, Sex = plant.position.info$Gender)
shoot.count.data <- data.frame(pp = paste("pp", shoot.countEst$Plant.Position,sep=""), shootEst.all = shoot.countEst$shootEst.all, shoot.countEst.no18 = shoot.countEst$shootEst.no18, Galls.found = shoot.countEst$Galls.found)

# create list of all dataframes
pp_dfs = list(shoot.count.data, pp_net_data_trans,pp_gall.ptoid.comm, pp_gall.height, pp_gall.width, pp_gall.sizePC1, pp_gall.sizePC2, pp.info)
pp_net_data_merge <- join_all(pp_dfs, by="pp", type = "left")

str(pp_net_data_merge)
# STILL NEED TO DECIDE WHETHER TO OMIT THESE OR NOT. SURELY THEY NEED TO BE APART OF THE ABUNDANCE DATA, BUT WHAT ABOUT GALL SURVIVAL? PROBABLY NOT.
pp_net_data_merge[which(pp_net_data_merge$Galls.found == 0),] # Right now, 77 is a mystery (field notes suggest no galls, but there is a gall collection (possibly from survey #1???)).  Right now, I've decided to keep it in the dataset.

pp_net_data_merge <- transform(pp_net_data_merge, connect_wt_NAzero = connect_wt)
pp_net_data_merge$connect_wt_NAzero[which(is.nan(pp_net_data_merge$connect_wt_NAzero))] <- 0 # since these values had galls but no observed parasitoids, I believe the connectance should be zero...

scatterplotMatrix(pp_net_data_merge[ ,c("vLG.pupa","connect_wt_NAzero")])

# Does gall size vary among genotypes? Need to double check the dcast to see how gall size is measured. I also may need to weight the variables or consider looking at the variance in gall size, or the proportion over a certain threshold value.
plot(height.vLG ~ Geno.ID, pp_net_data_merge)
vLG.height <- lm(height.vLG ~ Geno.ID, pp_net_data_merge)
summary(vLG.height)

# does the number of vLG.pupa vary among genotypes (STILL VERY CRUDE. NEED TO MAKE PROPORTIONAL)
plot(vLG.pupa ~ Geno.ID, pp_net_data_merge)
summary(glm(vLG.pupa ~ Geno.ID, pp_net_data_merge, family="poisson")) # no vLG on Genotype W?!?!?!


### explore data
pp_net_summary <- ddply(pp_net_data_merge, .(Geno.ID), summarize, mean_gallRich = mean(preyRich), mean_ptoidRich = mean(predRich))
max(pp_net_summary$mean_gallRich)
min(pp_net_summary$mean_gallRich)
max(pp_net_summary$mean_ptoidRich)
min(pp_net_summary$mean_ptoidRich)


plot(totRich ~ Geno.ID, pp_net_data_merge)
plot(preyRich ~ Geno.ID, pp_net_data_merge)
plot(predRich ~ Geno.ID, pp_net_data_merge)
plot(connect_wt ~ Geno.ID, pp_net_data_merge) # likely too few data points at the individual tree level to adequately examine this.
plot(connect_wt_NAzero ~ Geno.ID, pp_net_data_merge)
hist(pp_net_data_merge$connect_wt_NAzero)
hist(pp_net_data_merge$connect_wt)

totRich.glm <- glm(totRich ~ Geno.ID, pp_net_data_merge, family="poisson")
summary(totRich.glm)
anova(totRich.glm, test="Chisq")
plot(totRich.glm) # residuals a bit wonky with lm but the conclusions are qualitatively the same as the poisson. Also the R-squared version of glm (1 - residDev/NullDev) is virtually equivalent to the lm

plot(lm(preyRich ~ Geno.ID, pp_net_data_merge))
plot(lm(predRich ~ Geno.ID, pp_net_data_merge))
summary(lm(pred.preyRich ~ Geno.ID, pp_net_data_merge))
summary(lm(connect_wt ~ Geno.ID, pp_net_data_merge))

# weighted connectance model
connect_wt.lm <- lm(connect_wt_NAzero ~ Geno.ID, pp_net_data_merge, weights=g.total) # results weighted by the number of galls collected (not gall density!) on each tree.
summary(connect_wt.lm) # Genotype explains 40.5% of the variance
plot(connect_wt.lm) # scale-location residual plots looks MUCH better when weighted
anova(connect_wt.lm)

# scatterplots exploring relationships amon variables
scatterplotMatrix(pp_net_data_merge[ ,c("preyRich","predRich","totRich","pred.preyRich", "connect_wt_NAzero","g.total")])
scatterplotMatrix(na.omit(pp_net_data_merge[ ,c("V","G","LD_q","connect_wt","preyRich","predRich","totRich","pred.preyRich")]))


### Community composition models. Do I want to use Euclidean distance or a measure that reflects community composition more?
library(vegan)
zeros <- which(rowSums(pp_net_data_merge[ ,c("aSG","rG","rsLG","SG","twG","vLG")]) == 0)
zeros.ptoid <- which(rowSums(pp_net_data_merge[ ,c("Eury.tot","Tory.tot","Meso.tot","Eulo.tot","Lathro.tot","Lestodip","Platy","Mymarid")])==0)

gall.comm <- pp_net_data_merge[-zeros ,c("aSG","rG","rsLG","SG","twG","vLG")]
ptoid.comm <- pp_net_data_merge[-zeros.ptoid ,c("Eury.tot","Tory.tot","Meso.tot","Eulo.tot","Lathro.tot","Lestodip","Platy","Mymarid")]


adonis(gall.comm ~ Geno.ID, method = "horn", data = pp_net_data_merge[-zeros,] )
adonis(ptoid.comm ~ Geno.ID, method = "horn", data = pp_net_data_merge[-zeros.ptoid,])

gall.MDS <- metaMDS(gall.comm, distance = "horn", autotransform = FALSE)
plot(gall.MDS)

rda.gall.comm <- rda(decostand(gall.comm, method = "hell") ~ Geno.ID, data = pp_net_data_merge[-zeros,])
RsquareAdj(rda.gall.comm)
plot(rda.gall.comm, display = c("cn","sp"))
plot(rda.gall.comm) # is this evidence of something weird going on?

# promising that parasitoid community composition still varied among genotypes with a subset of the data
# Need to consider dropping genotypes that have less than 1 rep.
# Consider looking at Links instead of species, and how genotype influences this structure.
rda.ptoid.comm <- rda(decostand(ptoid.comm, method="norm") ~ Geno.ID, data = pp_net_data_merge[-zeros.ptoid,])
RsquareAdj(rda.ptoid.comm)
anova(rda.ptoid.comm, by="axis")
plot(rda.ptoid.comm, display=c("cn","sp"), choices=c(1,2))

### Plots for manuscript

# manage the data
connect_geno_mean <- ddply(pp_net_data_merge, .(Geno.ID), summarize, mean_connect = mean(connect_wt_NAzero))
connect_geno_mean.trans <- transform(connect_geno_mean, Geno.ID.ord = reorder(Geno.ID, mean_connect))


pp_net_data_merge_connect <- transform(pp_net_data_merge, Geno.ID.order = factor(Geno.ID, levels = levels(connect_geno_mean.trans$Geno.ID.ord)))
levels(pp_net_data_merge_connect$Geno.ID.order)[-27] # drop genotype
levels(pp_net_data_merge_connect$Geno.ID.order)[18] <- "C" # replace Genotype * with C, purely for aestheic purposes.

# this needs to come AFTER
levels(connect_geno_mean.trans$Geno.ID.ord)[-27]
levels(connect_geno_mean.trans$Geno.ID.ord)[18] <- "C"


# theme
theme <- theme_bw() + theme(axis.text=element_text(size=16), axis.title = element_text(size = 20, vjust = 0.3), axis.line=element_line(colour="black"), legend.key=element_blank(), panel.grid=element_blank(), panel.border=element_blank(), legend.title = element_blank()) 

# plot
ggplot(pp_net_data_merge_connect, aes(x = Geno.ID.order, y = connect_wt_NAzero)) + geom_jitter(aes(x = Geno.ID.order, y = connect_wt_NAzero), position = position_jitter(w = 0.1, h = 0), shape = 21, size = 4, color = "grey") + ylab("Connectance") + xlab("Genotype") + layer(data = connect_geno_mean.trans, mapping = aes(x = Geno.ID.ord, y = mean_connect), geom = "point", size = 10, shape = 21,fill = "red") +theme

# bubble plot to show which data points had the strongest weights
qplot(data = pp_net_data_merge_connect, x = Geno.ID.order, y = connect_wt_NAzero, size = g.total, position = position_jitter(w = 0.1, h = 0)) + theme + xlab("Genotype") + ylab("Connectance") + ggtitle("Size of points represent their weight in the linear model")


#######  Explore beta diversity of interaction networks
# source in packages from github on March 11, 2014
source('~/Documents/Genotype_Networks/betalink-master/R/metaweb.r')
source('~/Documents/Genotype_Networks/betalink-master/R/betalink.dist.r')
source('~/Documents/Genotype_Networks/betalink-master/R/betalink.r')
source('~/Documents/Genotype_Networks/betalink-master/R/measures.r')


# start off with qualitative metrics
qual.net.dissim <- betalink.dist(pp_net_new, bf=B01)

adonis(qual.net.dissim$WN ~ pp_net_data_merge$Geno.ID)
adonis(qual.net.dissim$OS ~ pp_net_data_merge$Geno.ID)
adonis(qual.net.dissim$S ~ pp_net_data_merge$Geno.ID)
adonis(qual.net.dissim$ST ~ pp_net_data_merge$Geno.ID)
adonis(qual.net.dissim$contrib ~ pp_net_data_merge$Geno.ID)
# so far, I only appeared to find differences in community composition, but none of the other metrics proposed by Poisot et al. 2012. These tree level networks are rather small. It may be better to see which traits are driving dissimilarity in the network at the genotype level, since these webs are likely  more robust. A quantitative metric might be better for the genotype level (maybe at the individual level?)

# What if I create a dataset of all of pairwise interactions at the individual level, and then look at the effect of Genotype? This may ensure that what I'm doing is appropriate.

##### Playing around. SOME OF THIS DATA NEEDS TO BE CLEANED UP.
# transform the data to focus on the parasitoid community only and melt it for easy data summaries
gall_net_all = data.frame(pp = gall_net$plant.position.new, Genotype = gall_net$Genotype, Gender = gall_net$Gender, gall.sp = gall_net$gall.sp, Eury.tot = gall_net$Eury.mal + gall_net$Eury.fem, Tory.tot = gall_net$Tory.mal + gall_net$Tory.fem, Meso.tot = gall_net$Mesopol + gall_net$Ptero.2, Eulo.tot = gall_net$Eulo.fem + gall_net$Eulo.mal, Lathro.tot = gall_net$Lathro.fem + gall_net$Lathro.mal, Lestodip = gall_net$Lestodip, Platy = gall_net$Platy, Mymarid = gall_net$Mymarid, vLG.surv = gall_net$vLG.pupa, rG.surv = gall_net$rG.wh.larv + gall_net$rG.or.larv, SG.surv = gall_net$SG.larv, aSG.surv = gall_net$aSG.larv, twG.surv = gall_net$twG.larv, Pontania.surv = gall_net$Pont.ad + gall_net$Pont.prep) # do I need to include exit.hole in this? May need to reassess rG white/orange larva situation
#  unclear = with(gall_net, moth.larv + unk.ptoid + unk.larv + moth.ad3 + incidental.sp.239 + diff.or.larv + unsure + nothing + inq.dam + exit.hole)
# g.total = rowSums(subset(gall_net, select = vLG.pupa:nothing))

gall_net_all_melt <- melt(gall_net_all, id.vars = c("Genotype","Gender","pp","gall.sp")) # preserve zeros for network analysis, because they identify herbivores that were present, just not connected to any parasitoids.

t <- dcast(gall_net_all_melt, pp ~ gall.sp + variable, sum)
zero_links <- which(colSums(t[,-1]) == 0)
zero_rows <- which(rowSums(t[,-1]) == 0)

t.sub <- t[-zero_rows,-(zero_links+1)]

t.sub.geno <- join(t.sub, pp.info, by = "pp")

adonis(t.sub.geno[ ,2:30] ~ Geno.ID, t.sub.geno, method = "horn")

rda.exp <- rda(decostand(t.sub.geno[ ,2:30], method="hell") ~ Geno.ID, t.sub.geno)
RsquareAdj(rda.exp)
#anova(rda.exp, by = "axis")
plot(rda.exp, display=c("sp","cn"))

cap.exp <- capscale(t.sub.geno[,2:30] ~ Geno.ID, t.sub.geno, distance = "horn")
summary(cap.exp)
plot(cap.exp, display=c("sp","cn"))

cca.exp <- cca(t.sub.geno[,c(2:20,22:30)] ~ Geno.ID, t.sub.geno) # try removing SG2_Platy
summary(cca.exp)
plot(cca.exp, display=c("sp"))
