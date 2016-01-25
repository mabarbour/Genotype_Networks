#test <- left_join(gall_net, plant.position_df)
#gall_net <- merge(gall_net, plant.position.info, by.x="plant.position", by.y="Plant.Position")
#gall_net <- transform(gall_net, plant.position.new = paste("pp",plant.position,sep=""))

### Data that I need
# plant position, genotype, count of observed gall-parasitoid species interactions
# plant position, genotype, mean and variance of gall size for each species
# plant position, genotype, specific gall-parsitoid interaction link (i.e. have to know for sure that this gall sp. came from this exact gall to match up size with interaction [mechanistic link]. Need to exclude all gall.id.unk for this calculation. Note that this should not matter for the mean/variance calculations because these galls didn't have any size measurements associated with them above since they don't require gall specific links.)

# is it okay to include galls for which had the category "nothing" into the mean/variance calculations???

#### Source in functions
source('~/Documents/Genotype_Networks/Rscripts/fxns_networks.R')

#### Load libraries
#library(plyr)
library(dplyr)

table(gall_net$point.count)
filter(gall_net, point.count == 1)

gall_net %>%
  select(plant.position, gall.sp, gall.id, g.height,point.count) %>%
  filter(gall.sp == "vLG")
select(gall_net, plant.position, gall.sp, gall.id, point.count)

### DID I REMOVE MOTH.AD1? (SEE LINE 68), incidental.sp.239? (see line 69)
# ADDRESS INQ.DAM, moth.ad3
# ADDRESS UNK.LARV CONCERN (LOOK AT LAB SPECIMENS)
# ADDRESS UNSURE, moth.larv,

hist(gall_net[which(gall_net$nothing == 1 & gall_net$gall.sp == "vLG"), "g.height" ])
hist(gall_net[which(gall_net$gall.sp == "vLG"), "g.height" ])
table(gall_net[ ,c("gall.sp","nothing")]) # not a huge problem for vLG and rsLG, but a problem for others...

# transform the data to focus on the parasitoid community only and melt it for easy data summaries
gall_net_trans = mutate(gall_net, Eury.tot = Eury.mal + Eury.fem, Tory.tot = Tory.mal + Tory.fem, Meso.tot = Mesopol + Ptero.2, Eulo.tot = Eulo.fem + Eulo.mal, Lathro.tot = Lathro.fem + Lathro.mal, Pont.surv = Pont.ad + Pont.prep, rG.surv = rG.wh.larv + rG.or.larv) 

test = tbl_df(gall_net_trans)
filter(test, gall.sp == "vLG", Platy == 1)

#ddply(gall_net_trans, .(plant.position.new, gall.sp), summarise, gall_count = count(gall.id), mean_size = mean(g.height, na.rm = TRUE), sd_size = sd(g.height, na.rm = TRUE), cv_size = sd_size/mean_size)

t <- summarySE(gall_net_trans, measurevar = "g.height", groupvars = c("plant.position.new","gall.sp","Gender","Genotype"), na.rm = TRUE, conf.interval = 0.95, .drop = TRUE)

library(ggplot2)
ggplot(subset(t, gall.sp == "vLG"), aes(x = g.height, y = sd, color = Genotype)) + geom_point()
ggplot(subset(t, gall.sp == "vLG"), aes(x = N, y = g.height, color = Genotype)) + geom_point()
ggplot(subset(t, gall.sp == "vLG"), aes(x = N, y = sd, color = Genotype)) + geom_point()
ggplot(subset(t, gall.sp == "rG"), aes(y = g.height, x = Genotype)) + geom_boxplot()
ggplot(subset(t, gall.sp == "vLG"), aes(y = sd, x = Genotype)) + geom_boxplot()
ggplot(subset(t, gall.sp == "vLG"), aes(y = sd/g.height, x = Genotype)) + geom_boxplot()

summary(lm(log(g.height) ~ Genotype, data = subset(t, gall.sp == "rG")))
summary(lm(log(sd) ~ Genotype, data = subset(t, gall.sp == "vLG")))
summary(lm(log(sd/g.height) ~ Genotype, data = subset(t, gall.sp == "vLG")))





# add another variable: counts up the total number of organisms observed in the galls.  Definitely one outlier with an rG with 14 organisms...
#gall_net$g.total <- rowSums(subset(gall_net, select = vLG.pupa:nothing))
#hist(gall_net$g.total, right = FALSE) # shows number of zeros, ones, twos, etc. Definitely worth noting that many galls had zero contents. Will need to account for whether there was an exit hole or not.

################## MAKE THE DATA TIDY. COLLAPSE GALL CONTENTS INTO ONE COLUMN...
### explore data
table(gall_net$gall.sp)

# Iteomyia community rank abundance plots
rankAbundancePlot(subset(gall_net, gall.sp == "vLG", select = vLG.pupa:exit.hole))

# Rabdophaga salicisbrassicoides rank abundance plots
rankAbundancePlot(subset(gall_net, gall.sp == "rG", select = vLG.pupa:exit.hole)) # note that there may be some "misidentifications", e.g. Pontania

# Pontania rank abundance plots
rankAbundancePlot(subset(gall_net, gall.sp == "rsLG", select = vLG.pupa:exit.hole))

# Relationship between gall height and gall community

vLG_comm <- subset(gall_net, gall.sp == "vLG" & nothing == 0, select = c("g.height","vLG.pupa","Platy","Tory.mal","Tory.fem","Mesopol","Ptero.2","Eulo.fem","Eulo.mal","Genotype")) # currently excluding a few species

vLG_comm$total <- rowSums(vLG_comm[ ,c("vLG.pupa","Platy","Tory.mal","Tory.fem","Mesopol","Ptero.2","Eulo.fem","Eulo.mal")]) # currently excluding 'nothing' categories
vLG_comm$Tory.all <- rowSums(vLG_comm[ ,c("Tory.mal","Tory.fem")])
vLG_comm$Mesopol.all <- rowSums(vLG_comm[ ,c("Mesopol","Ptero.2")])
vLG_comm$Euloph.all <- rowSums(vLG_comm[ ,c("Eulo.fem","Eulo.mal")])

vLG_comm <- with(vLG_comm, cbind(g.height, vLG_comm[ ,c("vLG.pupa","Platy","Tory.all","Mesopol.all","Euloph.all")]/total, total, Genotype)) # note that there are still interesting dynamics if you keeps males and females separate.
vLG_comm <- na.omit(vLG_comm)

vLG_comm_melt <- melt(vLG_comm, id.vars=c("g.height","total","Genotype"))

ggplot(data=vLG_comm_melt, aes(x=g.height, y=value, weight=total, color=variable)) + geom_point() + geom_smooth(se=TRUE) 

# explore distributions in gall height
qplot(data=subset(vLG_comm_melt, Genotype == c("J","X","*","Z","B","H","W")), x=g.height, color = Genotype, geom="density") # focuses on a subset of genotypes

ggplot(data=vLG_comm_melt, aes(x=Genotype, y=g.height)) + geom_boxplot()

# effects of gall height on Pontania community

rsLG_comm <- subset(gall_net, gall.sp == "rsLG" & nothing == 0 & Pont.larv == 0, select = c("g.height","Pont.ad","Pont.prep","Eury.fem","Eury.mal")) # currently excluding 'nothing' cateogires.

rsLG_comm$total <- rowSums(rsLG_comm[ ,c("Pont.ad","Pont.prep","Eury.fem","Eury.mal")])
rsLG_comm$Pont.all <- rowSums(rsLG_comm[ ,c("Pont.ad","Pont.prep")])

rsLG_comm <- with(rsLG_comm, cbind(g.height, rsLG_comm[ ,c("Pont.all","Eury.fem","Eury.mal")]/total, total))
rsLG_comm <- na.omit(rsLG_comm)

rsLG_comm_melt <- melt(rsLG_comm, id.vars=c("g.height","total"))

ggplot(data=rsLG_comm_melt, aes(x=g.height, y=value, weight=total, color=variable)) + geom_point() + geom_smooth(se=FALSE) 


#### Create Networks. Probably need to modify.
library(bipartite)

gall_net_trans = data.frame(pp = gall_net$plant.position.new, Genotype = gall_net$Genotype, Gender = gall_net$Gender, gall.sp = gall_net$gall.sp, Eury.tot = gall_net$Eury.mal + gall_net$Eury.fem, Tory.tot = gall_net$Tory.mal + gall_net$Tory.fem, Meso.tot = gall_net$Mesopol + gall_net$Ptero.2, Eulo.tot = gall_net$Eulo.fem + gall_net$Eulo.mal, Lathro.tot = gall_net$Lathro.fem + gall_net$Lathro.mal, Lestodip = gall_net$Lestodip, Platy = gall_net$Platy, Mymarid = gall_net$Mymarid)

gall_net_melt <- melt(gall_net_trans, id.vars = c("Genotype","Gender","pp","gall.sp")) # preserve zeros for network analysis, because they identify herbivores that were present, just not connected to any parasitoids.

#gall_net_melt <- melt(gall_net, id.vars = c("Genotype", "Gender", "plant.position.new","gall.sp"), measure.vars=c("vLG.pupa","rG.wh.larv","rG.or.larv","SG.larv","aSG.larv","twG.larv","Pont.ad","Pont.prep","Pont.larv","Eury.fem","Eury.mal","Platy","Tory.mal","Tory.fem","Mesopol","Ptero.2","Eulo.fem","moth.ad1","Lestodip","Lathro.fem","Lathro.mal","moth.larv","unk.ptoid","Eulo.mal","unk.larv","moth.ad3","incidental.sp.239","Mymarid","diff.or.larv","rsLG.solid","unsure","nothing","inq.dam","exit.hole"), na.rm=TRUE) # don't need to include "g.total", because I can use "margins" to sum it up. Using na.rm = TRUE when I did not substitute zeros for NAs. Replaced plant.position with plant.position.new which is a character for making later indexing easeir.

# analysis below isn't accurate at all.  Needs to be based on average interaction strengths and gall densities.

#genotype_gall_net <- cast(data = gall_net_melt, formula = Genotype ~ gall.sp, sum, subset = variable == "g.total")
#genotype_gall_matrix <- as.matrix(genotype_gall_net[ ,-1])
#rownames(genotype_gall_matrix) <- genotype_gall_net$Genotype
#colnames(genotype_gall_matrix) <- colnames(genotype_gall_net[ ,-1])

#visweb(genotype_gall_matrix, type ="diagonal")

#plotweb(genotype_gall_matrix)

#test <- computeModules(genotype_gall_matrix)
#plotModuleWeb(test)


# Gall-ptoid network. Note that these data are still not realistic, because they don't take into account densities or interaction strength.

#gall_ptoid_net <- cast(gall_net_melt, gall.sp ~ variable, sum, subset = variable %in% c("Eury.fem","Eury.mal","Platy","Tory.mal","Tory.fem","Mesopol","Ptero.2","Eulo.fem","moth.ad1","Lestodip","Lathro.fem","Lathro.mal","moth.larv","Eulo.mal", "moth.ad3","Mymarid","diff.or.larv"))
#gall_ptoid_net$gall.sp[1] <- "UNKNOWN"

#gall_ptoid_matrix <- as.matrix(gall_ptoid_net[ ,-1])
#colnames(gall_ptoid_matrix) <- colnames(gall_ptoid_net[ ,-1])
#rownames(gall_ptoid_matrix) <- gall_ptoid_net$gall.sp

#visweb(gall_ptoid_matrix, type="diagonal")
#visweb(gall_ptoid_matrix, type="nested")
#plotweb(gall_ptoid_matrix)

#test_gall_ptoids <- computeModules(gall_ptoid_matrix)
#plotModuleWeb(test_gall_ptoids)

### Create networks for each plant position
test <- cast(gall_net_melt, gall.sp ~ variable | pp, sum)#, margins=TRUE) # likely need to incorporate the density per interaction on a per shoot basis. PERHAPS I SHOULD ACTUALLY RETAIN THE ZEROS. AS THEY MAY BE VERY IMPORTANT IN CALCULATING CONNECTANCE...TRY CALCULATING MATRICES WITH AND WITHOUT ZEROS

# something weird going on with pp90

# NEED TO REMOVE COLSUMS = 0, BECAUSE THOSE PARASITOIDS WERE NEVER THERE, DON'T REMOVE ROWSUMS = 0, BECAUSE THOSE SPECIES WERE ACTUALLY OBSERVED.



# convert first column variable into rowname and remove it from the dataframe
#network_pp <- list()
#for(i in 1:length(test)){
# network_pp[[i]] <- networklevel(empty.col(test[[i]][,-1]), empty.web=FALSE) # I have already emptied the web of extra columns, because those parasitoids were not observed, but I'm preserving rows with no interactions because those gall species were observed.
#rownames(test[[i]]) <- test[[i]][,1]
#test[[i]] <- test[[i]][,-1]
}

#test1 <- list()
#for(i in 1:length(test)){
if(sum(test[[i]]) == 0){
  test1[[i]] <- 0
}
if(sum(test[[i]]) > 0){
  #rownames(test[[i]]) <- test[[i]][,1]
  test1[[i]] <- empty.col(test[[i]])
}
}

#test2 <- list()
#for(i in 1:length(test)){
if(any(dim(test[[i]]) > 1)){
  #rownames(test[[i]]) <- test[[i]][,1]
  test2[[i]] <- empty.col(test[[i]])
} else test2[[i]] <- 0
}

test3 <- list()
for(i in 1:length(test)){
  test3[[i]] <- empty.col(test[[i]])
}

### having trouble calculating this at the individual tree level...AFTER WORKING WITH THIS, IT SHOULD WORK THOUGH. I TRIED TESTING A VARIETY OF DIFFERENT MATRICES.

# calculate network level functions. Not completely working right now. Are there some plant positions with no data? This however would be the general correct method for indexes of interest from each value of the list
network_data <- vector()
for(i in 1:length(test3)){
  if(ncol(test3[[i]]) == 0){
    network_data[i] <- 0 # if there is a value of zero, then connectance is equal to zero
  }  else if(ncol(test3[[i]]) > 1 | nrow(test1[[i]]) > 1){
    tmp <- networklevel(test3[[i]],index = c("linkage density","weighted connectance"), empty.web=FALSE)
    network_data[i] <- tmp[2]
  }
}
#else if(dim(test1[[i]])[1] == 1 & dim(test1[[i]])[2] == 1){
#network_data[i] <- 1# need to modify this once I start taking into account per shoot density...
#}

t1 <- matrix(c(0,0,0,0), nrow=2, byrow=T)
t2 <- matrix(c(1,1,0,1,1,0), nrow=2, byrow=T)
networklevel(t1, index = c("linkage density","weighted connectance"), empty.web=FALSE)
networklevel(t2, index = c("linkage density","weighted connectance"), empty.web=FALSE)

t3 <- as.matrix(empty.col(test[[1]]))

# may need to use weighted least squares regression to account for small numbers of certain data points for some networks. Or would it be fine to just use the number of sampled shoots as a covariate? Note that I could divide this by the number of shoots sampled if necessary, because "networklevel" will take fractional data.

# Need to do some data clean up and some decisions made about what to include in the food web matrix
plotweb(test[[1]])
test[[1]]
