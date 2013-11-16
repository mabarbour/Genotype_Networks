
##### Data management
# set working directory and upload data
setwd(dir='~/Documents/Genotype_Networks/')
source('~/Documents/miscellaneous_R/community_explore.R')
library(reshape)
library(ggplot2)

Gall_Network_Data_Survey_2_2012 <- read.csv("~/Documents/Genotype_Networks/Gall_Network_Data_Survey_2_raw.csv", skip=1, stringsAsFactors = FALSE, strip.white = TRUE)

str(Gall_Network_Data_Survey_2_2012)
dimnames(Gall_Network_Data_Survey_2_2012)

plant.position.info <- read.csv("~/Documents/Genotype_Willow_Community/datasets/Willow Garden Positions.csv")
plant.position.info <- plant.position.info[ ,c("Plant.Position","Genotype","Gender")]

g.contents_matrix <- as.matrix(subset(Gall_Network_Data_Survey_2_2012, select = vLG.pupa:exit.hole))
str(g.contents_matrix)
g.contents_matrix[is.na(g.contents_matrix)] <- 0 # replaces all NAs in gall contentx matrix with zeros.  These are biologically meaningful zeros and do not represent missing data. May need to combine inq.dam category with exit hole

gall_net <- cbind.data.frame(Gall_Network_Data_Survey_2_2012[ ,c("survey","plant.position","gall.id","g.id.unk","gall.sp","g.sp.unsure","g.height","g.width","g.3meas","point.count", "exit.size")], g.contents_matrix)

gall_net <- merge(gall_net, plant.position.info, by.x="plant.position", by.y="Plant.Position")

gall_net$g.total <- rowSums(subset(gall_net, select = vLG.pupa:nothing)) # counts up the total number of organisms observed in the galls.  Definitely one outlier with an rG with 14 organisms...


### explore data

# Iteomyia community rank abundance plots
rankAbundancePlot(subset(gall_net, gall.sp == "vLG", select = vLG.pupa:exit.hole))

# Rabdophaga salicisbrassicoides rank abundance plots
rankAbundancePlot(subset(gall_net, gall.sp == "rG", select = vLG.pupa:exit.hole))

# Pontania rank abundance plots
rankAbundancePlot(subset(gall_net, gall.sp == "rsLG", select = vLG.pupa:exit.hole))

# effects of gall height on vLG community

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


#### Create Networks

library(bipartite)

gall_net_melt <- melt(gall_net, id.vars = c("Genotype", "Gender", "plant.position","gall.sp"), measure.vars=c("vLG.pupa","rG.wh.larv","rG.or.larv","SG.larv","aSG.larv","twG.larv","Pont.ad","Pont.prep","Pont.larv","Eury.fem","Eury.mal","Platy","Tory.mal","Tory.fem","Mesopol","Ptero.2","Eulo.fem","moth.ad1","Lestodip","Lathro.fem","Lathro.mal","moth.larv","unk.ptoid","Eulo.mal","unk.larv","moth.ad3","incidental.sp.239","Mymarid","diff.or.larv","rsLG.solid","unsure","nothing","inq.dam","exit.hole","g.total"))

# analysis below isn't accurate at all.  Needs to be based on average interaction strengths and gall densities.

genotype_gall_net <- cast(data = gall_net_melt, formula = Genotype ~ gall.sp, sum, subset = variable == "g.total")
genotype_gall_matrix <- as.matrix(genotype_gall_net[ ,-1])
rownames(genotype_gall_matrix) <- genotype_gall_net$Genotype
colnames(genotype_gall_matrix) <- colnames(genotype_gall_net[ ,-1])

visweb(genotype_gall_matrix, type ="diagonal")

plotweb(genotype_gall_matrix)

test <- computeModules(genotype_gall_matrix)
plotModuleWeb(test)


# Gall-ptoid network. Note that these data are still not realistic, because they don't take into account densities or interaction strength.

gall_ptoid_net <- cast(gall_net_melt, gall.sp ~ variable, sum, subset = variable %in% c("Eury.fem","Eury.mal","Platy","Tory.mal","Tory.fem","Mesopol","Ptero.2","Eulo.fem","moth.ad1","Lestodip","Lathro.fem","Lathro.mal","moth.larv","Eulo.mal", "moth.ad3","Mymarid","diff.or.larv"))
gall_ptoid_net$gall.sp[1] <- "UNKNOWN"

gall_ptoid_matrix <- as.matrix(gall_ptoid_net[ ,-1])
colnames(gall_ptoid_matrix) <- colnames(gall_ptoid_net[ ,-1])
rownames(gall_ptoid_matrix) <- gall_ptoid_net$gall.sp

visweb(gall_ptoid_matrix, type="diagonal")
visweb(gall_ptoid_matrix, type="nested")
plotweb(gall_ptoid_matrix)

test_gall_ptoids <- computeModules(gall_ptoid_matrix)
plotModuleWeb(test_gall_ptoids)
