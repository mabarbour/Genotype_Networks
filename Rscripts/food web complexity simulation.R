source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')

require(mvabund)
require(bipartite)

df <- as.data.frame(tree_level_interaxn_all_plants_traits_size) #%>%
  filter(Genotype != "U") # never collected any galls and thus gall-ptoid interactions so I removed it from the dataset. 
genotypes.noU <- levels(df$Genotype)[-21]

full.link.mvabund   <- mvabund(select(df, 
                                      willow_vLG = vLG_abund,
                                      willow_rG = rG_abund,
                                      willow_aSG = aSG_abund,
                                      willow_SG = SG_abund,
                                      vLG_Platy, vLG_Mesopol, vLG_Tory, vLG_Eulo, vLG_Mymarid,
                                      rG_Platy, rG_Mesopol, rG_Tory, rG_Eulo, rG_Lestodip,
                                      SG_Platy, aSG_Tory))
full.link.manyglm <- manyglm(full.link.mvabund ~ Genotype, data = df, family = "negative.binomial")


predict.allcombos <- function(manyglm.object,
                              factor.levels,
                              factor.index, # for extracting column names
                              predict.type = "response"){
  
  # load required libraries
  require(plyr)
  require(dplyr)
  
  predict.list <- list()
  samples.vector <- c()
  for(i in 2:length(factor.levels)){
    samples <- sort(as.character(sample(factor.levels, i)))
    tmp.data <- data.frame(samples)
    colnames(tmp.data) <- names(manyglm.object$data)[factor.index]
    predict.tmp <- round(predict(manyglm.object, newdata = tmp.data, type = predict.type), 1) # round to prevent pseudoaccuracy for each genotype
    predict.list[[i]] <- colMeans(predict.tmp) # take mean across genotypes
    samples.vector[i] <- paste(samples, collapse = " ") # track genotypes used in analysis
  }
  # create data frame from list
  predict.df <- ldply(predict.list) %>%
    mutate(Number.of.Genotypes = 2:length(factor.levels),
           Genotype.combo = samples.vector[-1])
}

N <- 100

full.link.list <- list()
for(i in 1:N){
  full.link.list[[i]] <- predict.allcombos(full.link.manyglm, 
                                           factor.levels = genotypes.noU,
                                           factor.index = 2)
} # takes a couple minutes with 50 reps
full.link.df <- ldply(full.link.list) %>%
  melt(id.vars = c("Number.of.Genotypes","Genotype.combo"))

## get predicted links for each genotype alone
single.pred.geno <- as.data.frame(
  round(  
    predict.manyglm(full.link.manyglm, 
                    newdata = data.frame(Genotype = genotypes.noU),
                    type = "response"), 
    1) # round to 1-digit to avoid pseudoaccuracy
  ) %>% 
  mutate(Genotype.combo = genotypes.noU,
         Number.of.Genotypes = 1) %>%
  melt(id.vars = c("Number.of.Genotypes","Genotype.combo"))

## combine predicted food webs with multiple genotypes with single genotype food webs
full.link.df <- rbind.data.frame(full.link.df, single.pred.geno)

## split up variable (trophic levels of pairwise interactions) for creating bipartite webs
full.link.split <- colsplit(full.link.df$variable, "_", names = c("lower","upper"))
full.link.split.df <- cbind(full.link.df, full.link.split)

## create bipartite food webs for willow-gall and gall-ptoid interactions. This enables me to accurately calculated linkage density. 
willow.gall.list <- cast(filter(full.link.split.df, lower == "willow"), lower ~ upper | Genotype.combo + Number.of.Genotypes, mean) # by taking the mean, I control for the fact that some of the replicate combinations may be duplicates.
gall.ptoid.list <- cast(filter(full.link.split.df, lower != "willow"), lower ~ upper | Genotype.combo + Number.of.Genotypes, mean) 

## this function calculates weighted linkage density, generality, and vulnerability of bipartite food webs
web.measures <- function(web.list){
  require(bipartite)
  
  web.list.measures <- list()
  genotype.number <- c()
  for(i in 1:dim(web.list)){
    genotype.number[i] <- names(web.list[[i]][1])
    tmp.web <- web.list[[i]][[1]]
    rownames(tmp.web) <- tmp.web$lower
    tmp.web <- as.matrix.data.frame(tmp.web[ ,-1])
    web.list.measures[[i]] <- networklevel(tmp.web, index = c("linkage density","generality","vulnerability"), weighted = TRUE)
  }
  web.list.measures.df <- ldply(web.list.measures)
  genotype.combos <- names(web.list)
  cbind.data.frame(web.list.measures.df, 
                   Number.of.Genotypes = as.numeric(genotype.number), 
                   Genotype.combos = genotype.combos)
}

## calculate measures separately for willow-galls and gall-parasitoids.
willow.gall.measures <- web.measures(willow.gall.list)
gall.ptoid.measures <- web.measures(gall.ptoid.list)

## join the results together and rename the column names
all.measures <- left_join(willow.gall.measures, gall.ptoid.measures, by = c("Genotype.combos","Number.of.Genotypes"))
colnames(all.measures) <- c("link.density.plant_gall", "generality.HL.plant_gall", "vulnerability.LL.plant_gall", "Number.of.Genotypes", "Genotype.combos", "link.density.gall_ptoid", "generality.HL.gall_ptoid", "vulnerability.LL.gall_ptoid")

## calculating total complexity in multiple ways gives the same results.
all.measures <- mutate(all.measures,
                       total_complexity = (link.density.plant_gall + link.density.gall_ptoid)/2,
                       total_complexity.check = (generality.HL.plant_gall + vulnerability.LL.plant_gall + generality.HL.gall_ptoid + vulnerability.LL.gall_ptoid)/4)

## save the results of the simulation as a dataframe.
write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity simulation.csv")

