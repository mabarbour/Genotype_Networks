
#### Simulation of the additive effects of genetic variation on food-web complexity

## source in required data ----
#source('network_management_tree_level.R')
#source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')
tree_level_interaxn_all_plants_traits_size <- read.csv("manuscript/Dryad_data_copy/empirical_data/tree_level_interaxn_all_plants_traits_size.csv")

## load required libraries
#require(bipartite)
library(reshape)
library(reshape2)

library(plyr)
library(dplyr)
library(tidyr)

## change Genotype * to C for plotting aesthetics
levels(tree_level_interaxn_all_plants_traits_size$Genotype)[1] <- "C"

## data frame of interactions
df.interactions <- tree_level_interaxn_all_plants_traits_size %>%
  tbl_df() %>%
  filter(Genotype != "U") %>% # never collected any galls and thus gall-ptoid interactions so I removed it from the dataset. 
  select(Genotype, 
         willow_vLG = vLG_abund,
         willow_rG = rG_abund,
         willow_aSG = aSG_abund,
         willow_SG = SG_abund,
         vLG_Platy, vLG_Mesopol, vLG_Tory, vLG_Eulo, vLG_Mymarid,
         rG_Platy, rG_Mesopol, rG_Tory, rG_Eulo, rG_Lestodip,
         SG_Platy, aSG_Tory)#,
         #vLG.size = vLG.height.mean)#,
         #Platy_abund, Tory_abund, Mesopol_abund, Eulo_abund, Lestodip_abund, Mymarid_abund)

## Function to sample X replicate plants for each genotype, and then sum up the interactions for each genotype. Sampling is done without replacement to make sure unique plant replicates are sampled.
geno.sample <- function(df, reps){
  df.group <- df %>% group_by(Genotype)
  tmp <- sample_n(df.group, size = reps, replace = FALSE)
  tmp.sum <- tmp %>%
    #select(-vLG.size) %>%
    summarise_each(funs(sum)) %>%
    mutate(plants.sampled = reps) 
  tmp.sum
  #vLG.size <- tmp %>% summarise(vLG.size = mean(vLG.size, na.rm = TRUE))
  
  #tmp.join <- left_join(tmp.sum, vLG.size)
}

## Function to sample different levels of genetic variation (i.e. number of genotypes) without replacement. Note also that this is virtually the same as the "geno.sample" function except this time we are sampling different numbers of genotypes.
diversity.sum <- function(df, number.of.genotypes){
  tmp <- sample_n(df, size = number.of.genotypes, replace = FALSE) %>% arrange(Genotype) # arranged so we can later identify and remove duplicate genotype combinations
  genotype.combo <- paste(tmp$Genotype, collapse = " ") # track genotypes sampled for each combination
  tmp.sum <- tmp %>%
    #select(-vLG.size) %>%
    summarise_each(funs(sum)) %>%
    mutate(genotypes.sampled = number.of.genotypes,
           genotype.combo = genotype.combo) %>%
    select(-Genotype)
  
  #tmp.sum$vLG.size.mean <- mean(tmp$vLG.size, na.rm = TRUE)
  #tmp.sum$vLG.size.sd <- sd(tmp$vLG.size, na.rm = TRUE)
  #tmp.sum$vLG.size.IQR <- IQR(tmp$vLG.size, na.rm = TRUE)
  #tmp.sum$vLG.size.mad <- mad(tmp$vLG.size, na.rm = TRUE)
  
  tmp.sum
}

## Sample different levels of genetic variation X times
df.reps <- 50 # replications of data frame for simulation. 40 for non-monoculture sim
sim.reps <- 1000 # replications of simulation
max.genotypes <- 25 # maximum number of genotypes in simulation

# create data frames for simulation
list.sim.food.web <- list()
for(k in 1:df.reps){
  list.sim.food.web[[k]] <- geno.sample(df = df.interactions, 
                                        reps = 4)   
}

df.sim.food.web <- list()
for(k in 1:df.reps){
  sim.food.web <- list()
  for(i in 1:sim.reps){
    food.web <- list()
    for(j in 1:max.genotypes){  
      food.web[[j]] <- diversity.sum(df = list.sim.food.web[[k]], number.of.genotypes = j)
    }
    sim.food.web[[i]] <- ldply(food.web) %>% mutate(sim.number = i)
    #browser()
  }
  tmp.df.sim.food.web <- ldply(sim.food.web) %>% 
    mutate(df.sim.number = k)
  #browser()
  # remove duplicate genotype combinations from the food webs
  dups.pos <- which(duplicated(tmp.df.sim.food.web$genotype.combo) == TRUE)
  df.sim.food.web[[k]] <- tmp.df.sim.food.web[-dups.pos, ]
}

## manage data sets for calculating food-web complexity
sim.food.web.df <- ldply(df.sim.food.web) %>%
  mutate(unique.sim = interaction(df.sim.number, sim.number, genotypes.sampled, sep = "_"))

#write.csv(sim.food.web.df, "~/Documents/Genotype_Networks/data/food web base data 50 reps of 1000 sims 4 reps.csv")

#require(plyr)
#require(dplyr)
#require(tidyr)
#require(reshape2)
#require(reshape)
#sim.food.web.df <- read.csv("~/Documents/Genotype_Networks/data/food web base data 50 reps of 1000 sims 4 reps.csv") %>%
 # select(-X)

sim.info <- sim.food.web.df %>%
  select(unique.sim, df.sim.number, sim.number, genotypes.sampled, genotype.combo, plants.sampled, willow_vLG, willow_rG, willow_aSG, willow_SG)#, 
         #vLG.size.mean, vLG.size.sd, vLG.size.IQR, vLG.size.mad) 

web.df <- sim.food.web.df %>%
  select(-genotype.combo, -genotypes.sampled, -df.sim.number, -sim.number, -plants.sampled) %>% #,
         #-vLG.size.mean, -vLG.size.sd, -vLG.size.IQR, -vLG.size.mad) %>%
  gather(unique.sim, variable)

## for monoculture sims only ----
#sim.food.web.df.mono <- ldply(list.sim.food.web) %>%
 # mutate(sim.number = as.character(rep(1:df.reps, each = max.genotypes)),
  #       unique.sim = interaction(sim.number, plants.sampled, Genotype, sep = "_"))

#sim.info.mono <- sim.food.web.df.mono %>%
 # select(unique.sim, sim.number, plants.sampled, Genotype)

#web.df.mono <- sim.food.web.df.mono %>%
 # select(-plants.sampled, -Genotype, -sim.number) %>%
  #gather(unique.sim, variable)

#web.df <- web.df.mono # only use for monoculture data
#sim.info <- sim.info.mono # only use for monoculture data

## split up variable (trophic levels of pairwise interactions) for creating bipartite webs ----
full.link.split <- colsplit(web.df$variable, "_", names = c("lower","upper"))
full.link.split.df <- cbind(web.df, full.link.split)

## create bipartite food webs for willow-gall and gall-ptoid interactions ----
willow.gall.list <- cast(filter(full.link.split.df, lower == "willow"), lower ~ upper | unique.sim) 
gall.ptoid.list <- cast(filter(full.link.split.df, lower != "willow"), lower ~ upper | unique.sim)  

## Remove NULL values from list, which are due to zero interactions being documented
willow.gall.list.noNULL <- willow.gall.list[!unlist(lapply(willow.gall.list, is.null))] 
gall.ptoid.list.noNULL <- gall.ptoid.list[!unlist(lapply(gall.ptoid.list, is.null))]

## this function calculates weighted linkage density, generality, and vulnerability of bipartite food webs ----
web.measures <- function(web.list){
  #require(bipartite)
  
  web.list.measures <- list()
  unique.sim <- c()
  
  for(i in 1:dim(web.list)){ 
    unique.sim[i] <- names(web.list[i])
    
    tmp.web <- web.list[[i]]
      rownames(tmp.web) <- tmp.web$lower
      web <- as.matrix.data.frame(tmp.web[ ,-1])
      web[is.na(web)] <- 0 # replaces NA with zeros in web, which does not affect the food-web indices I'm interested in.
      
      ## calculate linkage density, generality, and vulnerability. Code was taken from 'networklevel' function in bipartite, which for some reason, was incorrectly calculating linkage density, generality, and vulnerability in small webs.
      preytot.mat <- matrix(rep(colSums(web), NROW(web)), 
                            NROW(web), byrow = TRUE)
      preyprop.mat <- web/preytot.mat
      predtot.mat <- matrix(rep(rowSums(web), NCOL(web)), 
                            NROW(web), byrow = FALSE)
      predprop.mat <- web/predtot.mat
      
      H_Nk <- apply(preyprop.mat, 2, function(x) -sum(x * log(x), na.rm = TRUE))
      H_Pk <- apply(predprop.mat, 1, function(x) -sum(x * log(x), na.rm = TRUE))
      n_Nk <- ifelse(colSums(web) != 0, exp(H_Nk), 0)
      n_Pk <- ifelse(rowSums(web) != 0, exp(H_Pk), 0)
      
      vulnerability <- sum(rowSums(web)/sum(web) * n_Pk)
      generality <- sum(colSums(web)/sum(web) * n_Nk)
      linkage.density <- 0.5 * (vulnerability + generality) # LD_q
    
    web.list.measures[[i]] <- data.frame(linkage.density, vulnerability, generality)
  }
  web.list.measures.df <- ldply(web.list.measures)
  cbind.data.frame(unique.sim, web.list.measures.df)
}

## calculate measures separately for willow-galls and gall-parasitoids.
willow.gall.measures <- web.measures(willow.gall.list.noNULL)
gall.ptoid.measures <- web.measures(gall.ptoid.list.noNULL)

## join the results together and rename the column names
all.measures <- left_join(sim.info, willow.gall.measures, by = "unique.sim") %>%
  left_join(., gall.ptoid.measures, by = "unique.sim") %>%
  rename(., 
         link.density.plant_gall = linkage.density.x,
         vulnerability.plant_gall = vulnerability.x,
         generality.plant_gall = generality.x,
         link.density.gall_ptoid = linkage.density.y,
         vulnerability.gall_ptoid = vulnerability.y,
         generality.gall_ptoid = generality.y)# %>%
all.measures$total_complexity <- rowMeans(all.measures[ ,c("link.density.plant_gall","link.density.gall_ptoid")], na.rm = TRUE) 

## save the results of the simulation as a dataframe.
#write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity simulation 40 reps of 100 sims.csv") # it took ~5 min to run this simulation at 5 reps, so perhapse 40 min at 40 reps?
#write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity simulation 1000 monosims 4 reps.csv")
#write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity simulation 1 reps of 100 sims 4 reps TEST.csv")
write.csv(all.measures, "food web complexity simulation 50 reps of 1000 sims 4 reps.csv") 
#write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity simulation 50 reps of 1000 sims 4 reps.csv") 

