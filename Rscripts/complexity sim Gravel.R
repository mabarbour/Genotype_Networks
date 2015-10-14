# In order to get an expectation of how my data is behaving, I could randomly shuffle everything and see what patterns I observe. Alternatively, I could create different scenarios to see what is going on (ref notebook). Let's try shuffling things up first...

# PERHAPS BOOTSTRAPPING IS THE WAY TO GO. NEED TO THINK ABOUT THE ASSUMPTIONS BEHIND THIS A BIT MORE.
# bootstrapping relies on one crucial assumption: the data accurately represents the true population...
# add a calculation of functional diversity and see if this has a stronger relationship with food-web complexity.

# Gravel simulation


source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')

#require(mvabund)
require(bipartite)
require(tidyr)
#require(FD)

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
         #Platy_abund, Tory_abund, Mesopol_abund, Eulo_abund, Lestodip_abund, Mymarid_abund)

## Function to sample X replicate plants for each genotype, and then sum up the interactions for each genotype. Sampling is done without replacement to make sure unique plant replicates are sampled.
geno.sample <- function(df, reps){
  df.group <- df %>% group_by(Genotype)
  tmp <- sample_n(df.group, size = reps, replace = FALSE)
  tmp.sum <- tmp %>%
    summarise_each(funs(sum)) %>%
    mutate(plants.sampled = reps) 
}

## Function to sample different levels of genetic variation (i.e. number of genotypes) without replacement. Essentially, this peforms a sample-based rarefaction of the data, where samples correspond to unique genotypes. Note also that this is virtually the same as the "geno.sample" function except this time we are sampling different numbers of genotypes.
diversity.sum <- function(df, number.of.genotypes){
  tmp <- sample_n(df, size = number.of.genotypes, replace = FALSE) %>% arrange(Genotype) # arranged so we can later identify and remove duplicate genotype combinations
  genotype.combo <- paste(tmp$Genotype, collapse = " ") # track genotypes sampled for each combination
  tmp.sum <- tmp %>%
    summarise_each(funs(sum)) %>%
    mutate(genotypes.sampled = number.of.genotypes,
           genotype.combo = genotype.combo) %>%
    select(-Genotype)
}

## Sample different levels of genetic variation X times
df.reps <- 40 # replications of data frame for simulation
sim.reps <- 100 # replications of simulation
max.genotypes <- 25 # maximum number of genotypes in simulation
#df.sim <- geno.sample(df.interactions, reps = 4) # data frame for simulation

list.sim.food.web <- list()
for(k in 1:df.reps){
  list.sim.food.web[[k]] <- geno.sample(df = df.interactions, reps = 4) # create data frame for simulation  
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
 


#sim.food.web <- list()
#for(i in 1:sim.reps){
 # food.web <- list()
  #for(j in 1:max.genotypes){
   # food.web[[j]] <- diversity.sum(df = df.sim, number.of.genotypes = j)
  #}
  #sim.food.web[[i]] <- ldply(food.web) %>%
   # mutate(sim.number = i)
#}

sim.food.web.df <- ldply(df.sim.food.web) %>%
  mutate(unique.sim = interaction(df.sim.number, sim.number, genotypes.sampled, sep = "_"))

sim.info <- sim.food.web.df %>% 
  select(unique.sim, df.sim.number, sim.number, genotypes.sampled, genotype.combo, plants.sampled)

web.df <- sim.food.web.df %>%
  select(-genotype.combo, -genotypes.sampled, -df.sim.number, -sim.number, -plants.sampled) %>%
  gather(unique.sim, variable)

## split up variable (trophic levels of pairwise interactions) for creating bipartite webs
full.link.split <- colsplit(web.df$variable, "_", names = c("lower","upper"))
full.link.split.df <- cbind(web.df, full.link.split)

## create bipartite food webs for willow-gall and gall-ptoid interactions. This enables me to accurately calculated linkage density. 
willow.gall.list <- cast(filter(full.link.split.df, lower == "willow"), lower ~ upper | unique.sim) 
gall.ptoid.list <- cast(filter(full.link.split.df, lower != "willow"), lower ~ upper | unique.sim)  

## this function calculates weighted linkage density, generality, and vulnerability of bipartite food webs
web.measures.Gravel <- function(web.list){
  require(bipartite)
  
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
willow.gall.list.noNULL <- willow.gall.list[!unlist(lapply(willow.gall.list, is.null))] # Remove NULL values from list, which are due to zero interactions being documented
gall.ptoid.list.noNULL <- gall.ptoid.list[!unlist(lapply(gall.ptoid.list, is.null))]

willow.gall.measures <- web.measures.Gravel(willow.gall.list.noNULL)
gall.ptoid.measures <- web.measures.Gravel(gall.ptoid.list.noNULL)


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
write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity simulation 40 reps of 100 sims.csv") # it took ~5 min to run this simulation at 5 reps, so perhapse 40 min at 40 reps?


#### EVERYTHING BELOW IS OLD, BUT MAYBE USEFUL LATER
## Rarefaction function ----
rarefact <- function(df){
  rows <- nrow(df)
  
  samp.summary.list <- list()
  for(i in 1:rows){
    tmp <- sample_n(df, size = i, replace = FALSE)
    
    samp.summary.list[[i]] <- tmp %>%
      summarise_each(funs(mean, sum)) %>% 
      select(-Genotype_mean, -Genotype_sum) %>%
      mutate(reps.sampled = i)
  }
  samp.summary.df <- ldply(samp.summary.list)
}

## Resampling function ----
get.info <- function(df, sample.size, max.genotypes){
  sample.info <- list()
  for(i in 1:max.genotypes){
    
    ## sample i genotypes
    tmp <- filter(df, 
                  Genotype %in% sample(df$Genotype, i)) %>%
      group_by(Genotype)
    
    ## sample "sample.size" willows of each genotype
    tmp.samp <- sample_n(tmp, 
                         size = sample.size, 
                         replace = FALSE)
    
    ## calculate the average frequency of each willow-gall and gall-parasitoid interaction
    samp.summary <- tmp.samp %>%
      ungroup() %>%
      summarise_each(funs(sum)) %>% # funs(mean,sum)
      select(-Genotype) #-Genotype_mean, -Genotype_sum
    
    ## calculate total gall abundance and gall-parasitoid interactions from the averaged interactions above.
    #total.gall.abund <- rowSums(select(samp.summary, willow_vLG_sum, willow_rG_sum, willow_aSG_sum, willow_SG_sum))#rowSums(select(samp.summary, starts_with("willow"))) 
    
    #total.gall.rich <- rowSums(select(samp.summary, willow_vLG_sum, willow_rG_sum, willow_aSG_sum, willow_SG_sum) > 0)
    
    #total.gall_ptoid.abund <- rowSums(select(samp.summary, vLG_Platy_sum, vLG_Mesopol_sum, vLG_Tory_sum, vLG_Eulo_sum, vLG_Mymarid_sum, rG_Platy_sum, rG_Mesopol_sum, rG_Tory_sum, rG_Eulo_sum, rG_Lestodip_sum, SG_Platy_sum, aSG_Tory_sum)) #rowSums(select(samp.summary, vLG_Platy:aSG_Tory))
    
    #total.gall_ptoid.rich <- rowSums(select(samp.summary, vLG_Platy_sum, vLG_Mesopol_sum, vLG_Tory_sum, vLG_Eulo_sum, vLG_Mymarid_sum, rG_Platy_sum, rG_Mesopol_sum, rG_Tory_sum, rG_Eulo_sum, rG_Lestodip_sum, SG_Platy_sum, aSG_Tory_sum) > 0)
    
    #total.ptoid.rich <- rowSums(select(samp.summary, Platy_abund_sum, Tory_abund_sum, Mesopol_abund_sum, Eulo_abund_sum, Lestodip_abund_sum, Mymarid_abund_sum) > 0)

    ## keep track of number of genotypes and number of plants sampled for each genotype
    genotype.richness <- i
    reps.sampled <- sample.size

    ## collate information for each sample
    sample.info[[i]] <- samp.summary %>%
      mutate(genotype.richness = genotype.richness,
             reps.sampled = reps.sampled)#,
             #total.gall.abund = total.gall.abund,
             #total.gall_ptoid.abund = total.gall_ptoid.abund,
             #total.gall.rich = total.gall.rich,
             #total.ptoid.rich = total.ptoid.rich,
             #total.gall_ptoid.rich = total.gall_ptoid.rich,
             #total.rich = total.gall.rich + total.ptoid.rich,
             #total.interact.abund = total.gall.abund + total.gall_ptoid.abund,
             #total.interact.rich = total.gall.rich + total.gall_ptoid.rich)
  }
  ## condense list into a data frame
  sample.info.df <- ldply(sample.info)   
}

## Run plant-based rarefaction simulation ----
df.rarefact.sim <- list()
for(j in 1:100){
  get.df <- rarefact(df = df.interactions)
  df.rarefact.sim[[j]] <- cbind.data.frame(sim.number = j, get.df)
}


## Run simulation ----
df.sample.list.sim <- list()
for(j in 1:3){ 
  sample.list <- list()
  for(i in 1:1000){
    get.df <- get.info(df = df.interactions, sample.size = j, max.genotypes = 25)
    sample.list[[i]] <- cbind.data.frame(sim.number = i, get.df)
  }
  df.sample.list.sim[[j]] <- ldply(sample.list)
}

list.from.simulation <- df.sample.list.sim

## collapse simulation into a data frame
df.sample.list <- ldply(list.from.simulation) %>%
  mutate(unique.sim = interaction(sim.number, 
                                  genotype.richness, 
                                  reps.sampled, sep = "_"))

sim.info <- select(df.sample.list, unique.sim, sim.number, 
                   genotype.richness, 
                   reps.sampled)

web.data <- df.sample.list %>%
  select(unique.sim, 
         # replacing 'mean' with 'sum'
         willow_vLG,# = willow_vLG_sum, 
         willow_rG,# = willow_rG_sum,
         willow_aSG,# = willow_aSG_sum,
         willow_SG,# = willow_SG_sum,
         vLG_Platy,# = vLG_Platy_sum, 
         vLG_Mesopol,# = vLG_Mesopol_sum,
         vLG_Tory,# = vLG_Tory_sum, 
         vLG_Eulo,# = vLG_Eulo_sum,
         vLG_Mymarid,# = vLG_Mymarid_sum, 
         rG_Platy,# = rG_Platy_sum,
         rG_Mesopol,# = rG_Mesopol_sum, 
         rG_Tory,# = rG_Tory_sum, 
         rG_Eulo,# = rG_Eulo_sum, 
         rG_Lestodip,# = rG_Lestodip_sum,
         SG_Platy,# = SG_Platy_sum, 
         aSG_Tory) %>% # = aSG_Tory_sum
  gather(unique.sim, interaction_freq)

# not necessary for rarefaction simulation
#covariates <- select(df.sample.list, unique.sim, 
                     #total.gall.abund, total.gall_ptoid.abund, 
                     #total.gall.rich, total.ptoid.rich, 
                     #total.gall_ptoid.rich, total.interact.abund,
                     #total.interact.rich)

## split up variable (trophic levels of pairwise interactions) for creating bipartite webs
full.link.split <- colsplit(web.data$variable, "_", names = c("lower","upper"))
full.link.split.df <- cbind(web.data, full.link.split)

## create bipartite food webs for willow-gall and gall-ptoid interactions. This enables me to accurately calculated linkage density. 
willow.gall.list <- cast(filter(full.link.split.df, lower == "willow"), lower ~ upper | unique.sim) 
gall.ptoid.list <- cast(filter(full.link.split.df, lower != "willow"), lower ~ upper | unique.sim)  

## this function calculates weighted linkage density, generality, and vulnerability of bipartite food webs
web.measures.Gravel <- function(web.list){
  require(bipartite)
  
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
willow.gall.measures <- web.measures.Gravel(willow.gall.list)
gall.ptoid.measures <- web.measures.Gravel(gall.ptoid.list)

## join the results together and rename the column names
all.measures <- left_join(sim.info, willow.gall.measures, by = "unique.sim") %>%
  left_join(., gall.ptoid.measures, by = "unique.sim") %>%
  rename(., 
         link.density.plant_gall = link.density.x,
         vulnerability.plant_gall = vulnerability.x,
         generality.plant_gall = generality.x,
         link.density.gall_ptoid = link.density.y,
         vulnerability.gall_ptoid = vulnerability.y,
         generality.gall_ptoid = generality.y) %>%
  mutate(total_complexity = (link.density.plant_gall + link.density.gall_ptoid)/2)

# remove duplicate combinations of genotypes.
dups.pos <- which(duplicated(all.measures$genotype.combo) == TRUE)
all.measures.nodups <- all.measures[-dups.pos, ]

## save the results of the simulation as a dataframe.
write.csv(all.measures.nodups, "~/Documents/Genotype_Networks/data/food web complexity simulation TESTs.csv")

#write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity simulation 500sims 4reps.csv")
#write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity simulation 1000sims 4reps.csv")
#write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity 100 sample-based rarefaction.csv")
#write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity 1000sims 1-3reps.csv")
#write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity simulation 2500sims 4reps.csv")
#write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity simulation monocultures 1000sims 1-4reps.csv")

## how would non-additive processes affect food-web complexity?
willow_gall <- matrix(c(0,0,0,0), ncol = 4)

tot <- matrix(rep(colSums(willow_gall), NROW(willow_gall)), 
       NROW(willow_gall), byrow = TRUE)
willow_gall/tot

networklevel(willow_gall, index = c("linkage density","generality","vulnerability"), weighted = TRUE)


#write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity simulation.csv")

#trait.samp <- as.matrix(tmp.samp[ ,18:57])
#labels <- paste(rownames(tmp.samp), tmp.samp$Genotype, sep = "_")
#rownames(trait.samp) <- labels
#trait.dist <- dist(decostand(trait.samp, method = "standardize"), method = "euclidean")
#trait.dbFD <- dbFD(trait.dist, messages = F, w.abun = F,
#                  stand.x = TRUE, calc.FRic = F, calc.FGR = F, calc.FDiv = F)

## no longer using below
## Random webs ----
random.info <- function(df, sample.size, max.genotypes){
  sample.info <- list()
  for(i in 1:max.genotypes){
    df.new <- as.data.frame(swap.web(1, df[ ,2:17])) # only interactions.
    colnames(df.new) <- colnames(df)[2:17]
    df.new <- mutate(df.new,
                     Genotype = df$Genotype)
    
    tmp <- filter(df.new, 
                  Genotype %in% sample(df.new$Genotype, i)) %>%
      group_by(Genotype)
    
    tmp.samp <- sample_n(tmp, 
                         size = sample.size, 
                         replace = TRUE)
    
    
    samp.summary <- tmp.samp %>%
      ungroup() %>%
      summarise_each(funs(mean.narm = mean(., na.rm = TRUE))) %>% # 
      select(-Genotype)
    
    total.gall.abund <- rowSums(select(samp.summary, starts_with("willow"))) # take sum across samples
    total.gall_ptoid.abund <- rowSums(select(samp.summary, vLG_Platy:aSG_Tory))
    
    genotype.richness <- i
    reps.sampled <- sample.size
    
    sample.info[[i]] <- samp.summary %>%
      mutate(genotype.richness = genotype.richness,
             reps.sampled = reps.sampled,
             #functional.divergence,
             total.gall.abund = total.gall.abund,
             total.gall_ptoid.abund = total.gall_ptoid.abund)
  }
  sample.info.df <- ldply(sample.info)   
}

#, #%>%#,
#Total_Area:HCH.tremulacin.__A220nm,
#water_content:N_imputed)# %>%
# group_by(Genotype) # this will ensure that the sample_n function subsamples each genotype

#geno.avgs <- df %>%
# group_by(Genotype) %>%
#summarise_each(funs(mean.narm = mean(., na.rm = TRUE), n())) 
#geno.avgs$avg.interaction <- rowMeans(geno.avgs[ ,2:17])
#geno.avgs$interaction.rich <- rowSums(geno.avgs[ ,2:17]>0)

#plot(geno.avgs$willow_vLG_n, geno.avgs$interaction.rich)
#cor.test(geno.avgs$avg.interaction, geno.avgs$willow_vLG_n)

#df.new <- geno.avgs
#df.new <- select(df, Genotype, willow_vLG, willow_rG, willow_aSG,
#                vLG_Platy, vLG_Mesopol, vLG_Tory, rG_Tory, aSG_Tory)



