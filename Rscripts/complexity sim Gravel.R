# In order to get an expectation of how my data is behaving, I could randomly shuffle everything and see what patterns I observe. Alternatively, I could create different scenarios to see what is going on (ref notebook). Let's try shuffling things up first...

# PERHAPS BOOTSTRAPPING IS THE WAY TO GO. NEED TO THINK ABOUT THE ASSUMPTIONS BEHIND THIS A BIT MORE.
# bootstrapping relies on one crucial assumption: the data accurately represents the true population...
# add a calculation of functional diversity and see if this has a stronger relationship with food-web complexity.

# Gravel simulation

source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')

#require(mvabund)
require(bipartite)
require(tidyr)
require(FD)

df <- tree_level_interaxn_all_plants_traits_size %>%
  tbl_df() %>%
  filter(Genotype != "U") %>% # never collected any galls and thus gall-ptoid interactions so I removed it from the dataset. 
  select(Genotype, 
         willow_vLG = vLG_abund,
         willow_rG = rG_abund,
         willow_aSG = aSG_abund,
         willow_SG = SG_abund,
         vLG_Platy, vLG_Mesopol, vLG_Tory, vLG_Eulo, vLG_Mymarid,
         rG_Platy, rG_Mesopol, rG_Tory, rG_Eulo, rG_Lestodip,
         SG_Platy, aSG_Tory)#, #%>%#,
         #Total_Area:HCH.tremulacin.__A220nm,
         #water_content:N_imputed)# %>%
 # group_by(Genotype) # this will ensure that the sample_n function subsamples each genotype

geno.avgs <- df %>%
  group_by(Genotype) %>%
  summarise_each(funs(mean.narm = mean(., na.rm = TRUE), n())) 
geno.avgs$avg.interaction <- rowMeans(geno.avgs[ ,2:17])
geno.avgs$interaction.rich <- rowSums(geno.avgs[ ,2:17]>0)

plot(geno.avgs$willow_vLG_n, geno.avgs$interaction.rich)
cor.test(geno.avgs$avg.interaction, geno.avgs$willow_vLG_n)

#df.new <- geno.avgs
df.new <- select(df, Genotype, willow_vLG, willow_rG, willow_aSG,
                 vLG_Platy, vLG_Mesopol, vLG_Tory, rG_Tory, aSG_Tory)

#df <- df[complete.cases(df), ]

#sample.size <- 4
#max.genotypes <- 25

## Gravel sample function 1 ----
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
    
    #trait.samp <- as.matrix(tmp.samp[ ,18:57])
    #labels <- paste(rownames(tmp.samp), tmp.samp$Genotype, sep = "_")
    #rownames(trait.samp) <- labels
    #trait.dist <- dist(decostand(trait.samp, method = "standardize"), method = "euclidean")
    #trait.dbFD <- dbFD(trait.dist, messages = F, w.abun = F,
     #                  stand.x = TRUE, calc.FRic = F, calc.FGR = F, calc.FDiv = F)
    
    ## calculate the average frequency of each willow-gall and gall-parasitoid interaction
    samp.summary <- tmp.samp %>%
      ungroup() %>%
      summarise_each(funs(mean)) %>% # try without mean.narm
      select(-Genotype)

    total.gall.abund <- rowSums(select(samp.summary, starts_with("willow"))) 
    total.gall_ptoid.abund <- rowSums(select(samp.summary, vLG_Platy:aSG_Tory))
    
    #functional.divergence <- trait.dbFD$FEve
    
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
    
    #trait.samp <- as.matrix(tmp.samp[ ,18:57])
    #labels <- paste(rownames(tmp.samp), tmp.samp$Genotype, sep = "_")
    #rownames(trait.samp) <- labels
    #trait.dist <- dist(decostand(trait.samp, method = "standardize"), method = "euclidean")
    #trait.dbFD <- dbFD(trait.dist, messages = F, w.abun = F,
    #                  stand.x = TRUE, calc.FRic = F, calc.FGR = F, calc.FDiv = F)
    
    samp.summary <- tmp.samp %>%
      ungroup() %>%
      summarise_each(funs(mean.narm = mean(., na.rm = TRUE))) %>% # 
      select(-Genotype)
    
    total.gall.abund <- rowSums(select(samp.summary, starts_with("willow"))) # take sum across samples
    total.gall_ptoid.abund <- rowSums(select(samp.summary, vLG_Platy:aSG_Tory))
    
    #functional.divergence <- trait.dbFD$FEve
    
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

## run simulation ----
df.sample.list <- list()
for(j in 1:5){ #c(4,12,16,20)
  sample.list <- list()
  for(i in 1:100){
    #df.new <- as.data.frame(shuffle.web(1, df[ ,2:17])) # only interactions.  first used swap.web, try r2dtable
    #colnames(df.new) <- colnames(df)[2:17]
    #df.new <- mutate(df.new,
                  #   Genotype = df$Genotype)
    #get.df <- random.info(df, sample.size = j, max.genotypes)
    get.df <- get.info(df.new, sample.size = j, max.genotypes)
    sample.list[[i]] <- cbind.data.frame(sim.number = i, get.df)
  }
  df.sample.list[[j]] <- ldply(sample.list)
}

#sample.list <- list()
#for(i in 1:10){
 # get.df <- get.info(df, sample.size, max.genotypes)
  #sample.list[[i]] <- cbind.data.frame(sim.number = i, get.df)
#}
df.sample.list <- ldply(df.sample.list) %>%
  mutate(unique.sim = interaction(sim.number, genotype.richness, reps.sampled, sep = "_"))

sim.info <- select(df.sample.list, unique.sim, sim.number, genotype.richness, reps.sampled)

web.data <- df.sample.list %>%
  select(unique.sim, willow_vLG:aSG_Tory) %>%
  gather(unique.sim, interaction_freq)

covariates <- select(df.sample.list, unique.sim, 
                     #functional.divergence, 
                     total.gall.abund, total.gall_ptoid.abund)

## split up variable (trophic levels of pairwise interactions) for creating bipartite webs
full.link.split <- colsplit(web.data$variable, "_", names = c("lower","upper"))
full.link.split.df <- cbind(web.data, full.link.split)
#full.link.split.df$combo <- with(full.link.split.df, interaction(factor(replicate), factor(genotype.richness), sep = "_"))

## create bipartite food webs for willow-gall and gall-ptoid interactions. This enables me to accurately calculated linkage density. 
willow.gall.list <- cast(filter(full.link.split.df, lower == "willow"), lower ~ upper | unique.sim) 
gall.ptoid.list <- cast(filter(full.link.split.df, lower != "willow"), lower ~ upper | unique.sim) # 

## this function calculates weighted linkage density, generality, and vulnerability of bipartite food webs
web.measures.Gravel <- function(web.list){
  require(bipartite)
  
  web.list.measures <- list()
  unique.sim <- c()

  for(i in 1:dim(web.list)){ 
    unique.sim[i] <- names(web.list[i])

    tmp.web <- web.list[[i]]
    rownames(tmp.web) <- tmp.web$lower
    tmp.web <- as.matrix.data.frame(tmp.web[ ,-1])
    
    if(sum(tmp.web, na.rm = TRUE) > 0){
      web.list.measures[[i]] <- networklevel(tmp.web, index = c("linkage density","generality","vulnerability"), weighted = TRUE)
    } else {
      web.list.measures[[i]] <- c(NA,NA,NA)
    } 
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
  left_join(., covariates, by = "unique.sim")

colnames(all.measures)[5:10] <- c("link.density.plant_gall", "generality.HL.plant_gall", "vulnerability.LL.plant_gall", "link.density.gall_ptoid", "generality.HL.gall_ptoid", "vulnerability.LL.gall_ptoid")

## calculating total complexity in multiple ways gives the same results.
all.measures <- mutate(all.measures,
                       total_complexity = (link.density.plant_gall + link.density.gall_ptoid)/2,
                       total_complexity.check = (generality.HL.plant_gall + vulnerability.LL.plant_gall + generality.HL.gall_ptoid + vulnerability.LL.gall_ptoid)/4,
                       plants_sampled = genotype.richness*reps.sampled)

#hist(all.measures$genotype.richness)
#hist(all.measures$plants_sampled)

#plot(total_complexity ~ genotype.richness, all.measures)

#library(visreg)
#visreg(lm(total_complexity ~ plants_sampled + genotype.richness, all.measures))

#summary(lm(total_complexity ~ total.gall.abund + genotype.richness, filter(all.measures, total.gall.abund < 1000)))
#visreg(lm(total_complexity ~ total.gall.abund, filter(all.measures, total.gall.abund < 1000)))


## save the results of the simulation as a dataframe.
write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity simulation Gravel.csv")
#write.csv(all.measures, "~/Documents/Genotype_Networks/data/food web complexity simulation.csv")


