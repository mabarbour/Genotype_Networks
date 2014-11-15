# Total Gall-Parasitoid Network Analysis
library(reshape2)
library(reshape)
library(dplyr)
library(bipartite)

source('~/Documents/ggnet/bipartite_plot_info.R')
source('~/Documents/ggnet/ggnet_bipartite.R')

# upload molten gall network data
gall_net_melt <- read.csv("~/Documents/Genotype_Networks/data/gall_network_data.csv")
gall_net_melt <- tbl_df(gall_net_melt)

total_gall_parasitoid_network <- gall_net_melt %>%
  filter(gall_contents %in% c("aSG.larv", "Pont.ad", "Pont.prep", "rG.larv", "SG.larv", "vLG.pupa", "Eulo.fem", "Eulo.mal", "Eury.fem", "Eury.mal", "Lathro.fem", "Lathro.mal", "Lestodip", "Mesopol", "Mymarid", "Platy", "Ptero.2", "Tory.fem", "Tory.mal")) 

total_gall_parasitoid_network <- mutate(total_gall_parasitoid_network, gall_contents_collapse = revalue(gall_contents, c("Pont.ad" = "Pont.surv", "Pont.prep" = "Pont.surv", "Eulo.fem" = "Eulo", "Eulo.mal" = "Eulo", "Eury.fem" = "Eury", "Eury.mal" = "Eury", "Lathro.fem" = "Lathro", "Lathro.mal" = "Lathro", "Ptero.2" = "Mesopol", "Tory.fem" = "Tory", "Tory.mal" = "Tory") ))

# shoot estimates on a genotype level
shootEsts_df <- gall_net_melt %>%
  mutate(plant.position = as.factor(plant.position)) %>%
  group_by(Gender, Genotype, plant.position) %>%
  summarise(shootEst.no18 = mean(shootEst.no18), shootEst.all = mean(shootEst.all), galls.found = mean(Galls.found)) # needed to take the mean because these shoot estimates are duplicated throughout the molten dataframe
shootEsts_df <- mutate(shootEsts_df, plant.position = as.character(plant.position))

shootEsts_genotype_df <- shootEsts_df %>%
  group_by(Genotype) %>%
  summarise(n = n(), shootEst.no18 = sum(shootEst.no18), shootEst.all = sum(shootEst.all))

genotype_gall_network <- dcast(total_gall_parasitoid_network, Genotype ~ gall.sp, sum)

rownames(genotype_gall_network) <- genotype_gall_network$Genotype
genotype_gall_network <- select(genotype_gall_network, -Genotype)
genotype_gall_network_density <- genotype_gall_network/shootEsts_genotype_df$shootEst.no18*1000
interaxns_focus_density <- interaxns_focus/shootEsts_genotype_df$shootEst.no18*1000

low.abun = round(runif(dim(genotype_gall_network_density)[1],1,40)) #create
names(low.abun) <- rownames(genotype_gall_network_density)
plotweb(genotype_gall_network_density, text.rot=90, high.abun=low.abun, col.interaction="purple", 
        y.lim=c(0,4.5), high.lablength=0, arrow="up", method="normal", 
        y.width.high=0.05)

plotweb(interaxns_focus_density, y.width.low=0.05, y.width.high=0.1, method="normal", 
        add=TRUE, low.y=1.7,high.y=2.7, col.low="green", text.low.col="black", 
        low.lab.dis=0, arrow="down", adj.low=c(0.5,1.1), low.lablength=4, 
        high.lablength=0)

plotweb(Safariland,y.width.low=0.05, y.width.high=0.1, method="normal", 
        add=TRUE, low.y=2.95, high.y=3.95, col.low="green", text.low.col="black", 
        low.lab.dis=0, arrow="down", adj.low=c(0.5,1.1), low.lablength=4)


plotweb(genotype_gall_network_density, col.interaction = cb_palette, method = "normal")
plotweb(interaxns_focus_density, col.interaction = cb_palette, add = TRUE, method = "normal", low.y = 1.7, high.y = 2.7)
plotweb(genotype_gall_network, low.abun = shootEsts_genotype_df$shootEst.no18/10 - rowSums(genotype_gall_network), abuns.type = "additional", low.abun.col = "white")

plotweb2(genotype_gall_network_density, interaxns_focus_density, arrow = "no", low.abun.col2 = "white", col.interaction = c("blue"), col.interaction2 = "red")
plotweb2(genotype_gall_network, interaxns_focus, low.abun.col2 = "white", low.abun = shootEsts_genotype_df$shootEst.no18/10 - rowSums(genotype_gall_network), low.abun.col = "white")
total_gall_parasitoid_network <- dcast(total_gall_parasitoid_network, gall.sp ~ gall_contents_collapse, sum)

rownames(total_gall_parasitoid_network) <- total_gall_parasitoid_network$gall.sp
total_gall_parasitoid_network <- select(total_gall_parasitoid_network, -gall.sp)

interaxns_focus <- select(total_gall_parasitoid_network, -c(aSG.larv, Pont.surv, rG.larv, SG.larv, vLG.pupa))

surviving.larva <- with(total_gall_parasitoid_network, aSG.larv + Pont.surv + rG.larv + SG.larv + vLG.pupa)

# bipartite visualizations
library(png)
vLG <- readPNG("~/Desktop/figure_genotype_gall_network/Gall PNGs for Mateo/iteomyia gall (red volcano).png")
rG <-  readPNG("~/Desktop/figure_genotype_gall_network/Gall PNGs for Mateo/Rosette Gall.png")
rsLG <-  readPNG("~/Desktop/figure_genotype_gall_network/Gall PNGs for Mateo/Pontania Gall.png")

vLG <- rasterGrob(vLG, interpolate=TRUE)
rG <- rasterGrob(rG, interpolate=TRUE)
rsLG <- rasterGrob(rsLG, interpolate=TRUE)

cb_palette <- c("#56B4E9", "#999999", "#E69F00", "#009E73", "#CC79A7")

visweb(interaxns_focus, type = "diagonal")
plotweb(interaxns_focus, col.low = "black", col.high = "steelblue", arrow = "down", col.interaction = "grey")
plotweb(interaxns_focus, low.abun = surviving.larva - rowSums(interaxns_focus), abuns.type = "additional")

image.cex = 0.75
total.info <- bipartite_plot_info(interaxns_focus, order.type = "cca")
ggnet_bipartite(total.info, range = c(1,30))
ggnet_bipartite(total.info) + 
  annotation_custom(vLG, xmin = 2.75 - image.cex, xmax = 2.75 + image.cex, ymin = 0, ymax = 2) + 
  annotation_custom(rG, xmin = 4.5 - image.cex, xmax = 4.5 + image.cex, 0, 2)

hist(as.matrix(interaxns_focus[interaxns_focus > 0])) 
hist(log(as.matrix(interaxns_focus[interaxns_focus > 0]))) # degree distribution

dfun(interaxns_focus)
H2fun(interaxns_focus)

source('~/Documents/Genotype_Networks/Rscripts//best_QuaBiMo_function.R')
source('~/Documents/Genotype_Networks/Rscripts//null_model_analysis_WNODF_QuaBiMo.R')

best_partitions <- best_QuaBiMo(interaxns_focus, QuaBiMo_reps = 100) # 0.278704
plotModuleWeb(best_partitions$best_module_info[[1]])
QuaBiMo_Q = best_partitions$best_observed_value
N_null_webs = 1000

z_QuaBiMo_swap.web <- null_model_analysis_WNODF_QuaBiMo(web = interaxns_focus, 
                                                        observed_value = QuaBiMo_Q, 
                                                        type = "QuaBiMo",
                                                        N_null_webs = N_null_webs, 
                                                        null_model = "swap.web") # z = 0.074, P = 0.472


### experiment
plotweb(genotype_gall_network_density)
test.abund = c(rep(0,6),16,0)
web.info <- plotweb2_highlights(genotype_gall_network_density, method = "cca", interaxns_focus_density, col.interaction = c("red",rep("white",12),"red"), col.interaction2 = c("red"), col.prey = c(rep("white",23),"red"), col.pred = c(rep("white",4),"red"), col.pred2 = c(rep("white",6),"red","white"), bord.col.prey = c("black","black","red",rep("black",25)), high.abun2 = test.abund - colSums(interaxns_focus_density), high.abun.col2 = "green")
dim(web.info[[1]])
dim(web.info[[2]])
