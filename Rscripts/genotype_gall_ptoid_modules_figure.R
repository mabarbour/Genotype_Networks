##### Setup

# source in required functions
source('~/Documents/Genotype_Networks/Rscripts/matrix_highlights.R')
source('~/Documents/Genotype_Networks/Rscripts/plotweb2_highlights.R')
source('~/Documents/ggnet/bipartite_plot_info.R')
source('~/Documents/ggnet/tripartite_plot_info.R')
source('~/Documents/ggnet/ggnet_bipartite.R')

library(reshape)
library(reshape2)
library(bipartite)
library(plyr)
library(dplyr)
library(igraph)

#  upload general molten gall network data
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

# Sequence of genotypes and galls for plotting. Parasitoid can stay in their order as is
seq.cca.geno = c("*","P","Y","Q","H","R","W","A","N","J","D","F","I","Z","T","S","K","X","G","O","B","V","U","M","L","E")
seq.cca.galls = c("SG","rG","aSG","rsLG","vLG")

# genotype-gall food web adjusted to represent density per 1000 shoots
genotype_gall_network <- dcast(total_gall_parasitoid_network, Genotype ~ gall.sp, sum)
rownames(genotype_gall_network) <- genotype_gall_network$Genotype
genotype_gall_network <- select(genotype_gall_network, -Genotype)
genotype_gall_network_density <- genotype_gall_network/shootEsts_genotype_df$shootEst.no18*1000


genotype_gall_network_density <- genotype_gall_network_density[seq.cca.geno, seq.cca.galls]

# collapsing different willow genotypes to a single species
willow_gall_network <- as.matrix(colSums(genotype_gall_network_density)[c(4,2,1,3,5)]) # reorder to be consistent with cca plotting in genotype-gall-parasitoid network
willow_gall_network <- t(willow_gall_network)
rownames(willow_gall_network) <- "Salix hookeriana"

# gall-parasitoid food web adjusted to represent density per 1000 shoots
total_gall_parasitoid_network <- dcast(total_gall_parasitoid_network, gall.sp ~ gall_contents_collapse, sum)
rownames(total_gall_parasitoid_network) <- total_gall_parasitoid_network$gall.sp
total_gall_parasitoid_network <- select(total_gall_parasitoid_network, -gall.sp)
interaxns_focus <- select(total_gall_parasitoid_network, -c(aSG.larv, Pont.surv, rG.larv, SG.larv, vLG.pupa))
interaxns_focus_density <- interaxns_focus/shootEsts_genotype_df$shootEst.no18*1000

interaxns_focus_density <- interaxns_focus_density[seq.cca.galls, ]



# gall-parasitoid interactions occuring on each genotype. Already adjusted per density of 1000 shoots
genotype_gall_parasitoid_network <- read.csv('~/Documents/Genotype_Networks/data/genotype_gall_parasitoid_network_data.csv')
rownames(genotype_gall_parasitoid_network) <- genotype_gall_parasitoid_network$X 
genotype_gall_parasitoid_network <- select(genotype_gall_parasitoid_network, -X)

module.info <- read.csv('~/Documents/Genotype_Networks/data/module_info_genotype_gall_parasitoid_network.csv')
module.info <- arrange(module.info, Trophic, Module.ID, Vertex) # reorder for easy plotting
genotypes <- filter(module.info, Trophic == "Genotype")
gall.parasitoid <- filter(module.info, Trophic == "Gall-Parasitoid")

# reorder web based on modules
web <- genotype_gall_parasitoid_network
web <- web[as.character(genotypes$Vertex), # reorder within web
           as.character(gall.parasitoid$Vertex)] # reorder within web

# Module 1 highlighted
mod1.genotype.gall <- matrix_highlighter_2(genotype_gall_network_density, interaction.names = c("A_rG","H_rG","H_vLG","J_rG","J_vLG","K_rG","K_vLG","N_rG","N_vLG","O_rG","S_rG","S_vLG","Z_rG","Z_vLG"))


mod1.gall.ptoid <- matrix_highlighter_2(interaxns_focus_density, interaction.names = c("rG_Tory","vLG_Eulo","vLG_Mymarid","vLG_Tory"), replace.values = c(14,7,1,9))


# Module 2 highlighted
mod2.genotype.gall <- matrix_highlighter_2(genotype_gall_network_density, interaction.names = c("Y_SG"))

mod2.gall.ptoid <- matrix_highlighter_2(interaxns_focus_density, interaction.names = c("SG_Platy"), replace.values = 5)

# Module 3 highlighted
mod3.genotype.gall <- matrix_highlighter_2(genotype_gall_network_density, interaction.names = c("B_vLG", "D_vLG", "G_vLG", "P_vLG", "V_vLG","X_vLG", "X_rG"))


mod3.gall.ptoid <- matrix_highlighter_2(interaxns_focus_density, interaction.names = c("rG_Mesopol","vLG_Platy"), replace.values = c(1, 38))


# Module 4 highlighted
mod4.genotype.gall <- matrix_highlighter_2(genotype_gall_network_density, interaction.names = c("Q_rG", "R_rG", "W_rG"))

mod4.gall.ptoid <- matrix_highlighter_2(interaxns_focus_density, interaction.names = c("rG_Eulo","rG_Lestodip","rG_Platy"), replace.values = c(2, 3, 2))

# Module 5 highlighted
mod5.genotype.gall <- matrix_highlighter_2(genotype_gall_network_density, interaction.names = c("*_vLG", "E_vLG", "F_aSG", "F_rsLG","F_vLG","I_aSG","I_rsLG","I_vLG","L_rsLG","L_vLG","M_rsLG","M_vLG","T_vLG"))


mod5.gall.ptoid <- matrix_highlighter_2(interaxns_focus_density, interaction.names = c("aSG_Tory","rsLG_Eury","rsLG_Lathro","vLG_Mesopol"), replace.values = c(4,3,2,17))


# melt all of the module data
mod1.low.melt <- mod1.genotype.gall %>%
  melt() %>%
  filter(value > 0) %>%
  mutate(module = rep("module 1", 14))

mod1.high.melt <- mod1.gall.ptoid %>%
  melt() %>%
  filter(value > 0) %>%
  mutate(module = rep("module 1", 4))

mod2.low.melt <- mod2.genotype.gall %>%
  melt() %>%
  filter(value > 0) %>%
  mutate(module = rep("module 2", 1))

mod2.high.melt <- mod2.gall.ptoid %>%
  melt() %>%
  filter(value > 0) %>%
  mutate(module = rep("module 2", 1))

mod3.low.melt <- mod3.genotype.gall %>%
  melt() %>%
  filter(value > 0) %>%
  mutate(module = rep("module 3", 7))

mod3.high.melt <- mod3.gall.ptoid %>%
  melt() %>%
  filter(value > 0) %>%
  mutate(module = rep("module 3", 2))


mod4.low.melt <- mod4.genotype.gall %>%
  melt() %>%
  filter(value > 0) %>%
  mutate(module = rep("module 4", 3))

mod4.high.melt <- mod4.gall.ptoid %>%
  melt() %>%
  filter(value > 0) %>%
  mutate(module = rep("module 4", 3))

mod5.low.melt <- mod5.genotype.gall %>%
  melt() %>%
  filter(value > 0) %>%
  mutate(module = rep("module 5", 13))

mod5.high.melt <- mod5.gall.ptoid %>%
  melt() %>%
  filter(value > 0) %>%
  mutate(module = rep("module 5", 4))

# bind all of molten module data frames
mods <- rbind_all(list(mod1.low.melt, mod1.high.melt, 
                       mod2.low.melt, mod2.high.melt, 
                       mod3.low.melt, mod3.high.melt, 
                       mod4.low.melt, mod4.high.melt,
                       mod5.low.melt, mod5.high.melt))
mods <- dplyr::select(mods, X1.names = X1, X2.names = X2, value, module)

# turn into graph
mods.graph <- graph.edgelist(as.matrix(mods)[,1:2])
E(mods.graph)$weight <- mods$value
#E(mods.graph)$module.id <- mods$module # not necessary right now

mods.adj <- get.adjacency(mods.graph, sparse = F, attr = "weight")
vertex.order <- c("Y", # genotypes module 2
                  "Q","R","W", # genotypes module 4
                  "A","H","J","K","N","O","S","Z", # genotypes module 1
                  "B","D","G","P","V","X", # genotypes module 3
                  "*","E","F","I","L","M","T", # genotypes module 5
                  "SG","rG","vLG","aSG","rsLG",# galls
                  "Lestodip","Eulo","Mymarid","Platy","Mesopol","Tory","Eury",
                  "Lathro") # parasitoids

mods.adj <- mods.adj[vertex.order,
                     vertex.order] 
mods.graph.ordered <- graph.adjacency(mods.adj, weighted = TRUE)

# get plot info
mods.info <- tripartite_plot_info(mods.graph.ordered, edge.list.extra = mods)
#node.info <- mods.info[[2]]
#node.info <- dplyr::select(node.info, x.val = x, y.val = y)
#ggnet_bipartite(mods.info, weighted = TRUE, edge.color = 1, edge.color.palette = c("red3","darkorchid4", "steelblue","gold","darkgreen")) + geom_text(data = node.info, aes(x.val, y.val, label = vertex.names))

# Create empty ggplot2 theme
library(ggplot2)
new_theme_empty <- theme_bw()
new_theme_empty$line <- element_blank()
new_theme_empty$rect <- element_blank()
new_theme_empty$strip.text <- element_blank()
new_theme_empty$axis.text <- element_blank()
new_theme_empty$plot.title <- element_blank()
new_theme_empty$axis.title <- element_blank()
new_theme_empty$plot.margin <- structure(c(1, 1, 1, 1), unit = "lines", # c(0,0,-1,-1) # 1,1,1,1
                                         valid.unit = 3L, class = "unit")
 
# color schemes
col.mod1 <- "red3"
col.mod2 <- "grey"
col.mod3 <- "steelblue"
col.mod4 <- "gold"
col.mod5 <- "darkorchid4"
col.galls <- "navyblue"
col.ptoids <- "black"

plot.df.sub <- mods.info[[1]] %>%
  filter(Sequence == 500 | Sequence == 1) %>%
  reshape(idvar = "Group", timevar = "Sequence", direction = "wide") %>%
  mutate(weight.trans = Weight.1) # unscaled weights

nodes.info <- mods.info[[2]]
nodes.info <- mutate(nodes.info, 
                     vertex.weights.TLscaled = c(vertex.weights[1:25]/max(vertex.weights[1:25]), vertex.weights[26:30]/max(vertex.weights[26:30]), vertex.weights[31:38]/max(vertex.weights[31:38])),
                     vertex.names = factor(as.character(vertex.names), levels = unique(as.character(vertex.names)), ordered = TRUE),
                     vertex.colors = c(rep(col.mod2,1),
                                       rep(col.mod4,3),
                                       rep(col.mod1,8),
                                       rep(col.mod3,6),
                                       rep(col.mod5,7),
                                       rep(col.galls,5),
                                       rep(col.ptoids,8)),
                     number.of.colors = factor(c(rep(1,1), # necessary for plotting node colors
                                          rep(2,3),
                                          rep(3,8),
                                          rep(4,6),
                                          rep(5,7),
                                          rep(6,5),
                                          rep(7,8))))

nodes.info1 <- mutate(nodes.info, vertex.codes = c("Y","Q","R","W","A","H","J","K","N","O","S","Z","B","D","G","P","V","X","C","E","F","I","L","M","T",# replaced * with C for aesthetic purposes
                                      "Rabdophaga-S","Rabdophaga-B","Iteomyia-L","Cecidomyiid-S","Pontania-L",
                                      "Lestodiplosis","Eulophid","Mymarid","Platygaster","Mesopolobus","Torymus","Eurytoma","Lathrostizus"))  
nodes.info2 <- mutate(nodes.info, vertex.codes = c("Y","Q","R","W","A","H","J","K","N","O","S","Z","B","D","G","P","V","X","C","E","F","I","L","M","T",# replaced * with C for aesthetic purposes
                                                  "Rab.-S","Rab.-B","Iteo.","Cec.","Pont.",
                                                  "Lest.","Eulo.","Mym.","Platy.","Meso.","Tory.","Eury.","Lath.")) 

nodes.info <- nodes.info2

nodes.info.TL1 <- filter(nodes.info, vertex.trophic.level == 1)
nodes.info.TL2 <- filter(nodes.info, vertex.trophic.level == 2)
nodes.info.TL3 <- filter(nodes.info, vertex.trophic.level == 3)

plot.df.sub$edge.colors <- rep(NA, dim(plot.df.sub)[1])
plot.df.sub$edge.colors[which(plot.df.sub$Extra.1 == "module 1")] <- col.mod1
plot.df.sub$edge.colors[which(plot.df.sub$Extra.1 == "module 2")] <- col.mod2
plot.df.sub$edge.colors[which(plot.df.sub$Extra.1 == "module 3")] <- col.mod3
plot.df.sub$edge.colors[which(plot.df.sub$Extra.1 == "module 4")] <- col.mod4
plot.df.sub$edge.colors[which(plot.df.sub$Extra.1 == "module 5")] <- col.mod5

node.size = 12

net_modules_plot <- ggplot(plot.df.sub) + 
  geom_segment(data = plot.df.sub, aes(x = x.1, xend = x.500, y = y.1, yend = y.500),
                   color = plot.df.sub$edge.colors,
                   size = plot.df.sub$weight.trans/max(plot.df.sub$weight.trans)*10,
                   alpha = 0.75) +
  new_theme_empty + 
  geom_point(data = nodes.info.TL1, aes(x = x, y = y, fill = vertex.names),
             color = "black",
             shape = 21,
             size = node.size) + # node.size
  geom_point(data = nodes.info.TL2, aes(x = x, y = y),
             color = "black",
             fill = "white",
             shape = 24,
             size = node.size*2-5) + # 
  geom_point(data = nodes.info.TL3, aes(x = x, y = y),
             color = "black",
             fill = "white",
             shape = 25,
             size = node.size*2) + # 
  scale_fill_manual(values = nodes.info$vertex.colors, guide = "none") +
  geom_text(data = nodes.info.TL2, aes(x = x, y = y, label = vertex.codes),
            fontface = "italic",
            size = 4) +
  geom_text(data = nodes.info.TL3, aes(x = x, y = y, label = vertex.codes),
            fontface = "italic",
            size = 4) +
  geom_text(data = nodes.info.TL1, aes(x = x, y = y, label = vertex.codes))

#Code to override clipping
net_modules_plot
gt <- ggplot_gtable(ggplot_build(net_modules_plot))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

ggsave("~/Documents/Genotype_Networks/figures/genotype_gall_parasitoid_modules.png", height = 8, width = 11.5, units = "in")



########## Total network plots without different genotypes
gall.ptoid.melt <- as.matrix(interaxns_focus_density) %>%
  melt() %>%
  filter(value > 0) 

willow.gall <- data.frame(X1 = rep("willow",5), X2 = c("rG","vLG","rsLG","aSG","SG"), value = rep(1,5)) # not quantitative so value is okay

total.melt <- rbind.data.frame(willow.gall, gall.ptoid.melt)
total.melt <- dplyr::select(total.melt, X1.names = X1, X2.names = X2, value)

# turn into graph
total.graph <- graph.edgelist(as.matrix(total.melt)[,1:2])
E(total.graph)$weight <- total.melt$value
#E(total.graph)$module.id <- total$module # not necessary right now

total.adj <- get.adjacency(total.graph, sparse = F, attr = "weight")
vertex.order <- c("willow",
                  "SG","rG","vLG","aSG","rsLG",# galls
                  "Lestodip","Eulo","Mymarid","Platy","Mesopol","Tory","Eury",
                  "Lathro") # parasitoids

total.adj <- total.adj[vertex.order,
                     vertex.order] 
total.graph.ordered <- graph.adjacency(total.adj, weighted = TRUE)

gplot(total.adj)

# get plot info
total.info <- tripartite_plot_info(total.graph.ordered, edge.list.extra = total.melt)

ggnet_bipartite(total.info)

# Create empty ggplot2 theme
new_theme_empty <- theme_bw()
new_theme_empty$line <- element_blank()
new_theme_empty$rect <- element_blank()
new_theme_empty$strip.text <- element_blank()
new_theme_empty$axis.text <- element_blank()
new_theme_empty$plot.title <- element_blank()
new_theme_empty$axis.title <- element_blank()
new_theme_empty$plot.margin <- structure(c(1, 1, 1, 1), unit = "lines", # c(0,0,-1,-1) # 1,1,1,1
                                         valid.unit = 3L, class = "unit")

plot.df.sub <- total.info[[1]] %>%
  filter(Sequence == 500 | Sequence == 1) %>%
  reshape(idvar = "Group", timevar = "Sequence", direction = "wide") %>%
  mutate(weight.trans = Weight.1) # unscaled weights

nodes.info <- total.info[[2]]
nodes.info <- mutate(nodes.info, 
                     vertex.weights.TLscaled = c(vertex.weights[1:25]/max(vertex.weights[1:25]), vertex.weights[26:30]/max(vertex.weights[26:30]), vertex.weights[31:38]/max(vertex.weights[31:38])),
                     vertex.names = factor(as.character(vertex.names), levels = unique(as.character(vertex.names)), ordered = TRUE),
                     vertex.colors = c(rep(col.mod2,1),
                                       rep(col.mod4,3),
                                       rep(col.mod1,8),
                                       rep(col.mod3,6),
                                       rep(col.mod5,7),
                                       rep(col.galls,5),
                                       rep(col.ptoids,8)),
                     number.of.colors = factor(c(rep(1,1), # necessary for plotting node colors
                                                 rep(2,3),
                                                 rep(3,8),
                                                 rep(4,6),
                                                 rep(5,7),
                                                 rep(6,5),
                                                 rep(7,8))))

nodes.info <- mutate(nodes.info, vertex.codes = c("Y","Q","R","W","A","H","J","K","N","O","S","Z","B","D","G","P","V","X","C","E","F","I","L","M","T",# replaced * with C for aesthetic purposes
                                                  "Rabdophaga-S","Rabdophaga-B","Iteomyia-L","Cecidomyiid-S","Pontania-L",
                                                  "Lestodiplosis","Eulophid","Mymarid","Platygaster","Mesopolobus","Torymus","Eurytoma","Lathrostizus"))  

nodes.info.TL1 <- filter(nodes.info, vertex.trophic.level == 1)
nodes.info.TL2 <- filter(nodes.info, vertex.trophic.level == 2)
nodes.info.TL3 <- filter(nodes.info, vertex.trophic.level == 3)

plot.df.sub$edge.colors <- rep(NA, dim(plot.df.sub)[1])
plot.df.sub$edge.colors[which(plot.df.sub$Extra.1 == "module 1")] <- col.mod1
plot.df.sub$edge.colors[which(plot.df.sub$Extra.1 == "module 2")] <- col.mod2
plot.df.sub$edge.colors[which(plot.df.sub$Extra.1 == "module 3")] <- col.mod3
plot.df.sub$edge.colors[which(plot.df.sub$Extra.1 == "module 4")] <- col.mod4
plot.df.sub$edge.colors[which(plot.df.sub$Extra.1 == "module 5")] <- col.mod5

node.size = 12

net_modules_plot <- ggplot(plot.df.sub) + 
  geom_segment(data = plot.df.sub, aes(x = x.1, xend = x.500, y = y.1, yend = y.500),
               color = plot.df.sub$edge.colors,
               size = plot.df.sub$weight.trans/max(plot.df.sub$weight.trans)*10,
               alpha = 0.75) +
  new_theme_empty + 
  geom_point(data = nodes.info.TL1, aes(x = x, y = y, fill = vertex.names),
             color = "black",
             shape = 21,
             size = node.size) + 
  #  geom_point(data = nodes.info.TL2, aes(x = x, y = y),
  #             color = "black",
  #             fill = "white",
  #             shape = 21,
  #             size = node.size*4-5) + 
  #  geom_point(data = nodes.info.TL3, aes(x = x, y = y),
  #             color = "black",
  #             fill = "white",
  #             shape = 21,
  #             size = node.size*3) +
  scale_fill_manual(values = nodes.info$vertex.colors, guide = "none") +
  geom_text(data = nodes.info.TL2, aes(x = x, y = y, label = vertex.codes),
            fontface = "italic",
            size = 4) +
  geom_text(data = nodes.info.TL3, aes(x = x, y = y, label = vertex.codes),
            fontface = "italic",
            size = 4) +
  geom_text(data = nodes.info.TL1, aes(x = x, y = y, label = vertex.codes))

#Code to override clipping
net_modules_plot
gt <- ggplot_gtable(ggplot_build(net_modules_plot))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

ggsave("~/Documents/Genotype_Networks/figures/genotype_gall_parasitoid_modules.png", height = 8, width = 11.5, units = "in")

# create initial plot
gall.net.plot <- ggnet_bipartite(gall.net.info)

# make plot pretty. Save as .png, minimum width aspect ratio of 1100
gall.net.plot + 
  geom_point(data = gall.net.info[[2]], aes(x = x, y = y, fill = Module.ID, shape = Trophic), size = 12) +
  scale_fill_manual(values = cb_palette, guide = "none") +
  scale_shape_manual(values = c(23, 21), guide = "none")

# plotweb style
plotweb(web, method = "normal", 
        col.high = c(rep(col.mod1, 4), col.mod2, rep(col.mod3, 2), rep(col.mod4, 3), rep(col.mod5, 4)),
        col.low = c(rep(col.mod1, 8), col.mod2, rep(col.mod3, 6), rep(col.mod4, 3), rep(col.mod5, 7)))
