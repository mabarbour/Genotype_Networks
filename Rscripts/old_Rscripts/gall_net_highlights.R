source('~/Documents/Genotype_Networks/Rscripts/matrix_highlights.R')
source('~/Documents/Genotype_Networks/Rscripts/plotweb2_highlights.R')
source('~/Documents/ggnet/bipartite_plot_info.R')
source('~/Documents/ggnet/tripartite_plot_info.R')
source('~/Documents/ggnet/ggnet_bipartite.R')

library(reshape)
library(reshape2)
library(bipartite)
library(dplyr)
library(igraph)


########### upload datasets

###  upload general molten gall network data
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

### Sequence of genotypes and galls for plotting. Parasitoid can stay in their order as is
seq.cca.geno = c("*","P","Y","Q","H","R","W","A","N","J","D","F","I","Z","T","S","K","X","G","O","B","V","U","M","L","E")
seq.cca.galls = c("SG","rG","aSG","rsLG","vLG")

### genotype-gall food web adjusted to represent density per 1000 shoots
genotype_gall_network <- dcast(total_gall_parasitoid_network, Genotype ~ gall.sp, sum)
rownames(genotype_gall_network) <- genotype_gall_network$Genotype
genotype_gall_network <- select(genotype_gall_network, -Genotype)
genotype_gall_network_density <- genotype_gall_network/shootEsts_genotype_df$shootEst.no18*1000

genotype_gall_network_density <- genotype_gall_network_density[seq.cca.geno, seq.cca.galls]

### collapsing different willow genotypes to a single species
willow_gall_network <- as.matrix(colSums(genotype_gall_network_density)[c(4,2,1,3,5)]) # reorder to be consistent with cca plotting in genotype-gall-parasitoid network
willow_gall_network <- t(willow_gall_network)
rownames(willow_gall_network) <- "Salix hookeriana"

### gall-parasitoid food web adjusted to represent density per 1000 shoots
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


mod1.hl <- matrix_highlighter(web, low.names = c("A","H","J","K","N","O","S","Z"), high.names = c("rG_Tory","vLG_Eulo","vLG_Mymarid","vLG_Tory"))


mod1.genotype.gall <- matrix_highlighter_2(genotype_gall_network_density, interaction.names = c("A_rG","H_rG","H_vLG","J_rG","J_vLG","K_rG","K_vLG","N_rG","N_vLG","O_rG","S_rG","S_vLG","Z_rG","Z_vLG"))
mod1.genotype.gall <- mod1.genotype.gall[seq.cca.geno, ]

mod1.gall.ptoid <- matrix_highlighter_2(interaxns_focus_density, interaction.names = c("rG_Tory","vLG_Eulo","vLG_Mymarid","vLG_Tory"), replace.values = colSums(mod1.hl)[1:4])
mod1.gall.ptoid <- mod1.gall.ptoid[seq.cca.galls, ]

plotweb2(genotype_gall_network_density, interaxns_focus_density, low.abun.col2 = "white")

###### Make tripartite graphs

# Module 1 highlighted
plotweb2_highlights(genotype_gall_network_density, interaxns_focus_density, highlight.web = mod1.genotype.gall, highlight.web.high = mod1.gall.ptoid, text = TRUE)

png("~/Documents/Genotype_Networks/module_1_genotype_gall_parasitoid_web.png", width = 8, height = 8, units = "in", res = 300)
plotweb2_highlights(genotype_gall_network_density, interaxns_focus_density, highlight.web = mod1.genotype.gall, highlight.web.high = mod1.gall.ptoid, text = FALSE, node.highlight = c("A","H","J","K","N","O","S","Z","rG","vLG","Tory","Mymarid","Eulo"), bord.col.highlight = "red")
dev.off()

# Module 3 highlighted
mod3.genotype.gall <- matrix_highlighter_2(genotype_gall_network_density, interaction.names = c("B_vLG", "D_vLG", "G_vLG", "P_vLG", "V_vLG","X_vLG", "X_rG"))
mod3.genotype.gall <- mod3.genotype.gall[seq.cca.geno, ]

mod3.gall.ptoid <- matrix_highlighter_2(interaxns_focus_density, interaction.names = c("rG_Mesopol","vLG_Platy"), replace.values = c(1, 38))
mod3.gall.ptoid <- mod3.gall.ptoid[seq.cca.galls, ]

plotweb2_highlights(genotype_gall_network_density, interaxns_focus_density, highlight.web = mod3.genotype.gall, highlight.web.high = mod3.gall.ptoid, text = TRUE, node.highlight = c("B","D","X","V","G","vLG","rG","Platy","Mesopol"), bord.col.highlight = "steelblue", col.interaction.highlight = "steelblue")

png("~/Documents/Genotype_Networks/module_3_genotype_gall_parasitoid_web.png", width = 8, height = 8, units = "in", res = 300)
plotweb2_highlights(genotype_gall_network_density, interaxns_focus_density, highlight.web = mod3.genotype.gall, highlight.web.high = mod3.gall.ptoid, text = FALSE, node.highlight = c("B","D","X","V","G","vLG","rG","Platy","Mesopol"), bord.col.highlight = "steelblue", col.interaction.highlight = "steelblue")
dev.off()

# Module 5 highlighted
mod5.genotype.gall <- matrix_highlighter_2(genotype_gall_network_density, interaction.names = c("C_vLG", "E_vLG", "F_aSG", "F_rsLG","F_vLG","I_aSG","I_rsLG","I_vLG","L_rsLG","L_vLG","M_rsLG","M_vLG","T_vLG"))
mod5.genotype.gall <- mod5.genotype.gall[seq.cca.geno, ]

mod5.gall.ptoid <- matrix_highlighter_2(interaxns_focus_density, interaction.names = c("aSG_Tory","rsLG_Eury","rsLG_Lathro","vLG_Mesopol"), replace.values = c(4,3,2,17))
mod5.gall.ptoid <- mod5.gall.ptoid[seq.cca.galls, ]

plotweb2_highlights(genotype_gall_network_density, interaxns_focus_density, highlight.web = mod5.genotype.gall, highlight.web.high = mod5.gall.ptoid, text = TRUE, node.highlight = c("C","E","F","I","L","M", "T", "aSG", "rsLG", "vLG","Tory","Eury", "Lathro", "Mesopol"), bord.col.highlight = "green", col.interaction.highlight = "green")

png("~/Documents/Genotype_Networks/module_5_genotype_gall_parasitoid_web.png", width = 8, height = 8, units = "in", res = 300)
plotweb2_highlights(genotype_gall_network_density, interaxns_focus_density, highlight.web = mod5.genotype.gall, highlight.web.high = mod5.gall.ptoid, text = FALSE, node.highlight = c("C","E","F","I","L","M", "T", "aSG", "rsLG", "vLG","Tory","Eury", "Lathro", "Mesopol"), bord.col.highlight = "green", col.interaction.highlight = "green")
dev.off()

plotweb2_highlights(genotype_gall_network_density, interaxns_focus_density, highlight.web = mod5.genotype.gall, highlight.web.high = mod5.gall.ptoid, text = FALSE, node.highlight = c("C","E","F","I","L","M", "T", "aSG", "rsLG", "vLG","Tory","Eury", "Lathro", "Mesopol"), bord.col.highlight = "green", col.interaction.highlight = "green", col.interaction = "white", bord.col.interaction = "white", col.interaction2 = "white", bord.col.interaction2 = "white" )

########## try using ggnet to make plots of individual modules
melt.low.whole <- mutate(genotype_gall_network_density, Genotypes = rownames(genotype_gall_network_density))
melt.low.whole <- melt(melt.low.whole, id.vars = "Genotypes") %>%
  filter(value > 0) %>%
  select(X1 = Genotypes, X2 = variable, value)

melt.up.whole <- mutate(interaxns_focus_density, Galls = rownames(interaxns_focus_density))
melt.up.whole <- melt(melt.up.whole, id.vars = "Galls") %>%
  filter(value > 0) %>%
  select(X1 = Galls, X2 = variable, value)

melt.whole <- rbind_all(list(melt.low.whole, melt.up.whole))
#net.whole <- network(melt.whole, directed = TRUE, matrix.type = "edgelist")
#ggnet(net.whole, weight.method = "degree")
#library(GGally)
#ggnet(net.whole.i, top8.nodes = TRUE)

net.whole.i <- graph.edgelist(as.matrix(melt.whole)[,1:2])
E(net.whole.i)$weight <- melt.whole$value

t <- tripartite_plot_info(net.whole.i)

ggnet_bipartite(t, weighted = TRUE, scaled = function(x) log(x))

#library(GGally)
#ggnet(net.whole.i, segment.size = log(E(net.whole.i)$weight+1))
#debug(tripartite_plot_info(net.whole.i))
# can try picking different colors and then gradients for those colors based on their weights (so gradient is also correlated with size and color corresponds to trophic level)
blues <- brewer.pal(5, "Blues")
greens <- brewer.pal(, "Greens")
greys <- brewer.pal(8, "Greys")

#t[[2]]$y <- factor(round(t[[2]]$y))
#t[[2]]$y_cont <- as.numeric(t[[2]]$y)
#genos <- filter(t[[2]], y == 1)
#quantile(genos$vertex.weights, probs = seq(0,1,0.2))
#galls <- filter(t[[2]], y == 2)
#ptoids <- filter(t[[2]], y == 3)
#with(t[[2]], interaction(y, vertex.weights))

geno.colors <- colorRampPalette(c("lightgreen","darkgreen"))(26)
gall.colors <- colorRampPalette(c("lightblue","darkblue"))(5)
ptoid.colors <- colorRampPalette(c("lightgrey","darkgrey"))(8)

extras <- t[[2]]
extras <- mutate(extras, vertex.weights.TLscaled = c(vertex.weights[1:26]/max(vertex.weights[1:26]), vertex.weights[27:31]/max(vertex.weights[27:31]), vertex.weights[32:39]/max(vertex.weights[32:39])),
                 vertex.names = factor(as.character(vertex.names), levels = unique(as.character(vertex.names)), ordered = TRUE),
                 vertex.colors = c(geno.colors, gall.colors, ptoid.colors))

ggnet_bipartite(t, range = c(1/10,1)) +
  geom_point(data = extras, aes(x = x, y = y, fill = vertex.names), shape = 21, size = extras$vertex.weights.TLscaled*15) + scale_fill_manual(values = extras$vertex.colors)# + geom_text(data = extras, aes(x = x, y = y, label = vertex.names))

ggnet_polygon(t) +
  geom_point(data = extras, aes(x = x, y = y, fill = vertex.names), shape = 21, size = extras$vertex.weights.TLscaled*15) + scale_fill_manual(values = extras$vertex.colors)# + geom_text(data = extras, aes(x = x, y = y, label = vertex.names))

ggnet_bipartite(t, weighted = TRUE)

build1 + 
  geom_point(data = ptoids, aes(x = x, y = y, fill = ptoids$vertex.weights), shape = 22, size = log(ptoids$vertex.weights)*4) + scale_fill_gradient(low = "pink", high = "red") 

test <- get.adjacency(net.whole.i, sparse = F)
test2 <- data.frame(TrophInd(test))
test2 <- mutate(test2, vertices = rownames(test2), TL = round(TL))
filter(test2, TL == "1")
GenInd(test)

net.whole.comms <- fastgreedy.community(net.whole.i)

plot.igraph(net.whole.i)
library(igraph)
melt.low.mod5 <- melt(mod5.genotype.gall) %>%
  filter(value > 0)

melt.up.mod5 <- melt(mod5.gall.ptoid) %>%
  filter(value > 0)

melt.mod5 <- rbind_all(list(melt.low.mod5, melt.up.mod5))
unmelt.mod5 <- dcast(melt.mod5, X1 ~ X2, sum)
rownames(unmelt.mod5) <- unmelt.mod5$X1
unmelt.mod5 <- select(unmelt.mod5, -X1)
net.mod5 <- network(melt.mod5, matrix.type = "edgelist", directed = TRUE)
ggnet(net.mod5, weight.method = "degree", label.nodes = unique(c(melt.mod5$X1, melt.mod5$X2)))

# willow-gall-parasitoid plot with no labels
png("~/Documents/Genotype_Networks/willow_gall_parasitoid_web.png", width = 8, height = 8, units = "in", res = 300)
plotweb2_highlights(willow_gall_network, interaxns_focus_density)#, sequence = list(seq.pred = c("SG","rG","aSG","rsLG","vLG"), seq.prey = c("Salix hookeriana")), text = TRUE)
dev.off()

# genotype-gall-parasitoid
png("~/Documents/Genotype_Networks/genotype_gall_parasitoid_web_text.png", width = 8, height = 8, units = "in", res = 300)
plotweb2_highlights(genotype_gall_network_density, interaxns_focus_density, empty = TRUE, text = TRUE) # cca organization (happens once empty = TRUE)
dev.off()

png("~/Documents/Genotype_Networks/genotype_gall_parasitoid_web.png", width = 8, height = 8, units = "in", res = 300)
plotweb2_highlights(genotype_gall_network_density, interaxns_focus_density, empty = TRUE, text = FALSE) # cca organization (happens once empty = TRUE)
dev.off()

# community genetics representation of food web
png("~/Documents/Genotype_Networks/genotype_gall_parasitoid_web_no_interactions.png", width = 8, height = 8, units = "in", res = 300)
plotweb2_highlights(genotype_gall_network_density, interaxns_focus_density, empty = TRUE, text = FALSE,
                    col.interaction2 = "white",
                    bord.col.interaction2 = "white",
                    low.abun.col2 = "black",
                    bord.low.abun.col2 = "black")
dev.off()

plotweb2_highlights(genotype_gall_network_density,
                    interaxns_focus_density, 
                    empty = FALSE,
                    highlight.web = mod1.genotype.gall, 
                    highlight.web.high = mod1.gall.ptoid, 
                    node.highlight = c("rG","vLG","A","H","J","K","N","O","S","Z","Eulo","Mymarid","Tory"),
                    col.interaction = "grey90",
                    bord.col.interaction = "black",
                    bord.col.interaction.highlight = "black",
                    col.interaction2 = "grey90",
                    bord.col.interaction2 = "black",
                    bord.col.highlight = "red", 
                    bord.size.highlight = 8,
                    y_width = 0.15,
                    text = FALSE)


### Network with pie chart of node contribution to different modules

web.df <- mutate(web, Genotypes = rownames(web))
modules <- select(genotypes, Genotypes = Vertex, Module.ID = Module.ID)

web.df <- left_join(web.df, modules)

mod1 <- web.df %>%
  filter(Module.ID == 1) %>%
  select(Module.ID, rG_Tory:vLG_Tory) %>%
  mutate(vLG = vLG_Eulo + vLG_Mymarid + vLG_Tory,
         rG = rG_Tory,
         Tory = rG_Tory + vLG_Tory,
         Mymarid = vLG_Mymarid,
         Eulo = vLG_Eulo) %>%
  select(Module.ID, vLG:Eulo) %>%
  group_by(Module.ID) %>%
  summarise_each(funs(sum)) %>%
  melt(id.vars = "Module.ID")

mod2 <-  web.df %>%
  filter(Module.ID == 2) %>%
  select(Module.ID, SG_Platy) %>%
  mutate(SG = SG_Platy,
         Platy = SG_Platy) %>%
  select(Module.ID, SG:Platy) %>%
  group_by(Module.ID) %>%
  summarise_each(funs(sum)) %>%
  melt(id.vars = "Module.ID")

mod3 <- web.df %>%
  filter(Module.ID == 3) %>%
  select(Module.ID, rG_Mesopol:vLG_Platy) %>%
  mutate(rG = rG_Mesopol,
         vLG = vLG_Platy,
         Mesopol = rG_Mesopol,
         Platy = vLG_Platy) %>%
  select(Module.ID, rG:Platy) %>%
  group_by(Module.ID) %>%
  summarise_each(funs(sum)) %>%
  melt(id.vars = "Module.ID")

mod4 <- web.df %>%
  filter(Module.ID == 4) %>%
  select(Module.ID, rG_Eulo:rG_Platy) %>%
  mutate(rG = rG_Eulo + rG_Lestodip + rG_Platy,
         Eulo = rG_Eulo,
         Lestodip = rG_Lestodip,
         Platy = rG_Platy) %>%
  select(Module.ID, rG:Platy) %>%
  group_by(Module.ID) %>%
  summarise_each(funs(sum)) %>%
  melt(id.vars = "Module.ID")

mod5 <- web.df %>%
  filter(Module.ID == 5) %>%
  select(Module.ID, aSG_Tory:vLG_Mesopol) %>%
  mutate(aSG = aSG_Tory,
         rsLG = rsLG_Eury + rsLG_Lathro,
         vLG = vLG_Mesopol,
         Tory = aSG_Tory,
         Eury = rsLG_Eury,
         Lathro = rsLG_Lathro,
         Mesopol = vLG_Mesopol) %>%
  select(Module.ID, aSG:Mesopol) %>%
  group_by(Module.ID) %>%
  summarise_each(funs(sum)) %>%
  melt(id.vars = "Module.ID")

node.module.partitions <- rbind_all(list(mod1, mod2, mod3, mod4, mod5))

library(ggplot2)
theme_pie <- theme_classic() + theme(legend.position = "none",
                                     axis.line = element_blank(),
                                     axis.text = element_blank(),
                                     axis.title = element_blank(),
                                     axis.ticks = element_blank())

p.vLG <- ggplot(data = filter(node.module.partitions, variable == "vLG"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity") + coord_polar(theta = "y") + theme_pie

p.rG <- ggplot(data = filter(node.module.partitions, variable == "rG"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity") + coord_polar(theta = "y") + theme_pie

p.SG <- ggplot(data = filter(node.module.partitions, variable == "SG"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity") + coord_polar(theta = "y") + theme_pie



web.df.sum <- web.df %>%
  group_by(Module.ID) %>%
  select(-Genotypes) %>%
  summarise_each(funs(sum)) %>%
  melt(id.vars = "Module.ID") %>%
  mutate(rGs = starts_with("rG"))
  
  mutate(vLG_1 = vLG_Eulo + vLG_Mymarid + vLG_Tory,
            rG_1 = rG_Tory,
            Tory_1 = rG_Tory + vLG_Tory,
            Mymarid_1 = vLG_Mymarid,
            Eulo_1 = vLG_Eulo,
            SG_2 = SG_Platy,
            Platy_2 = SG_Platy,
            rG_3 = rG_Mesopol,
            vLG_3 = vLG_Platy,
            Mesopol_3 = rG_Mesopol,
            Platy_3 = vLG_Platy,
            rG_4 = rG_Eulo + rG_Lestodip + rG_Platy,
            Eulo_4 = rG_Eulo,
            Lestodip_4 = rG_Lestodip,
            Platy_4 = rG_Platy,
            aSG_5 = aSG_Tory,
            rsLG_5 = rsLG_Eury + rsLG_Lathro,
            vLG_5 = vLG_Mesopol,
            Tory_5 = aSG_Tory,
            Eury_5 = rsLG_Eury,
            Lathro_5 = rsLG_Lathro,
            Mesopol_5 = vLG_Mesopol) %>%
  select(Module.ID:Mesopol_5) %>%
  group_by(Module.ID) %>%
  summarise_each(funs(sum)) # note that this dataset include interactions outside of these modules as well

melt(web.df.sum, id.vars = "Module.ID")
