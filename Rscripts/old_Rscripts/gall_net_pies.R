# source in functions
source('~/Documents/ggnet/bipartite_plot_info.R')
source('~/Documents/ggnet/ggnet_bipartite.R')


#### load data

#  upload general molten gall network data
gall_net_melt <- read.csv("~/Documents/Genotype_Networks/data/gall_network_data.csv")
gall_net_melt <- tbl_df(gall_net_melt)

total_gall_parasitoid_network <- gall_net_melt %>%
  filter(gall_contents %in% c("aSG.larv", "Pont.ad", "Pont.prep", "rG.larv", "SG.larv", "vLG.pupa", "Eulo.fem", "Eulo.mal", "Eury.fem", "Eury.mal", "Lathro.fem", "Lathro.mal", "Lestodip", "Mesopol", "Mymarid", "Platy", "Ptero.2", "Tory.fem", "Tory.mal")) 

total_gall_parasitoid_network <- mutate(total_gall_parasitoid_network, gall_contents_collapse = revalue(gall_contents, c("Pont.ad" = "Pont.surv", "Pont.prep" = "Pont.surv", "Eulo.fem" = "Eulo", "Eulo.mal" = "Eulo", "Eury.fem" = "Eury", "Eury.mal" = "Eury", "Lathro.fem" = "Lathro", "Lathro.mal" = "Lathro", "Ptero.2" = "Mesopol", "Tory.fem" = "Tory", "Tory.mal" = "Tory") ))

# focus on gall-parasitoids
total_gall_parasitoid_network <- dcast(total_gall_parasitoid_network, gall.sp ~ gall_contents_collapse, sum)
rownames(total_gall_parasitoid_network) <- total_gall_parasitoid_network$gall.sp
total_gall_parasitoid_network <- select(total_gall_parasitoid_network, -gall.sp)
interaxns_focus <- select(total_gall_parasitoid_network, -c(aSG.larv, Pont.surv, rG.larv, SG.larv, vLG.pupa))
interaxns_focus_density <- interaxns_focus/shootEsts_genotype_df$shootEst.no18*1000

interaxns_focus_qual <- interaxns_focus_density
interaxns_focus_qual[interaxns_focus_qual > 0] <- 1

gall_ptoid_qual <- bipartite_plot_info(interaxns_focus_qual, order.type = "cca")
p.net <- ggnet_bipartite(gall_ptoid_qual, range = c(1/10,1))
p.net_text <- p.net + geom_text(data = gall_ptoid_qual[[2]], aes(x = x, y = y, label = vertex.names), size = 2)
ggsave("~/Documents/Genotype_Networks/figures/p.net_text.png", height = 4, width = 6, units = "in")


# color blind friednly color palette
cb_palette <- c("#56B4E9", "#999999", "#E69F00", "#009E73", "#CC79A7")
col.mod1 <- "#56B4E9"
col.mod2 <- "#999999"
col.mod3 <- "#E69F00"
col.mod4 <- "#009E73"
col.mod5 <- "#CC79A7"

#############

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

node.module.sums <- node.module.partitions %>%
  group_by(variable) %>%
  summarise(sum(log(value+1)))

library(ggplot2)
theme_pie <- theme_classic() + theme(legend.position = "none",
                                     plot.background = element_rect(fill = "transparent", color = NA),
                                     panel.background = element_rect(fill = "transparent", color = NA),
                                     axis.line = element_blank(),
                                     axis.text = element_blank(),
                                     axis.title = element_blank(),
                                     axis.ticks = element_blank(),
                                     plot.margin = unit(c(0,0,0,0), "in"))

# gall piecharts
p.vLG <- ggplot(data = filter(node.module.partitions, variable == "vLG"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_pie + scale_fill_manual(values = c(col.mod1, col.mod3, col.mod5))
ggsave("~/Documents/Genotype_Networks/figures/pie_vLG.png", bg = "transparent")

p.rG <- ggplot(data = filter(node.module.partitions, variable == "rG"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_pie + scale_fill_manual(values = c(col.mod1, col.mod3, col.mod4))
ggsave("~/Documents/Genotype_Networks/figures/pie_rG.png", bg = "transparent")

p.SG <- ggplot(data = filter(node.module.partitions, variable == "SG"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_pie + scale_fill_manual(values = c(col.mod2))
ggsave("~/Documents/Genotype_Networks/figures/pie_SG.png", bg = "transparent")

p.aSG <- ggplot(data = filter(node.module.partitions, variable == "aSG"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_pie + scale_fill_manual(values = c(col.mod5))
ggsave("~/Documents/Genotype_Networks/figures/pie_aSG.png", bg = "transparent")

p.rsLG <- ggplot(data = filter(node.module.partitions, variable == "rsLG"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_pie + scale_fill_manual(values = c(col.mod5))
ggsave("~/Documents/Genotype_Networks/figures/pie_rsLG.png", bg = "transparent")
  
# parasitoid piecharts
p.Platy <- ggplot(data = filter(node.module.partitions, variable == "Platy"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_pie + scale_fill_manual(values = c(col.mod2, col.mod3, col.mod4))
ggsave("~/Documents/Genotype_Networks/figures/pie_Platy.png", bg = "transparent")

p.Mesopol <- ggplot(data = filter(node.module.partitions, variable == "Mesopol"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_pie + scale_fill_manual(values = c(col.mod3, col.mod5))
ggsave("~/Documents/Genotype_Networks/figures/pie_Mesopol.png", bg = "transparent")

p.Tory <- ggplot(data = filter(node.module.partitions, variable == "Tory"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_pie + scale_fill_manual(values = c(col.mod1, col.mod5))
ggsave("~/Documents/Genotype_Networks/figures/pie_Tory.png", bg = "transparent")

p.Eulo <- ggplot(data = filter(node.module.partitions, variable == "Eulo"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_pie + scale_fill_manual(values = c(col.mod1, col.mod4))
ggsave("~/Documents/Genotype_Networks/figures/pie_Eulo.png", bg = "transparent")

p.Mymarid <- ggplot(data = filter(node.module.partitions, variable == "Mymarid"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_pie + scale_fill_manual(values = c(col.mod1))
ggsave("~/Documents/Genotype_Networks/figures/pie_Mymarid.png", bg = "transparent")

p.Lestodip <- ggplot(data = filter(node.module.partitions, variable == "Lestodip"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_pie + scale_fill_manual(values = c(col.mod4))
ggsave("~/Documents/Genotype_Networks/figures/pie_Lestodip.png", bg = "transparent")

p.Eury <- ggplot(data = filter(node.module.partitions, variable == "Eury"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_pie + scale_fill_manual(values = c(col.mod5))
ggsave("~/Documents/Genotype_Networks/figures/pie_Eury.png", bg = "transparent")

p.Lathro <- ggplot(data = filter(node.module.partitions, variable == "Lathro"), aes(x = factor(1), y = value, fill = factor(Module.ID))) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_pie + scale_fill_manual(values = c(col.mod5))
ggsave("~/Documents/Genotype_Networks/figures/pie_Lathro.png", bg = "transparent")

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

gall.net.info <- bipartite_plot_info(web = web, order.type = "normal")

# add some important node attributes
gall.net.info[[2]]$Module.ID <- factor(c(genotypes$Module.ID, gall.parasitoid$Module.ID))
gall.net.info[[2]]$Trophic <- factor(c(genotypes$Trophic, gall.parasitoid$Trophic))

# color blind friednly color palette
cb_palette <- c("#56B4E9", "#999999", "#E69F00", "#009E73", "#CC79A7")
col.mod1 <- "#56B4E9"
col.mod2 <- "#999999"
col.mod3 <- "#E69F00"
col.mod4 <- "#009E73"
col.mod5 <- "#CC79A7"

# create initial plot
gall.net.plot <- ggnet_bipartite(gall.net.info)