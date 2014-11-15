# source in functions
source('~/Documents/ggnet/bipartite_plot_info.R')
source('~/Documents/ggnet/ggnet_bipartite.R')


# load data
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

# make plot pretty. Save as .png, minimum width aspect ratio of 1100
gall.net.plot + 
  geom_point(data = gall.net.info[[2]], aes(x = x, y = y, fill = Module.ID, shape = Trophic), size = 12) +
  scale_fill_manual(values = cb_palette, guide = "none") +
  scale_shape_manual(values = c(23, 21), guide = "none")
  
# plotweb style
plotweb(web, method = "normal", 
        col.high = c(rep(col.mod1, 4), col.mod2, rep(col.mod3, 2), rep(col.mod4, 3), rep(col.mod5, 4)),
        col.low = c(rep(col.mod1, 8), col.mod2, rep(col.mod3, 6), rep(col.mod4, 3), rep(col.mod5, 7)))
