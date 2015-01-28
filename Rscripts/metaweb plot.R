## metaweb
library(vegan) # for cca function
library(igraph)
library(ggplot2)
library(grid)

## source in required files
source('~/Documents/Genotype_Networks/Rscripts/network_management_tree_level.R')
source('~/Documents/ggnet/bipartite_plot_info.R')
source('~/Documents/ggnet/tripartite_plot_info.R')
source('~/Documents/ggnet/ggnet_bipartite.R')

## create metaweb data
metaweb <- tree_interaxns_filter_add %>%
  filter(gall_contents_collapse %in% c("aSG.larv", "rG.larv", "SG.larv", "vLG.pupa", "Eulo", "Lestodip", "Mesopol", "Mymarid", "Platy", "Tory")) %>%
  dcast(gall.sp ~ gall_contents_collapse, sum, value.var = "value") %>%
  mutate(gall.surv = aSG.larv + rG.larv + SG.larv + vLG.pupa) %>%
  select(gall.sp, Eulo:Platy, Tory, gall.surv)
rownames(metaweb) <- metaweb$gall.sp

# get plot info
metaweb.info <- bipartite_plot_info(select(metaweb, Eulo:Tory), order.type = "cca")

interaction.df <- metaweb.info[[1]] %>%
  filter(Sequence == 500 | Sequence == 1) %>%
  reshape(idvar = "Group", timevar = "Sequence", direction = "wide") %>%
  mutate(weight.trans = Weight.1) # unscaled weights

nodeinfo.df <- metaweb.info[[2]] %>%
  mutate(guild = factor(c(rep("gall", 4),
                   "predator",
                   rep("Larval parasitoid", 3),
                   rep("Egg parasitoid", 2))),
         names = c("Cecidomyiid",
                   "Rab. (bud)",
                   "Iteomyia",
                   "Rab. (stem)",
                   "Lestodiplosis", "Torymus", "Eulophid", "Mesopolobus", "Platygaster", "Mymarid"))

new_theme_empty <- theme_bw()
new_theme_empty$line <- element_blank()
new_theme_empty$rect <- element_blank()
new_theme_empty$strip.text <- element_blank()
new_theme_empty$axis.text <- element_blank()
new_theme_empty$plot.title <- element_blank()
new_theme_empty$axis.title <- element_blank()
new_theme_empty$plot.margin <- structure(c(2, 2, 2, 2), unit = "lines", # c(0,0,-1,-1) # 1,1,1,1
                                         valid.unit = 3L, class = "unit")


metaweb.plot <- ggplot(interaction.df) + 
  geom_segment(data = interaction.df, aes(x = x.1, xend = x.500, y = y.1, yend = y.500),
               color = "grey",
               size = interaction.df$weight.trans/max(interaction.df$weight.trans)*25,
               alpha = 0.75)  +
  new_theme_empty + 
  geom_point(data = filter(nodeinfo.df, y > 1), aes(x = x, y = y, fill = vertex.names),
             color = "black",
             shape = 25,
             size = 30)  + 
  geom_text(data = filter(nodeinfo.df, y > 1), aes(x = x, y = y + 0.15, label = names), size = 8) +  
  geom_point(data = filter(nodeinfo.df, y == 1), aes(x = x, y = y, fill = vertex.names),
             color = "black",
             shape = 21,
             size = 30) + 
  scale_fill_brewer(palette = "Spectral", guide = "none") +
  geom_text(data = filter(nodeinfo.df, y == 1), aes(x = x, y = y - 0.15, label = names), size = 8)

#Code to override clipping
metaweb.plot
gt <- ggplot_gtable(ggplot_build(metaweb.plot))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)
