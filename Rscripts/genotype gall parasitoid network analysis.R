####################################################################################
# This R script analyzes the willow genotype-gall-parasitoid network
# Coded by: Matt Barbour
####################################################################################

##### Upload source files
source('~/Documents/vegan/R/nestednodf.R')
source('~/Documents/Genotype_Networks/Rscripts/null_model_analysis_WNODF_QuaBiMo.R')
source('~/Documents/Genotype_Networks/Rscripts/best_QuaBiMo_function.R')
source("~/Documents/Genotype_Networks/Rscripts/module.contribution.R")
source("~/Documents/Genotype_Networks/Rscripts/null_model_analysis_czvalues.R")
source("~/Documents/Genotype_Networks/Rscripts/czvalues_error.R")
source("~/Documents/Genotype_Networks/Rscripts/teammates.R")
source("~/Documents/Genotype_Networks/Rscripts/teammates_table.R")
source("~/Documents/miscellaneous_R/ggplot_themes.R")

##### Upload required libraries
library(plyr)
library(dplyr)
library(bipartite)
library(ggplot2)
library(gridExtra)

##### Upload required datasets
genotype_gall_parasitoid_network <- read.csv('~/Documents/Genotype_Networks/data/genotype_gall_parasitoid_network_data.csv')
rownames(genotype_gall_parasitoid_network) <- genotype_gall_parasitoid_network$X 
genotype_gall_parasitoid_network <- select(genotype_gall_parasitoid_network, -X)

#### Summary statistics
colSums(genotype_gall_parasitoid_network)/sum(colSums(genotype_gall_parasitoid_network)) # aSG_Tory, rG_Eulo, rG_Lestodip, rG_Mesopol, rG_Platy, SG_Platy, and vLG_Mymarid each make up less than 5% of the total interactions.

rowSums(genotype_gall_parasitoid_network)/sum(rowSums(genotype_gall_parasitoid_network)) # Genotypes G, P, H, Q, O, N, E, M, A, and J each make up less than 2% of the total interactions in the network

##### Null Model Analysis. Used 1000 reps of the 'swap_web' algorithm for both WNODF and QuaBiMo.

### Nestedness: weighted NODF (Almeida-Neto & Ulrich, 2010, Environmental Modelling Software).
# WNODF of network = 23.7
# Null models: z-score = -1.44, mean WNODF = 28.6, SD WNODF = 3.4
WNODF_info <- nestednodf(genotype_gall_parasitoid_network, order = T, weighted = T) # 23.69687
WNODF = WNODF_info$statistic[3] # 23.7
N_null_webs = 1000

z_WNODF_swap.web <- null_model_analysis_WNODF_QuaBiMo(web = genotype_gall_parasitoid_network, 
                                                      observed_value = WNODF,
                                                      type = "WNODF",
                                                      N_null_webs = N_null_webs, 
                                                      null_model = "swap.web")

### Modularity: QuaBiMo algorithm (Dormann & Strauss 2014, Methods in Ecology and Evolution)
# Note that this null model will take a couple of hours to run for 1000 iterations.
# Using the default number of steps, but searching for the best partition at least 100 times results in a stable modularity configuration. Boosting the number of steps to 1e10 does not seem to improve this detection, and appears to actually be worse than the method I'm using.
# Q of network = 0.33
# Null models: z-score = 2.41, mean Q = 0.28, SD Q = 0.02

best_partitions <- best_QuaBiMo(genotype_gall_parasitoid_network, QuaBiMo_reps = 100) # 0.327593
plotModuleWeb(best_partitions$best_module_info[[1]])
QuaBiMo_Q = best_partitions$best_observed_value
N_null_webs = 1000

# Null model takes a while to run which is why it is temporarily coded out
# z_QuaBiMo_swap.web <- null_model_analysis_WNODF_QuaBiMo(web = genotype_gall_parasitoid_network, 
#                                                        observed_value = QuaBiMo_Q, 
#                                                        type = "QuaBiMo",
#                                                        N_null_webs = N_null_webs, 
#                                                        null_model = "swap.web")

##### Roles of nodes in the network. 
# Note that instead of using the z-value proposed by Guimera & Amaral 2005, which was originally developed for qualitative, one-mode networks, I calculated the 'percent module contribution'. This is an intuitive measure and is simply the percentage of weighted interactions occuring in a module that are contributed by a particular node. I have retained the same calculation of 'c' since it makes sense and is similarly based off the percentage of a nodes links to other modules in the network.

Qcz_values.list <- czvalues_error(web = genotype_gall_parasitoid_network, QuaBiMo_reps = 100)

mean(Qcz_values.list$Q_values) # mean Q-value for 100 reps is 0.32

colMeans(Qcz_values.list$c_higher)
apply(Qcz_values.list$c_higher, 2, function(x) sd(x, na.rm = TRUE))

colMeans(Qcz_values.list$z_higher, na.rm = TRUE)
mean(colMeans(Qcz_values.list$z_higher, na.rm = TRUE))
range(colMeans(Qcz_values.list$z_higher, na.rm = TRUE))
apply(Qcz_values.list$z_higher, 2, function(x) sd(x, na.rm = TRUE))

plot(colMeans(Qcz_values.list$z_higher, na.rm = TRUE) ~ colMeans(Qcz_values.list$c_higher), xlab = "c", ylab = "z", type = "n")
text(y = colMeans(Qcz_values.list$z_higher, na.rm = TRUE), x = colMeans(Qcz_values.list$c_higher), labels = colnames(Qcz_values.list$z_higher))

colMeans(Qcz_values.list$c_lower)
colMeans(Qcz_values.list$z_lower, na.rm = TRUE)
mean(colMeans(Qcz_values.list$z_lower, na.rm = TRUE))
median(colMeans(Qcz_values.list$z_lower, na.rm = TRUE))
range(colMeans(Qcz_values.list$z_lower, na.rm = TRUE))
apply(Qcz_values.list$z_lower, 2, function(x) sd(x, na.rm = TRUE))

plot(colMeans(Qcz_values.list$z_lower, na.rm = TRUE) ~ colMeans(Qcz_values.list$c_lower), xlab = "c", ylab = "z", type = "n")
text(y = colMeans(Qcz_values.list$z_lower, na.rm = TRUE), x = colMeans(Qcz_values.list$c_lower), labels = colnames(Qcz_values.list$z_lower))

node.labels <- c("Ceci-Tory","Rabd.B-Eulo","Rabd.B-Lest","Rabd.B-Meso","Rabd.B-Plat","Rabd.B-Tory","Pont-Eury","Pont_Lath","Rabd.S-Plat","Iteo-Eulo","Iteo-Meso","Iteo-Myma","Iteo-Plat","Iteo-Tory",LETTERS[c(3,1:2,4:20,22:26)]) # for plotting aesthetics

cz_values.df <- data.frame(node.id = c(colnames(Qcz_values.list$c_higher), colnames(Qcz_values.list$c_lower)), 
                           node.set = c(rep("higher", 14), rep("lower",25)),
                           node.labels = node.labels,
                           c = c(colMeans(Qcz_values.list$c_higher), colMeans(Qcz_values.list$c_lower)), 
                           z = c(colMeans(Qcz_values.list$z_higher, na.rm = TRUE), colMeans(Qcz_values.list$z_lower, na.rm = TRUE)))

write.csv(cz_values.df, "~/Documents/Genotype_Networks/data/cz_values.df.csv")

# Null model analysis of role of nodes in network. Use 1000 null models of 'swap.web' (same protocol for null model analysis of modularity)
null_model_czvalues.df <- null_model_analysis_czvalues(web = genotype_gall_parasitoid_network, N_null_webs = 1000, null_model = "swap.web")

null.high.c.95quant <- apply(null_model_czvalues.df$high.null_c.values, 1, function(x) quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))

null.high.z.95quant <- apply(null_model_czvalues.df$high.null_z.values, 1, function(x) quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))

null.low.c.95quant <- apply(null_model_czvalues.df$low.null_c.values, 1, function(x) quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))

null.low.z.95quant <- apply(null_model_czvalues.df$low.null_z.values, 1, function(x) quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))

null.cz.05.95.quantile.df <- data.frame(null.low.c.05.mean = mean(null.low.c.95quant[1, ]),
                                        null.low.c.95.mean = mean(null.low.c.95quant[2, ]),
                                        null.low.z.05.mean = mean(null.low.z.95quant[1, ]),
                                        null.low.z.95.mean = mean(null.low.z.95quant[2, ]),
                                        null.high.c.05.mean = mean(null.high.c.95quant[1, ]),
                                        null.high.c.95.mean = mean(null.high.c.95quant[2, ]),
                                        null.high.z.05.mean = mean(null.high.z.95quant[1, ]),
                                        null.high.z.95.mean = mean(null.high.z.95quant[2, ]),
                                        null.low.c.05.sd = sd(null.low.c.95quant[1, ]),
                                        null.low.c.95.sd = sd(null.low.c.95quant[2, ]),
                                        null.low.z.05.sd = sd(null.low.z.95quant[1, ]),
                                        null.low.z.95.sd = sd(null.low.z.95quant[2, ]),
                                        null.high.c.05.sd = sd(null.high.c.95quant[1, ]),
                                        null.high.c.95.sd = sd(null.high.c.95quant[2, ]),
                                        null.high.z.05.sd = sd(null.high.z.95quant[1, ]),
                                        null.high.z.95.sd = sd(null.high.z.95quant[2, ]))

write.csv(null.cz.05.95.quantile.df, "~/Documents/Genotype_Networks/data/null.cz.05.95.quantile.df.csv")

null.cz.95.summary <- data.frame(c = with(null.cz.05.95.quantile.df, c(null.high.c.95.mean, null.low.c.95.mean)),
                                 z = with(null.cz.05.95.quantile.df, c(null.high.z.95.mean, null.low.z.95.mean)), 
                                 node.set = c("higher","lower"))

##### Plots of c-z values
theme_cz <- theme_facet + theme(axis.title = element_blank(),
                                  axis.text = element_text(size = 18))

cz.low <- ggplot(cz_values.df, aes(x = c, y = z)) + 
  geom_text(data = subset(cz_values.df, node.set == "lower"), aes(label = node.labels)) +
  theme_cz +   
  geom_hline(data = subset(null.cz.95.summary, node.set == "lower"), aes(yintercept = z), linetype = "dashed") + 
  geom_vline(data = subset(null.cz.95.summary, node.set == "lower"), aes(xintercept = c), linetype = "dashed") +
  coord_cartesian(xlim = c(-0.1,1.05)) +
  theme(plot.margin = unit(c(0,1,0.5,0.5),"lines")) +
  annotate("text", x = -0.05, y = 2.1, label = "(B)", size = 8)

cz.high <- ggplot(cz_values.df, aes(x = c, y = z)) + 
  geom_text(data = subset(cz_values.df, node.set == "higher"), aes(label = node.labels), fontface = "italic") +
  theme_cz + 
  geom_hline(data = subset(null.cz.95.summary, node.set == "higher"), aes(yintercept = z), linetype = "dashed") + 
  geom_vline(data = subset(null.cz.95.summary, node.set == "higher"), aes(xintercept = c), linetype = "dashed") +
  coord_cartesian(xlim = c(-0.1,1.05)) +
  scale_y_continuous(breaks = c(-1,0,1)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  annotate("text", x = -0.05, y = 1.4, label = "(A)", size = 8)


gp1<- ggplot_gtable(ggplot_build(cz.high))
gp2<- ggplot_gtable(ggplot_build(cz.low))
maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
gp1$widths[2:3] <- maxWidth
gp2$widths[2:3] <- maxWidth


# consider moving abundance to top, then richness, followed by evenness.
grid.arrange(gp1, gp2, left = textGrob("z", gp = gpar(fontsize = 30, fontface = "italic"), rot = 90, vjust = 0.5), sub = textGrob("c", gp = gpar(fontsize = 30, fontface = "italic"), vjust = 0))

##### Identify frequency of module asssociates for important nodes. 
vLG_Platy_team <- teammates_table(web = genotype_gall_parasitoid_network, focal = "vLG_Platy", QuaBiMo_reps = 100) 
vLG_Mesopol_team <- teammates_table(web = genotype_gall_parasitoid_network, focal = "vLG_Mesopol", QuaBiMo_reps = 100)
vLG_Tory_team <- teammates_table(web = genotype_gall_parasitoid_network, focal = "vLG_Tory", QuaBiMo_reps = 100)
rG_Tory_team <- teammates_table(web = genotype_gall_parasitoid_network, focal = "rG_Tory", QuaBiMo_reps = 2)
I_team <- teammates_table(web = genotype_gall_parasitoid_network, focal = "I", node.location = "row", QuaBiMo_reps = 2) 
