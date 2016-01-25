# Bipartite network with tapered, intensity curved edges

# load libraries
library(sna)
library(ggplot2)
library(Hmisc)
library(reshape2)
require(plyr)
require(dplyr)
require(vegan)



# load data
genotype_gall_parasitoid_network <- read.csv('~/Documents/Genotype_Networks/data/genotype_gall_parasitoid_network_data.csv')
rownames(genotype_gall_parasitoid_network) <- genotype_gall_parasitoid_network$X 
genotype_gall_parasitoid_network <- select(genotype_gall_parasitoid_network, -X)

module.info <- read.csv('~/Documents/Genotype_Networks/data/module_info_genotype_gall_parasitoid_network.csv')
module.info <- arrange(module.info, Trophic, Module.ID, Vertex) # reorder for easy plotting
genotypes <- filter(module.info, Trophic == "Genotype")
gall.parasitoid <- filter(module.info, Trophic == "Gall-Parasitoid")

web <- genotype_gall_parasitoid_network
web <- web[as.character(genotypes$Vertex), # reorder within web
           as.character(gall.parasitoid$Vertex)] # reorder within web

web <- as.matrix(web)
attr(web, "modules") <- list(genotypes$Module.ID,
                             gall.parasitoid$Module.ID)

# Requires an ordered matrix with a "modules" attributes giving the module ID for each node of the bipartite network

bipartite_info <- function(web){

  # layout for bipartite network. Evenly spaces out trophic level with the fewest number of nodes
  layout.bipartite <- matrix(nrow = sum(dim(web)), ncol = 2)
  layout.bipartite[ ,1] <- c(1:max(dim(web)), 
                             seq(1, max(dim(web)), length.out = min(dim(web)))) 
  layout.bipartite[ ,2] <- c(rep(which.max(dim(web)), max(dim(web))), rep(which.min(dim(web)), min(dim(web))))
  
  # name a replicate of the web "adjacencyMatrix" with unique, numbered values for each node. Necessary for accurate plotting.
  adjacencyMatrix <- web
  rownames(adjacencyMatrix) <- 1:dim(web)[1] 
  colnames(adjacencyMatrix) <- (dim(web)[1]+1):sum(dim(web)) 

  # get layout coordinates for graph
  layoutCoordinates <- gplot(adjacencyMatrix, 
                             gmode = "twomode", 
                             coord = layout.bipartite,
                             jitter = FALSE, 
                             label = c(rownames(web), colnames(web))) # double check that network is plotting correctly.

  layout_df <- data.frame(layoutCoordinates, 
                          vertex.names = c(dimnames(web)[[which.max(dim(web))]],
                                           dimnames(web)[[which.min(dim(web))]]),
                          trophic = factor(c(rep(which.max(dim(web)), max(dim(web))), 
                                             rep(which.min(dim(web)), min(dim(web))))),
                          module = factor(c(attr(web, which = "modules")[[which.max(dim(web))]],
                                            attr(web, which = "modules")[[which.min(dim(web))]]))
  )
     
  adjacencyList <- melt(adjacencyMatrix)  # Convert to list of ties only
  adjacencyList <- adjacencyList[adjacencyList$value > 0, ]

  # Function to generate paths between each connected node. Adapted from the following     gist -> http://is-r.tumblr.com/post/38459242505/beautiful-network-diagrams-with-ggplot2
  edgeMaker <- function(whichRow, len = 100, curved = TRUE){
    fromC <- layoutCoordinates[adjacencyList[whichRow, 1], ]  # Origin
    toC <- layoutCoordinates[adjacencyList[whichRow, 2], ]  # Terminus
  
    # Add curve:
    graphCenter <- colMeans(layoutCoordinates)  # Center of the overall graph
    bezierMid <- c(fromC[1], toC[2])  # A midpoint, for bended edges
    distance1 <- sum((graphCenter - bezierMid)^2)
    if(distance1 < sum((graphCenter - c(toC[1], fromC[2]))^2)){
      bezierMid <- c(toC[1], fromC[2])
     }  # To select the best Bezier midpoint
    bezierMid <- (fromC + toC + bezierMid) / 3  # Moderate the Bezier midpoint
    if(curved == FALSE){bezierMid <- (fromC + toC) / 2}  # Remove the curve
  
    edge <- data.frame(bezier(c(fromC[1], bezierMid[1], toC[1]),  # Generate
                            c(fromC[2], bezierMid[2], toC[2]),  # X & y
                            evaluation = len))  # Bezier path coordinates
    edge$Sequence <- 1:len  # For size and colour weighting in plot
    edge$Group <- paste(adjacencyList[whichRow, 1:2], collapse = ">")
    edge$Weight <- adjacencyList[whichRow, 3]
    return(edge)
  }
  # Generate a (curved) edge path for each pair of connected nodes
  allEdges <- lapply(1:nrow(adjacencyList), edgeMaker, len = 500, curved = TRUE)
  allEdges <- do.call(rbind, allEdges)  # a fine-grained path ^, with bend ^

  out <- list(allEdges, layout_df)
}

test <- bipartite_info(web)

bipartite_network_plot <- function(bipartite.info, size = 15, range = c(3/10, 3), shape_values = c(21,24)){
  # Empty ggplot2 theme
  new_theme_empty <- theme_bw()
  new_theme_empty$line <- element_blank()
  new_theme_empty$rect <- element_blank()
  new_theme_empty$strip.text <- element_blank()
  new_theme_empty$axis.text <- element_blank()
  new_theme_empty$plot.title <- element_blank()
  new_theme_empty$axis.title <- element_blank()
  new_theme_empty$plot.margin <- structure(c(0, 0, -1, -1), unit = "lines",
                                           valid.unit = 3L, class = "unit")
  
  # ggplot
  ggplot(bipartite.info[[1]]) + 
    geom_path(aes(x = x, y = y, group = Group, # edges with gradient
                  color = -Sequence*Weight, # add taper
                  size = Sequence*Weight)) + # add taper
    geom_point(data = bipartite.info[[2]],
               aes(x = x, y = y, fill = module, shape = trophic), size = size) +
    scale_fill_discrete(guide = "none") + 
    scale_shape_manual(values = shape_values, guide = "none") +
    scale_colour_gradient(low = gray(0), high = gray(9/10), guide = "none") +
    scale_size(range = range, guide = "none") + 
    new_theme_empty
}

test.plot <- bipartite_network_plot(bipartite.info = test)
test.plot + scale_fill_manual(values = c("black","red","blue","grey","orange"), guide = "none")
cb_palette <- c("#56B4E9", "#999999", "#E69F00", "#009E73", "#CC79A7")


n
zp1 <- zp1 + geom_point(data = layout_df,  # Add nodes
                        aes(x = x, y = y, fill = Module.ID, shape = Trophic), 
                        color = "black", size = 15) +
  scale_fill_manual(values = cb_palette, guide = "none") + 
  scale_shape_manual(values = c(23, 21), guide = "none")
#zp1 <- zp1 + geom_text(data = layout_df,
#  aes(x = x, y = y, label = vertex.names.plot), size = 9)
zp1 <- zp1 + scale_colour_gradient(low = gray(0), high = gray(9/10), guide = "none")
zp1 <- zp1 + scale_size(range = c(5/10, 3), guide = "none")  # Customize taper c(1/10,1)
zp1 <- zp1 + new_theme_empty  # Clean up plot
print(zp1) # save as png, maintain aspect ratio and have width of at least 1000
# Looks better when saved as a PNG:
#ggsave("genotype_gall-parasitoid_network.png", zp1, h =2, w = 3, units = "in", type = "cairo-png", path = "~/Documents/Genotype_Networks/")


test <- bipartite_plot_info(web)
test.plot <- ggnet_bipartite(test)
test.plot + geom_point(data = test[[2]], aes(x = x, y = y, fill = trophic, shape = trophic), size = 12, guide = "none") + scale_shape_manual(values = c(21,23), guide = "none")
