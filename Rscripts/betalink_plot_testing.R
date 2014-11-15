source('~/Documents/ggnet/bipartite_plot_info.R')
source('~/Documents/ggnet/ggnet_bipartite.R')

betalink_bipartite_plot_info <- function(m1, m2, order.type = "cca") {
  m1.info <- bipartite_plot_info(m1, order.type = order.type)
  m2.info <- bipartite_plot_info(m2, order.type = order.type)
  
  m1.m2.bind <- rbind(m1.info[[2]], m2.info[[2]])
  m1.m2.bind$webID <- factor(c(rep(1, sum(dim(m1))), 
                               rep(2, sum(dim(m2)))))
  
  up.1 = rownames(m1)
  up.2 = rownames(m2)
  lo.1 = colnames(m1)
  lo.2 = colnames(m2)
  up.shared = up.1[up.1%in%up.2]
  lo.shared = lo.1[lo.1%in%lo.2]
  up.un.1 = up.1[!(up.1%in%up.shared)]
  up.un.2 = up.2[!(up.2%in%up.shared)]
  lo.un.1 = lo.1[!(lo.1%in%lo.shared)]
  lo.un.2 = lo.2[!(lo.2%in%lo.shared)]
  
  shared = c(up.shared, lo.shared)
  unshared = c(up.un.1, up.un.2, lo.un.1, lo.un.2)
  
  m1.m2.bind$species.turnover.colors <- rep(NA, dim(m1.m2.bind)[1])
  m1.m2.bind$species.turnover.colors[which(m1.m2.bind$vertex.names %in% shared)] <- "shared"
  m1.m2.bind$species.turnover.colors[which(m1.m2.bind$vertex.names %in% unshared)] <- "unshared"
  out <- list(m1.info, m2.info, m1.m2.bind)
}

comp1 <- betalink_bipartite_plot_info(m1 = genotype_net_matrices[[1]], m2 = genotype_net_matrices[[3]], order.type = "normal")

m1.plot <- ggnet_bipartite(comp1[[1]]) + 
  geom_point(data = filter(comp1[[3]], webID == 1), 
             aes(x = x, y = y, fill = species.turnover.colors), size = 12, shape = 21) + 
  scale_fill_manual(values = c("gray", "steelblue"), guide = "none") + 
  geom_text(data = filter(comp1[[3]], webID == 1), 
            aes(x = x, y = y, label = vertex.names))

m2.plot <- ggnet_bipartite(comp1[[2]]) +
  geom_point(data = filter(comp1[[3]], webID == 2), 
             aes(x = x, y = y, fill = species.turnover.colors), size = 12, shape = 21) + 
  scale_fill_manual(values = c("gray", "steelblue"), guide = "none") + 
  geom_text(data = filter(comp1[[3]], webID == 2), 
            aes(x = x, y = y, label = vertex.names))

grid.arrange(m1.plot, m2.plot, ncol = 2)

# may need to consider melting dataframes for comparison...

test1 <- bipartite_plot_info(genotype_net_matrices[[1]], order.type = "cca")
#test1[[2]]$trophic <- factor(c(1,1,2,2,2))

test2 <- bipartite_plot_info(genotype_net_matrices[[3]], order.type = "cca")
test2.plot <- ggnet_bipartite(test2)

test1.plot <- ggnet_bipartite(test1, range = c(min(genotype_net_matrices[[1]]), max(genotype_net_matrices[[1]])))
test1.plot + 
  geom_point(data = test1[[2]], aes(x = x, y = y, shape = vertex.names, fill = trophic), size = 15) +
  scale_shape_manual(values = 21:25, guide = "none") + 
  scale_fill_manual(values = c("steelblue", "red"), guide = "none")

library(gridExtra)
grid.arrange(test1.plot, test2.plot, ncol = 2)



betalink.plot(genotype_net_matrices[[1]], genotype_net_matrices[[3]])
source('~/Documents/betalink/R/metaweb.r')

test.bind <- rbind(test1[[2]], test2[[2]])
test.bind$webID <- factor(c(rep(1, dim(test1[[2]])[1]), 
                            rep(2, dim(test2[[2]])[1])))

m1 <- genotype_net_matrices[[1]]
m2 <- genotype_net_matrices[[2]]

which(test.bind$vertex.names == up.shared | test.bind$vertex.names == lo.shared)

unshared <- c(up.un.1, up.un.2, lo.un.1, lo.un.2)
which(test.bind$vertex.names %in% unshared)

betalink.plot.test = function(m1,m2,by.unique=TRUE){
  up.1 = rownames(m1)
  up.2 = rownames(m2)
  lo.1 = colnames(m1)
  lo.2 = colnames(m2)
  up.shared = up.1[up.1%in%up.2]
  lo.shared = lo.1[lo.1%in%lo.2]
  up.un.1 = up.1[!(up.1%in%up.shared)]
  up.un.2 = up.2[!(up.2%in%up.shared)]
  lo.un.1 = lo.1[!(lo.1%in%lo.shared)]
  lo.un.2 = lo.2[!(lo.2%in%lo.shared)]
  m = metaweb(list(m1,m2))$web
  if(by.unique){
    oCol = rep(0,ncol(m))
    oCol[colnames(m)%in%lo.un.1] = c(1:length(lo.un.1))
    oCol[colnames(m)%in%lo.shared] = max(oCol)+c(1:length(lo.shared))
    oCol[colnames(m)%in%lo.un.2] = max(oCol)+c(1:length(lo.un.2))
    oRow = rep(0,nrow(m))
    oRow[rownames(m)%in%up.un.1] = c(1:length(up.un.1))
    oRow[rownames(m)%in%up.shared] = max(oRow)+c(1:length(up.shared))
    oRow[rownames(m)%in%up.un.2] = max(oRow)+c(1:length(up.un.2))
  } else {
    oCol = ncol(m)-rank(colSums(m),ties.method='r')+1
    oRow = nrow(m)-rank(rowSums(m),ties.method='r')+1
  }
  out = list(oCol, oRow)
}
t <- betalink.plot.test(m1, m2)
oRow = t[[2]]
oCol = t[[1]]
ms <- vector()

for(ro in 1:nrow(m)){
  current.up = rownames(m)[ro]
  for(co in 1:ncol(m)){
    current.low = colnames(m)[co]
    c.BG = 'lightgrey'
    if((current.up%in%up.un.1) | (current.low%in%lo.un.1)){c.BG = 'palegreen'}
    if((current.up%in%up.un.2) | (current.low%in%lo.un.2)){c.BG = 'skyblue'}
    if((current.up%in%up.shared) & (current.low%in%lo.shared)){
      if(m1[current.up,current.low]!=m2[current.up,current.low]){
        c.BG = 'orange'
      }
    }
    if(m[current.up,current.low] == 0){c.BG = 'white'}
    rect(oRow[ro]-0.9,oCol[co]-0.9,oRow[ro]-0.1,oCol[co]-0.1,col=c.BG,border=NA)
  }
}