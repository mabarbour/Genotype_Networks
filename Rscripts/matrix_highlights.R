# matrix highlighter function for plotweb2_highlights

matrix_highlighter <- function(web, high.names = "", low.names = "") {
  
  web <- as.matrix(web)
  
  grab.colnames <- which(colnames(web) %in% high.names)
  grab.rownames <- which(rownames(web) %in% low.names)
  
  web.sub <- web
  web.sub[ ,-grab.colnames] <- 0
  web.sub[-grab.rownames, ] <- 0
  
  out <- web.sub
}


matrix_highlighter_2 <- function(web, interaction.names = "", replace.values = NULL) {
  require(dplyr)
  
  if(!is.null(replace.values)){
    replace.df <- data.frame(interaction.names = interaction.names, replace.values = replace.values)
  }
  
  web <- data.frame(web)
  web$low.nodes <- factor(rownames(web))
  
  web.melt <- melt(web)
  web.melt <- mutate(web.melt, interactions = paste(low.nodes, variable, sep = "_"))
  
  web.melt.sub <- web.melt
  web.melt.sub[-which(web.melt.sub$interactions %in% interaction.names), 3] <- 0
  
  if(!is.null(replace.values)){
    interaction.order <- web.melt.sub$interactions[which(web.melt.sub$interactions %in% interaction.names)]
    replace.df <- mutate(replace.df, interaction.order = factor(interaction.order, ordered = TRUE))
    test <- match(replace.df$interaction.names, replace.df$interaction.order)
    #replace.df.sort <- sort_df(replace.df, "interaction.order")
    test2 <- web.melt.sub[which(web.melt.sub$interactions %in% interaction.names), ]
    test3 <- replace.df$replace.values[test]
    
    web.melt.sub <- web.melt.sub %>% 
      left_join(replace.df %>% select(interactions = interaction.names, replace.values),
                missing = 0) %>%
      mutate(value = ifelse(is.na(replace.values),yes = 0,no = replace.values))
  }
    #web.melt.sub[which(web.melt.sub$interactions %in% interaction.names), 3] <- replace.df$replace.values[test]
  #}

  web.melt.sub.matrix <- dcast(web.melt.sub, low.nodes ~ variable, value.var = "value")
  rownames(web.melt.sub.matrix) <- web.melt.sub.matrix$low.nodes
  web.melt.sub.matrix <- web.melt.sub.matrix[ ,-1]
  web.melt.sub.matrix <- as.matrix(web.melt.sub.matrix)
  
  out <- web.melt.sub.matrix
  #web.melt.expand <- expand.grid(web.melt)
}

#debug(matrix_highlighter_2)


