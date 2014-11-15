
plotweb2_highlights <- function (web, 
                                 web2, 
                                 highlight.web = NULL,
                                 highlight.web.high = NULL,
                                 method = "cca", 
                                 empty = FALSE, # when FALSE, method changes to "normal"
                                 labsize = 1, 
                                 ybig = 1, 
                                 y_width = 0.1,
                                 y_width_scale = 0,
                                 spacing = 0.05, 
                                 arrow = "no", 
                                 col.interaction = "grey80", 
                                 bord.col.interaction = "black",
                                 col.pred = "grey10", 
                                 col.prey = "grey10",
                                 lab.space = 1, 
                                 lablength = NULL, # buggy if not set to null
                                 sequence = NULL, 
                                 low.abun = NULL, 
                                 low.abun.col = "green", 
                                 high.abun = NULL, 
                                 high.abun.col = "red", 
                                 method2 = "cca", 
                                 empty2 = TRUE, 
                                 spacing2 = 0.05, 
                                 arrow2 = "no", 
                                 col.interaction2 = "grey80", 
                                 bord.col.interaction2 = "black",
                                 col.pred2 = "grey10", 
                                 col.prey2 = NA, 
                                 lablength2 = NULL, # buggy if not set to null
                                 sequence.pred2 = NULL, 
                                 low.abun2 = NULL, 
                                 low.abun.col2 = "white", 
                                 bord.low.abun.col2 = "black",
                                 high.abun2 = NULL,
                                 high.abun.col2 = "white",
                                 col.interaction.highlight = "red",
                                 bord.col.highlight = "black",
                                 bord.col = "black",
                                 bord.col.interaction.highlight = "black", 
                                 bord.size.interaction.highlight = 1,
                                 node.highlight = "",
                                 bord.size = 1,
                                 bord.size.highlight = 5,
                                 font.type = 3,
                                 text = FALSE) 
{
  require(vegan)
  if (empty) 
    web <- empty(web)
  else method <- "normal"
  web <- as.matrix(web)
  meths <- c("normal", "cca")
  meths.match <- pmatch(method, meths)
  if (is.na(meths.match)) 
    stop("Choose plot-method: normal/cca.\n")
  if (meths.match == 2) {
    ca <- cca(web)
    web <- web[order(summary(ca)$sites[, 1], decreasing = TRUE), 
               order(summary(ca)$species[, 1], decreasing = TRUE)]
  }
  if (!is.null(sequence)) {
    cs <- sequence$seq.pred %in% colnames(web)
    rs <- sequence$seq.prey %in% rownames(web)
    web <- web[sequence$seq.prey[rs], sequence$seq.pred[cs]]
  }
  websum <- sum(web)
  difff <- diffh <- 0
  if (!is.null(low.abun)) {
    lowfreq = rowSums(web)
    dummy <- lowfreq
    for (i in 1:length(low.abun)) {
      ind <- which(names(low.abun)[i] == names(dummy))
      lowfreq[ind] <- lowfreq[ind] + low.abun[i]
    }
    difff = (lowfreq - rowSums(web))/websum
  }
  if (!is.null(high.abun)) {
    highfreq = colSums(web)
    dummy <- highfreq
    for (i in 1:length(high.abun)) {
      ind <- which(names(high.abun)[i] == names(dummy))
      highfreq[ind] <- highfreq[ind] + high.abun[i]
    }
    diffh = (highfreq - colSums(web))/websum
  }
  if (is.null(high.abun)) 
    pred_prop <- colSums(web)/websum
  else pred_prop <- highfreq/websum
  if (is.null(low.abun)) 
    prey_prop <- rowSums(web)/websum
  else prey_prop <- lowfreq/websum
  n.pred <- length(pred_prop)
  n.prey <- length(prey_prop)
  pred_x <- 0
  pred_xold <- -1
  pred_versatz <- 0
  pred_y <- 1.5
  prey_x <- 0
  prey_xold <- -1
  prey_versatz <- 0
  prey_y <- 0.5
  if (length(colnames(web)) == 0) 
    colnames(web) <- colnames(web, do.NULL = FALSE)
  if (length(rownames(web)) == 0) 
    rownames(web) <- rownames(web, do.NULL = FALSE)
  if (!is.null(lablength)) 
    colnames(web) <- substr(colnames(web), 1, lablength)
  if (!is.null(lablength)) 
    rownames(web) <- substr(rownames(web), 1, lablength)
  par(mai = c(0.2, 0.2, 0.2, 0.2))
  #browser()
  pred_spacing = (n.prey - 1)/(n.pred - 1)
  prey_spacing = (n.pred - 1)/(n.prey - 1)
  pred_spacing <- pred_spacing * spacing
  prey_spacing <- prey_spacing * spacing
  #if (n.pred > n.prey & n.prey == 1){ # added this if statment to permit plots with a single prey species
   # prey_spacing <- 0 
  #}
  if (n.pred > n.prey) # added this to be distinct from above statement
    prey_spacing <- pred_spacing * (n.pred - 1)/(n.prey - 
                                                   1)
  else pred_spacing <- prey_spacing * (n.prey - 1)/(n.pred - 
                                                      1)
  if (n.pred == 1)
    pred_spacing <- 1 #* spacing
  if (n.prey == 1)
    prey_spacing <- 1 #* spacing
    pred_spacing <- spacing # adds some space in between predator species so they can be distinguished.
  if (!is.null(low.abun)) 
    pred_spacing <- pred_spacing + sum(difff)/n.pred
  if (!is.null(high.abun)) 
    prey_spacing <- prey_spacing + sum(diffh)/n.prey
  #browser()
  wleft = 0
  wright = (max(n.pred, n.prey)) * min(prey_spacing, pred_spacing) + 
    1 + max(sum(diffh), sum(difff))
  wup <- 2.6 + y_width + lab.space * 0.05
  wdown <- 0.4 - y_width - lab.space * 0.05
  plot(0, type = "n", xlim = range(wleft, wright), ylim = range(wdown/ybig, 
                                                                wup * ybig), axes = FALSE, xlab = "", ylab = "")

  pred_x = 0
  hoffset <- 0
  links <- 0
  rechts <- 0
  hoehe <- strheight(colnames(web)[1], cex = 0.6)
  for (i in 1:n.pred) {
    # attempt to add prey border highlight
    if(names(pred_prop[i]) %in% node.highlight)
      bord.col.mid = c(bord.col.highlight, bord.size.highlight)
    else {
      bord.col.mid = c(bord.col, bord.size)
    }
    rect(pred_x, pred_y - y_width, pred_x + pred_prop[i], 
         pred_y + y_width, col = col.pred, border = bord.col.mid[1], lwd = bord.col.mid[2])
    if (!is.null(high.abun)) {
      rect(pred_x + pred_prop[i] - diffh[i], pred_y - y_width, 
           pred_x + pred_prop[i], pred_y + y_width, col = high.abun.col, border = bord.col.mid[1], lwd = bord.col.mid[2])
    }
    breite <- strwidth(colnames(web)[i], cex = 0.6 * labsize)
    links <- pred_x + pred_prop[i]/2 - breite/2
    if (links < rechts && i > 1) 
      hoffset = hoffset + hoehe
    else {
      rechts <- pred_x + pred_prop[i]/2 + breite/2
      hoffset <- 0
    }
    if(text == TRUE){
      text(pred_x + pred_prop[i]/2, pred_y + y_width + hoehe + 
             hoffset, colnames(web)[i], cex = 0.6 * labsize, offset = 0, font = font.type)
    }
      pred_x <- pred_x + pred_prop[i] + pred_spacing
  }

  prey_x <- 0
  links <- 0
  rechts <- 0
  hoehe <- strheight(rownames(web)[1], cex = 0.6)
  hoffset <- hoehe
  for (i in 1:n.prey) {
    # attempt to add prey border highlight
    if(names(prey_prop[i]) %in% node.highlight)
      bord.col.prey = c(bord.col.highlight, bord.size.highlight)
    else {
      bord.col.prey = c(bord.col, bord.size)
    }
    rect(prey_x, prey_y - y_width, prey_x + prey_prop[i], 
         prey_y + y_width, col = col.prey, border = bord.col.prey[1], lwd = bord.col.prey[2])
    if (!is.null(low.abun)) {
      rect(prey_x + prey_prop[i] - difff[i], prey_y - y_width, 
           prey_x + prey_prop[i], prey_y + y_width, col = low.abun.col)
    }
    breite <- strwidth(rownames(web)[i], cex = 0.6 * labsize)
    links <- prey_x + prey_prop[i]/2 - breite/2
    if (links < rechts && i > 1) 
      hoffset = hoffset + hoehe
    else {
      rechts <- prey_x + prey_prop[i]/2 + breite/2
      hoffset <- hoehe
    }
    if(text == TRUE){
    text(prey_x + prey_prop[i]/2, prey_y - y_width - hoffset, 
         rownames(web)[i], cex = 0.6 * labsize, offset = 0, font = font.type)
    }
    prey_x <- prey_x + prey_prop[i] + prey_spacing
  }

  pred_x <- 0
  zwischenweb <- web
  XYcoords1 <- matrix(ncol = 2, nrow = length(zwischenweb))
  for (i in 1:length(zwischenweb)) {
    XYcoords1[i, ] <- which(zwischenweb == max(zwischenweb), 
                           arr.ind = TRUE)[1, ]
    zwischenweb[XYcoords1[i, 1], XYcoords1[i, 2]] <- -1
  }
  y1 <- pred_y - y_width - y_width_scale
  y2 <- y1
  y3 <- prey_y + y_width + y_width_scale
  y4 <- y3
  for (p in 1:sum(web > 0)) {
    i <- XYcoords1[p, 1]
    j <- XYcoords1[p, 2]
    if (j == 1 & i == 1) 
      x1 <- 0
    else x1 <- (j - 1) * pred_spacing + cumsum(web)[(j - 
                                                       1) * nrow(web) + (i - 1)]/websum
    if (!is.null(high.abun) && j > 1) 
      x1 <- x1 + cumsum(diffh)[j - 1]
    x2 <- x1 + web[i, j]/websum
    if (arrow == "up" || arrow == "both") {
      x2 <- (x1 + x2)/2
      x1 <- x2
    }
    tweb <- t(web)
    if (j == 1 & i == 1) 
      x3 <- 0
    else x3 <- (i - 1) * prey_spacing + cumsum(tweb)[(i - 
                                                        1) * nrow(tweb) + (j - 1)]/websum
    if (!is.null(low.abun) && i > 1) 
      x3 <- x3 + cumsum(difff)[i - 1]
    x4 <- x3 + tweb[j, i]/websum
    if (arrow == "down" || arrow == "both") {
      x4 <- (x3 + x4)/2
      x3 <- x4
    }
    polygon(c(x1, x2, x4, x3), c(y1, y2, y4, y3), col = col.interaction, border = bord.col.interaction)
  }
  # attempt to highlight specific interactions in lower plots
  if(!is.null(highlight.web)){
  pred_x <- 0
  zwischenweb <- highlight.web
  XYcoords1 <- matrix(ncol = 2, nrow = length(zwischenweb))
  for (i in 1:length(zwischenweb)) {
    XYcoords1[i, ] <- which(zwischenweb == max(zwischenweb), 
                            arr.ind = TRUE)[1, ]
    zwischenweb[XYcoords1[i, 1], XYcoords1[i, 2]] <- -1
  }
  y1 <- pred_y - y_width - y_width_scale
  y2 <- y1
  y3 <- prey_y + y_width + y_width_scale
  y4 <- y3
  for (p in 1:sum(highlight.web > 0)) { # web
    i <- XYcoords1[p, 1]
    j <- XYcoords1[p, 2]
    if (j == 1 & i == 1) 
      x1 <- 0
    else x1 <- (j - 1) * pred_spacing + cumsum(web)[(j - 
                                                       1) * nrow(web) + (i - 1)]/websum
    if (!is.null(high.abun) && j > 1) 
      x1 <- x1 + cumsum(diffh)[j - 1]
    x2 <- x1 + web[i, j]/websum
    if (arrow == "up" || arrow == "both") {
      x2 <- (x1 + x2)/2
      x1 <- x2
    }
    tweb <- t(web)
    tweb2 <- t(highlight.web)
    if (j == 1 & i == 1) 
      x3 <- 0
    else x3 <- (i - 1) * prey_spacing + cumsum(tweb)[(i - 
                                                        1) * nrow(tweb) + (j - 1)]/websum
    if (!is.null(low.abun) && i > 1) 
      x3 <- x3 + cumsum(difff)[i - 1]
    x4 <- x3 + tweb2[j, i]/websum # if highlight interaction value is smaller, it will be plotted proportionately to the rest of the web.
    if (arrow == "down" || arrow == "both") {
      x4 <- (x3 + x4)/2
      x3 <- x4
    }
    polygon(c(x1, x2, x4, x3), c(y1, y2, y4, y3), col = col.interaction.highlight, border = bord.col.interaction.highlight, lwd = bord.size.interaction.highlight)
  }
  }
  # replot base for clear borders
  prey_x <- 0
  links <- 0
  rechts <- 0
  hoehe <- strheight(rownames(web)[1], cex = 0.6)
  hoffset <- hoehe
  for (i in 1:n.prey) {
    # attempt to add prey border highlight
    if(names(prey_prop[i]) %in% node.highlight)
      bord.col.prey = c(bord.col.highlight, bord.size.highlight)
    else {
      bord.col.prey = c(bord.col, bord.size)
    }
    rect(prey_x, prey_y - y_width, prey_x + prey_prop[i], 
         prey_y + y_width, col = col.prey, border = bord.col.prey[1], lwd = bord.col.prey[2])
    if (!is.null(low.abun)) {
      rect(prey_x + prey_prop[i] - difff[i], prey_y - y_width, 
           prey_x + prey_prop[i], prey_y + y_width, col = low.abun.col)
    }
    breite <- strwidth(rownames(web)[i], cex = 0.6 * labsize)
    links <- prey_x + prey_prop[i]/2 - breite/2
    if (links < rechts && i > 1) 
      hoffset = hoffset + hoehe
    else {
      rechts <- prey_x + prey_prop[i]/2 + breite/2
      hoffset <- hoehe
    }
    if(text == TRUE){
    text(prey_x + prey_prop[i]/2, prey_y - y_width - hoffset, 
         rownames(web)[i], cex = 0.6 * labsize, offset = 0, font = font.type)
    }
    prey_x <- prey_x + prey_prop[i] + prey_spacing
  }
  #browser() # setting lablength different from null screws up the function.
  #
  web2 <- as.matrix(web2)
  for (i in 1:dim(web)[2]) {
    dn <- dimnames(web)[[2]][i]
    if (is.na(match(dn, dimnames(web2)[[1]]))) {
      dummy <- matrix(rep(0, dim(web2)[2]), nrow = 1)
      rownames(dummy) <- dn
      web2 <- rbind(web2, dummy)
    }
  }
  web2 <- web2[order(dimnames(web2)[[1]]), ]
  web2 <- web2[rank(dimnames(web)[[2]]), ]
  difff <- diffh <- 0
  dummy <- colSums(web)
  lowfreq = rowSums(web2)
  for (i in 1:length(dummy)) {
    dummy[i] <- dummy[i] - lowfreq[which(names(lowfreq) == 
                                           names(dummy[i]))]
  }
  low.abun2 <- dummy
  lowfreq = lowfreq + low.abun2
  difff = low.abun2/websum
  if (!is.null(high.abun2)) {
    highfreq = colSums(web2)
    dummy <- highfreq
    for (i in 1:length(high.abun2)) {
      ind <- which(names(high.abun2)[i] == names(dummy))
      highfreq[ind] <- highfreq[ind] + high.abun2[i]
    }
    diffh = (highfreq - colSums(web2))/websum
  }
  if (is.null(high.abun2)) 
    pred_prop <- colSums(web2)/websum
  else pred_prop <- highfreq/websum
  if (is.null(low.abun2)) 
    prey_prop <- rowSums(web2)/websum
  else prey_prop <- lowfreq/websum
  n.pred <- length(pred_prop)
  n.prey <- length(prey_prop)
  pred_x <- 0
  pred_xold <- -1
  pred_versatz <- 0
  pred_y <- 2.5
  prey_x <- 0
  prey_xold <- -1
  prey_versatz <- 0
  prey_y <- 1.5
  if (length(colnames(web2)) == 0) 
    colnames(web2) <- colnames(web2, do.NULL = FALSE)
  if (length(rownames(web2)) == 0) 
    rownames(web2) <- rownames(web2, do.NULL = FALSE)
  if (!is.null(lablength2)) 
    colnames(web2) <- substr(colnames(web2), 1, lablength2)
  if (!is.null(lablength2)) 
    rownames(web2) <- substr(rownames(web2), 1, lablength2)
  prey_spacing = pred_spacing
  if(pred_spacing == 0) # trying ifelse statement to avoid overplotting predators outside the plot space.
    pred_spacing = 0
  else{
    pred_spacing = (n.prey - 1)/(n.pred - 1)
  }
  pred_spacing <- pred_spacing * spacing2
  if (n.pred < n.prey) 
    pred_spacing <- prey_spacing * (n.prey - 1)/(n.pred - 
                                                   1)
  if (!is.null(low.abun2)) 
    pred_spacing <- pred_spacing + sum(difff)/n.pred
  pred_x = 0
  hoffset <- 0
  links <- 0
  rechts <- 0
  hoehe <- strheight(colnames(web2)[1], cex = 0.6)
  for (i in 1:n.pred) {
    # attempt to add predator border highlight
    if(names(pred_prop[i]) %in% node.highlight)
      bord.col.pred = c(bord.col.highlight, bord.size.highlight)
    else {
      bord.col.pred = c(bord.col, bord.size)
    }
    rect(pred_x, pred_y - y_width, pred_x + pred_prop[i], 
         pred_y + y_width, col = col.pred2, border = bord.col.pred[1], lwd = bord.col.pred[2])
    if (!is.null(high.abun2)) {
      rect(pred_x + pred_prop[i] - diffh[i], pred_y - y_width, 
           pred_x + pred_prop[i], pred_y + y_width, col = high.abun.col2)
    }
    breite <- strwidth(colnames(web2)[i], cex = 0.6 * labsize)
    links <- pred_x + pred_prop[i]/2 - breite/2
    if (links < rechts && i > 1) 
      hoffset = hoffset + hoehe
    else {
      rechts <- pred_x + pred_prop[i]/2 + breite/2
      hoffset <- 0
    }
    if(text == TRUE){
    text(pred_x + pred_prop[i]/2, pred_y + y_width + hoehe + 
           hoffset, colnames(web2)[i], cex = 0.6 * labsize, 
         offset = 0, font = font.type)
    }
    pred_x <- pred_x + pred_prop[i] + pred_spacing
  }
  prey_x <- 0
  links <- 0
  rechts <- 0
  hoehe <- strheight(rownames(web2)[1], cex = 0.6)
  hoffset <- hoehe
  for (i in 1:n.prey) {
    # attempt to add prey border highlight
    if(names(prey_prop[i]) %in% node.highlight)
      bord.col.prey = c(bord.col.highlight, bord.size.highlight)
    else {
      bord.col.prey = c(bord.col, bord.size)
    }
    if (!is.null(low.abun2)) {
      rect(prey_x + prey_prop[i] - difff[i], prey_y, prey_x + 
             prey_prop[i], prey_y + y_width, col = low.abun.col2, border = bord.low.abun.col2 )
    } 
    
    rect(prey_x, prey_y - y_width, prey_x + prey_prop[i], 
         prey_y + y_width, col = col.prey2, border = bord.col.prey[1], lwd = bord.col.prey[2])

    breite <- strwidth(rownames(web)[i], cex = 0.6 * labsize)
    links <- prey_x + prey_prop[i]/2 - breite/2
    if (links < rechts && i > 1) 
      hoffset = hoffset + hoehe
    else {
      rechts <- prey_x + prey_prop[i]/2 + breite/2
      hoffset <- hoehe
    }
    if(text == TRUE){
    text(prey_x + prey_prop[i]/2, prey_y - y_width - hoffset, 
         rownames(web2)[i], cex = 0.6 * labsize, offset = 0, font = font.type)
    }
    prey_x <- prey_x + prey_prop[i] + prey_spacing
  }
  pred_x <- 0
  zwischenweb <- web2
  XYcoords2 <- matrix(ncol = 2, nrow = length(zwischenweb))
  for (i in 1:length(zwischenweb)) {
    XYcoords2[i, ] <- which(zwischenweb == max(zwischenweb), 
                           arr.ind = TRUE)[1, ]
    zwischenweb[XYcoords2[i, 1], XYcoords2[i, 2]] <- -1
  }
  y1 <- pred_y - y_width - y_width_scale
  y2 <- y1
  y3 <- prey_y + y_width + y_width_scale
  y4 <- y3
  if (sum(web2 > 0)) {
    for (p in 1:sum(web2 > 0)) {
      i <- XYcoords2[p, 1]
      j <- XYcoords2[p, 2]
      if (j == 1 & i == 1) 
        x1 <- 0
      else x1 <- (j - 1) * pred_spacing + cumsum(web2)[(j - 
                                                          1) * nrow(web2) + (i - 1)]/websum
      if (!is.null(high.abun2) && j > 1) 
        x1 <- x1 + cumsum(diffh)[j - 1]
      x2 <- x1 + web2[i, j]/websum
      if (arrow == "up" || arrow == "both") {
        x2 <- (x1 + x2)/2
        x1 <- x2
      }
      tweb <- t(web2)
      if (j == 1 & i == 1) 
        x3 <- 0
      else x3 <- (i - 1) * prey_spacing + cumsum(tweb)[(i - 
                                                          1) * nrow(tweb) + (j - 1)]/websum
      if (!is.null(low.abun2) && i > 1) 
        x3 <- x3 + cumsum(difff)[i - 1]
      x4 <- x3 + tweb[j, i]/websum
      if (arrow == "down" || arrow == "both") {
        x4 <- (x3 + x4)/2
        x3 <- x4
      }
      polygon(c(x1, x2, x4, x3), c(y1, y2, y4, y3), col = col.interaction2, border = bord.col.interaction2)
    }
  }
  #browser()
  # code below attempts to integrate interaction strengths of a third matrix
  if(!is.null(highlight.web.high)){
  pred_x <- 0
  #zwischenweb <- web2
  zwischenweb <- highlight.web.high
  web3 <- highlight.web.high
  XYcoords2 <- matrix(ncol = 2, nrow = length(zwischenweb))
  for (i in 1:length(zwischenweb)) {
    XYcoords2[i, ] <- which(zwischenweb == max(zwischenweb), 
                            arr.ind = TRUE)[1, ]
    zwischenweb[XYcoords2[i, 1], XYcoords2[i, 2]] <- -1
  }
  y1 <- pred_y - y_width - y_width_scale
  y2 <- y1
  y3 <- prey_y + y_width + y_width_scale
  y4 <- y3
  if (sum(web3 > 0)) { # web2 # web3
    for (p in 1:sum(web3 > 0)) { #web2 # web3
      i <- XYcoords2[p, 1]
      j <- XYcoords2[p, 2]
      if (j == 1 & i == 1) 
        x1 <- 0
      else x1 <- (j - 1) * pred_spacing + cumsum(web2)[(j - 
                                                          1) * nrow(web2) + (i - 1)]/websum
      if (!is.null(high.abun2) && j > 1) 
        x1 <- x1 + cumsum(diffh)[j - 1]
      x2 <- x1 + web3[i, j]/websum # adding web 3 so the portion of highlighted area is smaller
      if (arrow == "up" || arrow == "both") {
        x2 <- (x1 + x2)/2
        x1 <- x2
      }
      tweb <- t(web2)
      tweb3 <- t(web3)
      if (j == 1 & i == 1) 
        x3 <- 0
      else x3 <- (i - 1) * prey_spacing + cumsum(tweb)[(i - 
                                                          1) * nrow(tweb) + (j - 1)]/websum
      if (!is.null(low.abun2) && i > 1) 
        x3 <- x3 + cumsum(difff)[i - 1]
      x4 <- x3 + tweb3[j, i]/websum # changed to tweb3
      if (arrow == "down" || arrow == "both") {
        x4 <- (x3 + x4)/2
        x3 <- x4
      }
      polygon(c(x1, x2, x4, x3), c(y1, y2, y4, y3), col = col.interaction.highlight, border = bord.col.interaction.highlight, lwd = bord.size.interaction.highlight)
    }
  }
  }
  # replot middle for clear border
  prey_x <- 0
  links <- 0
  rechts <- 0
  hoehe <- strheight(rownames(web2)[1], cex = 0.6)
  hoffset <- hoehe
  for (i in 1:n.prey) {
    # attempt to add prey border highlight
    if(names(prey_prop[i]) %in% node.highlight)
      bord.col.prey = c(bord.col.highlight, bord.size.highlight)
    else {
      bord.col.prey = c(bord.col, bord.size)
    }
    if (!is.null(low.abun2)) {
      rect(prey_x + prey_prop[i] - difff[i], prey_y, prey_x + 
             prey_prop[i], prey_y + y_width, col = low.abun.col2, border = bord.low.abun.col2 )
    } 
    
    rect(prey_x, prey_y - y_width, prey_x + prey_prop[i], 
         prey_y + y_width, col = col.prey2, border = bord.col.prey[1], lwd = bord.col.prey[2])
    
    breite <- strwidth(rownames(web)[i], cex = 0.6 * labsize)
    links <- prey_x + prey_prop[i]/2 - breite/2
    if (links < rechts && i > 1) 
      hoffset = hoffset + hoehe
    else {
      rechts <- prey_x + prey_prop[i]/2 + breite/2
      hoffset <- hoehe
    }
    if(text == TRUE){
    text(prey_x + prey_prop[i]/2, prey_y - y_width - hoffset, 
         rownames(web2)[i], cex = 0.6 * labsize, offset = 0, font = font.type)
    }
    prey_x <- prey_x + prey_prop[i] + prey_spacing
  }
  # replot top for clear borders
  pred_x = 0
  hoffset <- 0
  links <- 0
  rechts <- 0
  hoehe <- strheight(colnames(web2)[1], cex = 0.6)
  for (i in 1:n.pred) {
    # attempt to add predator border highlight
    if(names(pred_prop[i]) %in% node.highlight)
      bord.col.pred = c(bord.col.highlight, bord.size.highlight)
    else {
      bord.col.pred = c(bord.col, bord.size)
    }
    rect(pred_x, pred_y - y_width, pred_x + pred_prop[i], 
         pred_y + y_width, col = col.pred2, border = bord.col.pred[1], lwd = bord.col.pred[2])
    if (!is.null(high.abun2)) {
      rect(pred_x + pred_prop[i] - diffh[i], pred_y - y_width, 
           pred_x + pred_prop[i], pred_y + y_width, col = high.abun.col2)
    }
    breite <- strwidth(colnames(web2)[i], cex = 0.6 * labsize)
    links <- pred_x + pred_prop[i]/2 - breite/2
    if (links < rechts && i > 1) 
      hoffset = hoffset + hoehe
    else {
      rechts <- pred_x + pred_prop[i]/2 + breite/2
      hoffset <- 0
    }
    if(text == TRUE){
    text(pred_x + pred_prop[i]/2, pred_y + y_width + hoehe + 
           hoffset, colnames(web2)[i], cex = 0.6 * labsize, 
         offset = 0, font = font.type)
    }
    pred_x <- pred_x + pred_prop[i] + pred_spacing
  }
}
