#' @rdname plot.logicDT
#' @name plot.logicDT
#' @importFrom graphics par plot lines rect points text strwidth
#' @importFrom grDevices hcl.colors gray
#' @export
fancy.plot <- function(x, cdot = FALSE, ...) {
  model <- x
  # Preparing tree for plotting
  pet <- model$pet
  splits <- pet[[1]]
  splits_bin_or_cont <- pet[[2]]
  split_points <- pet[[3]]
  preds <- pet[[4]]
  pl.model <- pet[[7]]

  X <- model$X
  Z <- model$Z
  X2 <- getDesignMatrix(X, model$disj)
  Y <- model$y
  N <- length(Y)
  X_split_pos <- pet[[2]] == 0 & pet[[1]] > 0
  Z_split_pos <- pet[[2]] == 1
  leaf_pos <- pet[[1]] <= 0
  X_vars <- pet[[1]][X_split_pos]
  Z_vars <- pet[[1]][Z_split_pos]
  splits[X_split_pos] <- colnames(X2)[X_vars]
  splits[Z_split_pos] <- colnames(Z)[Z_vars]
  splits[leaf_pos] <- "<leaf>"

  number_of_nodes <- length(splits)
  x_layers <- rep(0, number_of_nodes)
  y_layers <- rep(0, number_of_nodes)
  node_inds <- rep(0, number_of_nodes)
  node_inds[1] <- 1
  obs_inds <- vector(mode = "list", length = number_of_nodes)
  obs_inds[[1]] <- as.integer(1:N)
  right_nodes_y <- integer()
  right_nodes_x <- integer()
  right_nodes <- integer()
  right_edges <- matrix(ncol = 2, nrow = 0)
  for(i in 1:number_of_nodes) {
    if(splits[i] != "<leaf>") {
      y_layers[i+1] <- y_layers[i] + 1
      right_nodes_y <- c(right_nodes_y, y_layers[i] + 1)

      right_nodes <- c(right_nodes, i)

      node_inds[i+1] <- 2 * node_inds[i]

      if(splits_bin_or_cont[i] == 0) {
        submat <- cbind(X2[obs_inds[[i]], splits[i], drop=FALSE], obs_inds[[i]])
        submat <- submat[submat[,1] == 0,]
        obs_inds[[i+1]] <- as.integer(submat[,2])
      } else {
        submat <- cbind(Z[obs_inds[[i]], splits[i], drop=FALSE], obs_inds[[i]])
        submat <- submat[submat[,1] <= split_points[i],]
        obs_inds[[i+1]] <- as.integer(submat[,2])
      }
    } else {
      if(i+1 <= number_of_nodes) {
        y_layers[i+1] <- right_nodes_y[length(right_nodes_y)]
        right_nodes_y <- right_nodes_y[-length(right_nodes_y)]

        j <- right_nodes[length(right_nodes)]
        right_edges <- rbind(right_edges, c(j, i+1))

        node_inds[i+1] <- 2 * node_inds[j] + 1

        if(splits_bin_or_cont[j] == 0) {
          submat <- cbind(X2[obs_inds[[j]], splits[j], drop=FALSE], obs_inds[[j]])
          submat <- submat[submat[,1] == 1,]
          obs_inds[[i+1]] <- as.integer(submat[,2])
        } else {
          submat <- cbind(Z[obs_inds[[j]], splits[j], drop=FALSE], obs_inds[[j]])
          submat <- submat[submat[,1] > split_points[j],]
          obs_inds[[i+1]] <- as.integer(submat[,2])
        }

        right_nodes <- right_nodes[-length(right_nodes)]
      }
    }
  }

  max_y <- max(y_layers)
  y_coords <- max_y - y_layers

  # Code inspired by rpart and splinetree

  # Plot tree skeleton
  y <- (1+max_y - y_layers) / max(y_layers, 4)
  x <- rep(0.0, number_of_nodes)
  x[leaf_pos] <- seq(sum(leaf_pos))
  left.child <- match(2 * node_inds, node_inds)
  right.child <- match(2 * node_inds + 1, node_inds)
  temp <- split(seq(node_inds)[!leaf_pos], y_layers[!leaf_pos])
  for(i in rev(temp)) {
    x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])
  }
  cxy <- par("cxy")
  margin.x <- 0.1
  margin.y <- 0.15
  temp1 <- range(x) + diff(range(x)) * c(-margin.x, margin.x)
  # Reserve enough space for leaf annotations at the bottom (leaf means)
  # Also some space at the top for the first horizontal edge
  temp2 <- range(y) + diff(range(y)) * c(-margin.y, 0) - c(cxy[2], -cxy[2]/5)
  plot(temp1, temp2, type = "n", axes = FALSE, xlab = "", ylab = "",
                 ylim = temp2, xlim = temp1,
                 xaxs="i", yaxs="i")
  is.left <- (node_inds %% 2L == 0L)
  node.left <- node_inds[is.left]
  parent <- match(node.left/2L, node_inds)
  sibling <- match(node.left + 1L, node_inds)
  branch <- 1
  temp <- (x[sibling] - x[is.left]) * (1 - branch)/2
  xx <- rbind(x[is.left], x[is.left] + temp,
              x[sibling] - temp, x[sibling], NA)
  yy <- rbind(y[is.left], y[parent], y[parent], y[sibling], NA)
  branch.col <- 1
  branch.lty <- 1
  branch.lwd <- 1
  lines(c(xx), c(yy), col = branch.col, lty = branch.lty,
                  lwd = branch.lwd)

  # Plot subplot frames
  leaves <- which(leaf_pos)
  limx <- list()
  limy <- list()
  for(i in 1:length(leaves)) {
    limx[[i]] <- x[leaves[[i]]]
    # 40% of (horizontal) white space between leaves for deep trees
    # or 7.5% of the whole horizontal plotting space for shallow trees
    if(0.4 < diff(range(x)) * 0.075) {
      limx[[i]] <- limx[[i]] + c(-0.4, 0.4)
    } else {
      limx[[i]] <- limx[[i]] + c(-1, 1) * diff(range(x)) * 0.075
    }
    limy[[i]] <- y[leaves[i]] + c(-2, 0) * diff(range(y)) * 0.075
    rect(limx[[i]][1], limy[[i]][1], limx[[i]][2], limy[[i]][2],
                   density = -1, col = gray(0.9))
  }

  if(!is.null(pl.model)) {
    # Plot 4pL curves and observed data
    # colors <- rainbow(length(leaves), v = 0.9)
    # colors <- rainbow(number_of_nodes, v = 0.9)
    colors <- hcl.colors(length(leaves), palette = "Dark 3")

    plotList.x <- list()
    plotList.y <- list()
    for(i in 1:length(leaves)) {
      plotList.x[[i]] <- Z[obs_inds[[leaves[i]]], 1]
      plotList.y[[i]] <- Y[obs_inds[[leaves[i]]]]
    }

    plotList2 <- list()
    for(i in 1:length(leaves)) {
      newx <- seq(min(Z), max(Z), length = 101)
      params <- pl.model[[leaves[i]]]
      if(is.null(params)) {
        plotList2[[i]] <- rep(preds[leaves[i]], length(newx))
      } else {
        plotList2[[i]] <- predict(params, newx)
      }
    }
    subYLim <- range(plotList.y, na.rm = TRUE)
    for(i in 1:length(leaves)) {
      # (Linearly) Scale y and x to the subplot boxes
      pts.y <- plotList.y[[i]]
      pts.y <- (pts.y - subYLim[1])/diff(subYLim)
      pts.y <- pts.y * diff(limy[[i]]) + limy[[i]][1]
      pts.x <- plotList.x[[i]]
      pts.x <- (pts.x - range(plotList.x)[1])/diff(range(plotList.x))
      pts.x <- pts.x * diff(limx[[i]]) + limx[[i]][1]
      if(!model$y_bin) {
        points(x = pts.x, pts.y, cex = 0.1, type = "p", col = "black", pch = 20)
      }

      pts <- plotList2[[i]]
      pts <- (pts - subYLim[1])/diff(subYLim)
      pts <- pts * diff(limy[[i]]) + limy[[i]][1]
      points(x = seq(limx[[i]][1], limx[[i]][2], length.out = length(pts)),
                       pts, lwd = 2, type = "l", col = colors[i])
    }

    cxy <- par("cxy")
    # Print leaf means
    for(i in 1:length(leaves)) {
      text(mean(limx[[i]]), limy[[i]][1] - 0.5 * cxy[2], col = "black", labels = sprintf("%.2f", preds[leaves[i]]))
    }
  } else {
    # Plot leaf means
    for(i in 1:length(leaves)) {
      text(mean(limx[[i]]), mean(limy[[i]]), col = "black", labels = sprintf("%.2f", preds[leaves[i]]))
    }
  }

  # Plot split variables
  splits2 <- printable.splits(splits, splits_bin_or_cont, X, cdot = cdot)
  cxy <- par("cxy"); cin <- par("cin")
  for(i in which(!leaf_pos)) {
    # Scale text appropriately
    available.width <- x[right.child[i]] - x[left.child[i]]
    # Text width in coordinates, inches -> points
    w <- strwidth(splits2[[i]], units = "inches") * cxy[1]/cin[1]
    cex.text <- min(available.width/w, 1)
    text(x[i], y[i] - 0.5 * cxy[2L], labels = splits2[[i]], cex = cex.text)
  }

  # Plot split points
  for(i in which(!leaf_pos)) {
    if(splits_bin_or_cont[i] == 0) {
      left_label <- "0 ="
      right_label <- "= 1"
    } else {
      left_label <- bquote(.(sprintf("%.2f", split_points[i])) >= "")
      right_label <- bquote("" > .(sprintf("%.2f", split_points[i])))
    }

    i2 <- right_edges[right_edges[,1] == i, 2]
    inds <- c(i+1, i2)
    for(j in inds) {
      if(j == i+1) {
        label <- left_label
        mult <- 1
        angle <- 90
      } else {
        label <- right_label
        mult <- -1
        angle <- -90
      }
      x_label <- x[j]
      y_label <- (y[i] + y[j])/2

      orth <- 0.5 * (par("cin")[2])#strheight("0", units = "inches")
      orth <- orth*diff(temp1) / (par("pin")[1])

      text(x_label + mult * orth, y_label,
                     labels = label, col = "black", srt = angle)
    }
  }
}

# x11(width=15, height=10)
# fancy.plot(model)

