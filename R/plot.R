#' Plot a logic decision tree
#'
#' This function plots a logicDT model on the active graphics device.
#'
#' There are two plotting modes:
#' \itemize{
#'   \item \code{fancy = FALSE} which draws a tree with direct edges between
#'   the nodes. Leaves are represented by their prediction value which is
#'   obtained by the (observed) conditional mean.
#'   \item \code{fancy = TRUE} plots a tree similar to those in the \code{rpart}
#'   (Therneau and Atkinson, 2019) and \code{splinetree} (Neufeld and Heggeseth,
#'   2019) \code{R} packages. The trees are drawn in an angular manner and
#'   if leaf regression models were fitted, appropriate plots of the fitted
#'   curves are depicted in the leaves. Otherwise, the usual prediction values
#'   are shown.
#' }
#'
#' @param x An object of the class \code{logicDT}
#' @param fancy Should the fancy mode be used for plotting? Default is
#'   \code{TRUE}.
#' @param x_scaler Scaling factor on the horizontal axis for deeper trees,
#'   i.e., \code{x_scaler = 0.5} means that the horizontal distance
#'   between two adjacent nodes is halved for every vertical level.
#' @param margin_scaler Margin factor. Smaller values lead to smaller
#'   margins.
#' @param cex Scaling factor for the plotted text elements.
#' @param cdot Should a centered dot be used instead of a logical and
#'   for depicting interactions?
#' @param ... Arguments passed to fancy plotting function
#' @return No return value, called for side effects
#'
#' @references
#' \itemize{
#'   \item Therneau, T. & Atkinson, B. (2019). rpart: Recursive Partitioning
#'   and Regression Trees. \url{https://CRAN.R-project.org/package=rpart}
#'   \item Neufeld, A. & Heggeseth, B. (2019). splinetree: Longitudinal
#'   Regression Trees and Forests.
#'   \url{https://CRAN.R-project.org/package=splinetree}
#' }
#'
#' @importFrom graphics par plot lines rect text strwidth strheight
#' @export
plot.logicDT <- function(x, fancy = TRUE, x_scaler = 0.5, margin_scaler = 0.2,
                         cex = 1, cdot = FALSE, ...) {
  model <- x
  if(fancy) return(fancy.plot(model, cdot = cdot, ...))

  pet <- model$pet
  splits <- pet[[1]]
  splits_bin_or_cont <- pet[[2]]
  split_points <- pet[[3]]
  preds <- pet[[4]]

  X <- model$X
  Z <- model$Z
  X2 <- getDesignMatrix(X, model$disj)
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
  right_nodes_y <- integer()
  right_nodes_x <- integer()
  right_nodes <- integer()
  right_edges <- matrix(ncol = 2, nrow = 0)
  for(i in 1:number_of_nodes) {
    if(splits[i] != "<leaf>") {
      y_layers[i+1] <- y_layers[i] + 1
      x_layers[i+1] <- x_layers[i] - 1 * x_scaler^(y_layers[i+1])
      right_nodes_y <- c(right_nodes_y, y_layers[i] + 1)
      right_nodes_x <- c(right_nodes_x, x_layers[i] + 1 * x_scaler^(y_layers[i+1]))

      right_nodes <- c(right_nodes, i)
    } else {
      if(i+1 <= number_of_nodes) {
        y_layers[i+1] <- right_nodes_y[length(right_nodes_y)]
        right_nodes_y <- right_nodes_y[-length(right_nodes_y)]
        x_layers[i+1] <- right_nodes_x[length(right_nodes_x)]
        right_nodes_x <- right_nodes_x[-length(right_nodes_x)]

        right_edges <- rbind(right_edges, c(right_nodes[length(right_nodes)], i+1))
        right_nodes <- right_nodes[-length(right_nodes)]
      }
    }
  }

  max_y <- max(y_layers)
  y_coords <- max_y - y_layers
  x_margin <- 0.1 * (max(x_layers) - min(x_layers))
  y_margin <- 0.1 * (max(y_layers) - min(y_layers))
  x_range <- max(x_layers) - min(x_layers) + 2*x_margin
  y_range <- max(y_layers) + 2*y_margin
  plot(c(min(x_layers) - x_margin, max(x_layers) + x_margin), c(0-y_margin, max(y_layers)+y_margin), type = "n", axes = FALSE, xlab = "", ylab = "",
                 xaxs="i", yaxs="i")

  splits4 <- printable.splits(splits, splits_bin_or_cont, X, cdot = cdot)

  for(i in 1:number_of_nodes) {
    if(splits[i] != "<leaf>") {
      x <- x_layers[i]
      y <- max_y - y_layers[i]
      # txt <- as.character(splits[i])
      txt <- splits4[[i]]
      sw <- strwidth(txt, units = "user", cex = cex)
      sh <- strheight(txt, units = "user", cex = cex)

      lines(c(x_layers[i], x_layers[i+1]), c(y_coords[i], y_coords[i+1]))
      i2 <- right_edges[right_edges[,1] == i, 2]
      lines(c(x_layers[i], x_layers[i2]), c(y_coords[i], y_coords[i2]))

      sw <- 1.0 * sw + margin_scaler * x_scaler^max_y
      sh <- 1.0 * sh + margin_scaler
      frsz <- 0.0
      rect(x - sw/2 - frsz, y - sh/2 - frsz, x + sw/2 + frsz, y + sh/2 + frsz, col = "white")

      # min_sw <- strwidth(999, units = "user")
      # radius <- max(sw, min_sw) * 3/4
      # symbols(x, y, circles = radius, add = TRUE, inches = FALSE, bg = "white")

      text(x, y, labels = txt, cex = cex)

      inds <- c(i+1, i2)
      for(j in inds) {
        vec <- c(x_layers[j] - x_layers[i], y_coords[j] - y_coords[i])
        # Translate line vector to graphic units
        vec2 <- c(vec[1]/x_range, vec[2]/y_range) * par("pin")
        # Orthogonal for label placement
        orth <- c(vec2[2], -vec2[1])
        # Scale length such that the labels are placed right on the edges
        orth <- orth * 0.5*strheight("0", units = "inches", cex = cex)/sqrt(sum(orth^2))
        # Backtransformation
        orth <- c(orth[1]*x_range, orth[2]*y_range) / par("pin")
        # Tilt angle
        angle <- atan(vec2[2]/vec2[1]) * 180/pi
        if(splits_bin_or_cont[i] == 0) {
          left_label <- "0 ="
          right_label <- "= 1"
        } else {
          left_label <- bquote(.(sprintf("%.2f", split_points[i])) >= "")
          right_label <- bquote("" > .(sprintf("%.2f", split_points[i])))
        }
        if(j == i+1) {
          label <- left_label
          mult <- 1
        } else {
          label <- right_label
          mult <- -1
        }
        middle_point <- c(x_layers[i] + 0.5*vec[1], y_coords[i] + 0.5*vec[2])
        middle_point <- middle_point + mult * orth
        text(middle_point[1], middle_point[2],
                       labels = label, col = "black", srt = angle, cex = cex)
      }
    } else {
      x <- x_layers[i]
      y <- y_coords[i]
      txt <- sprintf("%.2f", preds[i])
      sw <- strwidth(txt, cex = cex)
      sh <- strheight(txt, cex = cex)
      sw <- 1.0 * sw + margin_scaler * x_scaler^max_y
      sh <- 1.0 * sh + margin_scaler
      frsz <- 0.0
      rect(x - sw/2 - frsz, y - sh/2 - frsz, x + sw/2 + frsz, y + sh/2 + frsz, col = "black")
      text(x, y, labels = txt, col = "white", cex = cex)
    }
  }
}

printable.splits <- function(splits, splits_bin_or_cont, X, cdot = FALSE) {
  splits2 <- strsplit(splits, "^", fixed=TRUE)
  splits3 <- lapply(splits2, function(x) ifelse(startsWith(x, "-"), substring(x, 2), x))
  splits3 <- lapply(splits3, function(x) ifelse(x == "<leaf>", "0", x))
  splitsm <- lapply(splits2, function(x) startsWith(x, "-"))
  splits3 <- lapply(seq_along(splits3), function(i) {
    if(!splits_bin_or_cont[i])
      colnames(X)[as.integer(splits3[[i]])]
    else
      splits3[[i]]
  })
  splits4 <- list()
  for(i in 1:length(splits3)) {
    splits4[[i]] <- splits3[[i]]
    if(length(splits3[[i]]) > 0) {
      for(j in 1:length(splits3[[i]])) {
        if(splitsm[[i]][j]) {
          splits4[[i]][j] <- paste(splits3[[i]][j], "^c", sep="")
          # splits4[[i]][j] <- paste("!", splits3[[i]][j], sep="")
          # splits4[[i]][j] <- tryCatch(parse(text=paste("symbol(\330)~ ", splits3[[i]][j], sep="")),
          #                             error=function(cond) paste("!", splits3[[i]][j], sep=""))
        }
      }
    }
  }
  if(cdot) {
    splits4 <- tryCatch(lapply(splits4, function(x) parse(text=paste(x, collapse = " %.% "))),
                        error=function(cond) lapply(splits4, function(x) parse(text=paste(x, collapse = " * "))))
  } else {
    splits4 <- tryCatch(lapply(splits4, function(x) parse(text=paste(x, collapse = " ~symbol(\331)~ "))),
                        error=function(cond) lapply(splits4, function(x) parse(text=paste(x, collapse = " ~AND~ "))))
  }
  splits4
}

# x11(width=15, height=10)
# plot.logicDT(model, fancy = FALSE, x_scaler = 0.525, margin_scaler = 0.3, cex = 0.5)

#' Plot calculated VIMs
#'
#' This function plots variable importance measures yielded by the function
#' \code{\link{vim}} in a dotchart.
#'
#' @param x An object of the class \code{vim}
#' @param p The number of most important terms which will be included
#'   in the plot. A value of 0 leads to plotting all terms.
#' @param ... Ignored additional parameters
#' @return No return value, called for side effects
#'
#' @importFrom graphics dotchart abline
#' @export
plot.vim <- function(x, p = 10, ...) {
  vim <- x$vims$vim
  names(vim) <- x$vims$var
  p <- min(p, length(vim))
  if(p > 0) vim <- vim[1:p]
  xlab <- paste(vim_def(x$scoring_rule), scoring_rule_name(x$scoring_rule))
  dotchart(sort(vim), xlab = xlab, xlim = c(0, max(vim)), main = vim_type_name(x$vim_type))
  abline(v = 0, lty = 2)
}

scoring_rule_name <- function(scoring_rule) {
  if(scoring_rule == "deviance") {
    "Deviance"
  } else if(scoring_rule == "brier") {
    "Brier Score"
  } else if(scoring_rule == "mis") {
    "Misclassification Rate"
  } else if(scoring_rule == "auc") {
    "Area under the Curve"
  } else if(scoring_rule == "nce") {
    "Normalized Cross Entropy"
  } else if(scoring_rule == "mse") {
    "Mean Squared Error"
  } else if(scoring_rule == "nrmse") {
    "Normalized Root Mean Squared Error"
  }
}

vim_def <- function(scoring_rule) {
  if(scoring_rule == "auc")
    "Mean Increase in the"
  else
    "Mean Decrease in the"
}

vim_type_name <- function(vim_type) {
  if(vim_type == "logic") {
    "Logic VIM"
  } else if(vim_type == "remove") {
    "Removal VIM"
  } else if(vim_type == "permutation") {
    "Permutation VIM"
  }
}
