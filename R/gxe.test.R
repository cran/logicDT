fit4plModelWithGroups <- function(y, Z, G) {
  if(any(!(y %in% c(0, 1)))) {
    y <- as.numeric(y)
  } else {
    y <- as.integer(y)
  }
  G <- as.integer(G)
  .Call(fit4plModelWithGroups_, y, Z, G)
}

predict.4pl.grouped <- function(object, Z, G, ...) {
  bcde <- object[[1]]
  y_bin <- object[[2]]
  beta <- object[[3]]
  n_groups <- object[[4]]
  ret <- bcde[2] + (bcde[3]-bcde[2])/(1+exp(bcde[1]*(Z-bcde[4])))
  G <- G + 1
  beta[n_groups] <- 0
  ret <- ret + beta[G]
  if(y_bin) {
    ret[ret > 1] <- 1
    ret[ret < 0] <- 0
  }
  return(ret)
}

#' @importFrom stats ave aggregate
fitLinearModelWithGroups <- function(y, Z, G) {
  y_bin <- !any(!(y %in% c(0, 1)))

  y.group.means <- ave(y, G, FUN = mean)
  Z.group.means <- ave(Z, G, FUN = mean)
  alpha <- sum((y.group.means - y) * Z)/sum((Z.group.means - Z) * Z)

  y.group.means <- aggregate(y, by = list(G), FUN = mean)
  Z.group.means <- aggregate(Z, by = list(G), FUN = mean)
  betas <- y.group.means[,2] - alpha * Z.group.means[,2]

  ret <- list(alpha = alpha, betas = betas, y_bin = y_bin)
  class(ret) <- "linear.grouped"
  return(ret)
}

predict.linear.grouped <- function(object, Z, G, ...) {
  alpha <- object$alpha
  betas <- object$betas
  y_bin <- object$y_bin
  ret <- alpha * Z + betas[G + 1]
  if(y_bin) {
    ret[ret > 1] <- 1
    ret[ret < 0] <- 0
  }
  return(ret)
}

fitLDAModelWithGroups <- function(y, Z, G) {
  pi.1 <- mean(y)
  pi.0 <- 1-mean(y)
  G.mat <- oneHotEncoder(G)
  X <- cbind(Z, G.mat)
  X.0 <- X[y == 0,]
  X.1 <- X[y == 1,]
  mu.0 <- colMeans(X.0)
  mu.1 <- colMeans(X.1)
  MU.0 <- matrix(mu.0, ncol = ncol(X.0), nrow = nrow(X.0), byrow = TRUE)
  MU.1 <- matrix(mu.1, ncol = ncol(X.1), nrow = nrow(X.1), byrow = TRUE)
  XTX.0 <- t(X.0 - MU.0) %*% (X.0 - MU.0)
  XTX.1 <- t(X.1 - MU.1) %*% (X.1 - MU.1)
  N <- length(y)
  sigma <- (XTX.0 + XTX.1) / (N - 2)
  siginv <- solve(sigma)

  alpha <- siginv %*% matrix(mu.1 - mu.0, ncol = 1)
  beta <- log(pi.1 / pi.0) - 0.5 * matrix(mu.1 + mu.0, nrow = 1) %*%
    siginv %*% matrix(mu.1 - mu.0, ncol = 1)
  beta <- as.numeric(beta)

  ret <- list(alpha = alpha, beta = beta)
  class(ret) <- "lda.grouped"
  return(ret)
}

predict.lda.grouped <- function(object, Z, G, ...) {
  alpha <- object$alpha
  beta <- object$beta

  G.mat <- oneHotEncoder(G)
  X <- cbind(Z, G.mat)

  ret <- as.numeric(X %*% alpha + beta)
  ret <- 1/(1 + exp(-ret))
  return(ret)
}

#' @importFrom stats model.matrix
oneHotEncoder <- function(G) {
  G.tmp <- as.factor(G)
  G.mat <- model.matrix(~0+G.tmp)
  G.mat <- G.mat[,-ncol(G.mat)]
  return(G.mat)
}

get.tree.scenario <- function(model, X) {
  pet <- model$pet
  splits <- pet[[1]]
  splits_bin_or_cont <- pet[[2]]
  split_points <- pet[[3]]
  preds <- pet[[4]]

  X.new <- as.matrix(X)
  mode(X.new) <- "integer"

  X <- model$X
  Z <- model$Z
  X_split_pos <- pet[[2]] == 0 & pet[[1]] > 0
  Z_split_pos <- pet[[2]] == 1
  leaf_pos <- pet[[1]] <= 0
  X_vars <- pet[[1]][X_split_pos]
  Z_vars <- pet[[1]][Z_split_pos]
  splits[leaf_pos] <- "<leaf>"

  number_of_nodes <- length(splits)
  right_nodes <- integer()
  right_edges <- matrix(ncol = 2, nrow = 0)
  left_edges <- matrix(ncol = 2, nrow = 0)
  for(i in 1:number_of_nodes) {
    if(splits[i] != "<leaf>") {
      right_nodes <- c(right_nodes, i)
      left_edges <- rbind(left_edges, c(i, i+1))
    } else {
      if(i+1 <= number_of_nodes) {
        right_edges <- rbind(right_edges, c(right_nodes[length(right_nodes)], i+1))
        right_nodes <- right_nodes[-length(right_nodes)]
      }
    }
  }

  left_edges_from <- left_edges[,1]
  left_edges_to <- left_edges[,2]
  right_edges_from <- right_edges[,1]
  right_edges_to <- right_edges[,2]
  dm <- getDesignMatrix(X.new, model$disj)
  splits[splits == "<leaf>"] <- 0
  splits <- as.integer(splits)

  ret <- vector(mode = "numeric", length = nrow(dm))

  for(i in 1:nrow(dm)) {
    node.ind <- 1
    while(TRUE) {
      go.right <- dm[i, splits[node.ind]]
      if(length(go.right) == 0) {
        ret[i] <- node.ind
        break
      } else if(go.right == 0) {
        start.point <- which(node.ind == left_edges_from)
        go.to <- left_edges_to[start.point]
        node.ind <- go.to
      } else if(go.right == 1) {
        start.point <- which(node.ind == right_edges_from)
        go.to <- right_edges_to[start.point]
        node.ind <- go.to
      }
    }
  }
  ret
}

#' Gene-environment interaction test
#'
#' Using a fitted \code{logicDT} model, a general GxE interaction
#' test can be performed.
#'
#' The testing is done by fitting one shared 4pL model
#' for all tree branches with different offsets, i.e., allowing main effects
#' of SNPs. This shared model is compared to the individual 4pL models fitted
#' in the \code{\link{logicDT}} procedure using a likelihood ratio test which
#' is asymptotically \eqn{\chi^2} distributed. The degrees of freedom are
#' equal to the difference in model parameters.
#' For regression tasks, alternatively, a F-test can be utilized.
#'
#' The shared 4pL model is given by
#' \deqn{Y = \tilde{f}(x, z, b, c, d, e, \beta_1, \ldots, \beta_{G-1})
#'   + \varepsilon = c + \frac{d-c}{1+\exp(b \cdot (x-e))}
#'   + \sum_{g=1}^{G-1} \beta_g \cdot 1(z = g) + \varepsilon}
#' with \eqn{z \in \lbrace 1, \ldots, G \rbrace} being a grouping variable,
#' \eqn{\beta_1, \ldots, \beta_{G-1}} being the offsets for the different
#' groups, and \eqn{\varepsilon} being a random error term.
#' Note that the last group \eqn{G} does not have an offset parameter, since
#' the model is calibrated such that the curve without any \eqn{\beta}'s
#' fits to the last group.
#'
#' The likelihood ratio test statistic is given by
#' \deqn{\Lambda = -2(\ell_{\mathrm{shared}} - \ell_{\mathrm{full}})}
#' for the log likelihoods of the shared and full 4pL models, respectively.
#' In the regression case, the test statistic can be calculated as
#' \deqn{\Lambda = N(\log(\mathrm{RSS}_{\mathrm{shared}}) -
#'   \log(\mathrm{RSS}_{\mathrm{full}}))}
#' with \eqn{\mathrm{RSS}} being the residual sum of squares for the
#' respective model.
#'
#' For regression tasks, the alternative F test statistic is given by
#' \deqn{f = \frac{\frac{1}{\mathrm{df}_1}(\mathrm{RSS}_{\mathrm{shared}} -
#'   \mathrm{RSS}_{\mathrm{full}})}
#'   {\frac{1}{\mathrm{df}_2} \mathrm{RSS}_{\mathrm{full}}}}
#' with
#' \deqn{\mathrm{df}_1 = \mathrm{Difference\ in\ the\ number\ of\ model\
#' parameters} = 3 \cdot n_{\mathrm{scenarios}} - 3,}
#' \deqn{\mathrm{df}_2 = \mathrm{Degrees\ of\ freedom\ of\ the\ full\ model}
#' = N - 4 \cdot n_{\mathrm{scenarios}},}
#' and \eqn{n_{\mathrm{scenarios}}} being the number of identified predictor
#' scenarios/groups by \code{\link{logicDT}}.
#'
#' Alternatively, if linear models were fitted in the supplied \code{logicDT}
#' model, shared linear models can be used to test for a GxE interaction.
#' For continuous outcomes, the shared linear model is given by
#' \deqn{Y = \tilde{f}(x, z, \alpha, \beta_1, \ldots, \beta_{G})
#'   + \varepsilon = \alpha \cdot x
#'   + \sum_{g=1}^{G} \beta_g \cdot 1(z = g) + \varepsilon.}
#' For binary outcomes, LDA (linear discriminant analysis) models are fitted.
#' In contrast to the 4pL-based test for binary outcomes, varying offsets for
#' the individual groups are injected to the linear predictor instead of to
#' the probability (response) scale.
#'
#' If only few samples are available and the asymptotics of likelihood ratio
#' tests cannot be justified, alternatively, a permutation test approach
#' can be employed by setting \code{perm.test = TRUE} and specifying an
#' appropriate number of random permutations via \code{n.perm}.
#' For this approach, computed likelihoods of the shared and (paired) full
#' likelihood groups are randomly interchanged approximating the null
#' distribution of equal likelihoods. A p-value can be computed by determining
#' the fraction of more extreme null samples compared to the original
#' likelihood ratio test statistic, i.e., using the fraction of higher
#' likelihood ratios in the null distribution than the original likelihood
#' ratio.
#'
#' @param model A fitted \code{logicDT} model with 4pL models in its leaves.
#' @param X Binary predictor data for testing the interaction effect.
#'   This can be equal to the training data.
#' @param y Response vector for testing the interaction effect.
#'   This can be equal to the training data.
#' @param Z Quantitative covariable for testing the interaction effect.
#'   This can be equal to the training data.
#' @param perm.test Should additionally permutation testing be performed?
#'   Useful if likelihood ratio test asymptotics cannot be justified.
#' @param n.perm Number of random permutations for permutation testing
#' @return A list containing
#'   \item{\code{p.chisq}}{The p-value of the chi-squared test statistic.}
#'   \item{\code{p.f}}{The p-value of the F test statistic.}
#'   \item{\code{p.perm}}{The p-value of the optional permutation test.}
#'   \item{\code{ll.shared}}{Log likelihood of the shared parameters 4pL model.}
#'   \item{\code{ll.full}}{Log likelihood of the full \code{logicDT} model.}
#'   \item{\code{rss.shared}}{Residual sum of squares of the shared parameters
#'     4pL model.}
#'   \item{\code{rss.full}}{Residual sum of squares of the full \code{logicDT}
#'     model.}
#' @importFrom stats pchisq pf
#' @export
gxe.test <- function(model, X, y, Z, perm.test = TRUE, n.perm = 10000) {
  n.scenarios <- sum(model$pet[[1]] == 0)
  y_bin <- model$y_bin
  N <- nrow(X)

  tree.scenarios <- get.tree.scenario(model, model$X)
  unique.scenarios <- unique(tree.scenarios)

  G <- match(tree.scenarios, unique.scenarios) - 1

  covariable_mode <- model$pet[[9]]
  if(covariable_mode < 2)
    return(NULL)

  if(covariable_mode == 2) {
    grouped.model <- fit4plModelWithGroups(model$y, model$Z, G)
    df <- abs(3 - 3 * n.scenarios)
    df2 <- abs(N - 4*n.scenarios)
  } else {
    if(!y_bin) {
      grouped.model <- fitLinearModelWithGroups(model$y, model$Z, G)
    } else {
      grouped.model <- fitLDAModelWithGroups(model$y, model$Z, G)
    }
    df <- abs(n.scenarios - 1)
    df2 <- abs(N - 2*n.scenarios)
  }
  tree.scenarios2 <- get.tree.scenario(model, X)
  G2 <- match(tree.scenarios2, unique.scenarios) - 1
  grouped.preds <- as.numeric(predict(grouped.model, Z, G2))

  full.preds <- predict(model, X = X, Z = Z, ensemble = FALSE)

  null.perm <- vector(mode = "numeric", length = n.perm)

  if(y_bin) {
    threshold <- log(0.01)

    ll.shared.vec <- log(y * grouped.preds + (1-y) * (1-grouped.preds))
    ll.shared.vec[ll.shared.vec < threshold] <- threshold
    ll.shared <- sum(ll.shared.vec)

    ll.full.vec <- log(y * full.preds + (1-y) * (1-full.preds))
    ll.full.vec[ll.full.vec < threshold] <- threshold
    ll.full <- sum(ll.full.vec)

    ll.shared.full <- cbind(ll.shared.vec, ll.full.vec)

    if(perm.test) {
      for(i in 1:n.perm) {
        perm <- sample(1:2, N, replace = TRUE)
        perm.shared <- ll.shared.full[cbind(1:N, perm)]
        perm.full <- ll.shared.full[cbind(1:N, 3-perm)]
        null.perm[i] <- -2 * (sum(perm.shared) - sum(perm.full))
      }
    }

    Q <- -2 * (ll.shared - ll.full)
    p.chisq <- 1 - pchisq(Q, df)

    if(perm.test) {
      p.perm <- mean(null.perm >= Q)
    } else {
      p.perm <- NULL
    }

    ret <- list(p.chisq = p.chisq, p.perm = p.perm,
                ll.full = ll.full, ll.shared = ll.shared)
  } else {
    rss.shared.vec <- (y - grouped.preds)^2
    rss.full.vec <- (y - full.preds)^2
    rss.shared <- sum(rss.shared.vec)
    rss.full <- sum(rss.full.vec)

    rss.shared.full <- cbind(rss.shared.vec, rss.full.vec)

    if(perm.test) {
      for(i in 1:n.perm) {
        perm <- sample(1:2, N, replace = TRUE)
        perm.shared <- rss.shared.full[cbind(1:N, perm)]
        perm.full <- rss.shared.full[cbind(1:N, 3-perm)]
        null.perm[i] <- N * log(sum(perm.shared)) - N * log(sum(perm.full))
      }
    }

    Q <- N * log(rss.shared) - N * log(rss.full)
    p.chisq <- 1 - pchisq(Q, df)

    if(perm.test) {
      p.perm <- mean(null.perm >= Q)
    } else {
      p.perm <- NULL
    }

    F <- ((rss.shared - rss.full)/df)/(rss.full/df2)
    p.f <- 1 - pf(F, df, df2)
    ret <- list(p.chisq = p.chisq, p.f = p.f, p.perm = p.perm,
                rss.full = rss.full, rss.shared = rss.shared)
  }
  return(ret)
}






