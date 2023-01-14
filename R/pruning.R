#' Pruning path of a logic decision tree
#'
#' Using a single fitted logic decision tree, the cost-complexity pruning path
#' containing the ideal subtree for a certain complexity penalty can be
#' computed.
#'
#' This is mainly a helper function for \code{\link{cv.prune}} and should only
#' be used by the user if manual pruning is preferred.
#' More details are given in \code{\link{cv.prune}}.
#'
#' @param pet A fitted logic decision tree. This can be extracted from a
#'   \code{\link{logicDT}} model, e.g., using \code{model$pet}.
#' @param y Training outcomes for potentially refitting regression models in
#'   the leaves. This can be extracted from a \code{\link{logicDT}} model,
#'   e.g., using \code{model$y}.
#' @param Z Continuous training predictors for potentially refitting regression
#'   models in the leaves. This can be extracted from a \code{\link{logicDT}}
#'   model, e.g., using \code{model$Z}. If no continuous covariable was used in
#'   fitting the model, \code{Z = model$Z = NULL} should be specified.
#' @return Two lists. The first contains the sequence of complexity penalties
#'   \eqn{alpha}. The second list contains the corresponding logic decision
#'   trees which can then be substituted in an already fitted
#'   \code{\link{logicDT}} model, e.g., using
#'   \code{model$pet <- result[[2]][[i]]} where \code{result} is the returned
#'   object from this function and \code{i} is the chosen tree index.
#' @export
prune.path <- function(pet, y, Z) {
  .Call(prune_, pet, y, Z)
}

#' Post-pruning using a fixed complexity penalty
#'
#' Using a fitted \code{\link{logicDT}} model and a fixed complexity penalty
#' \code{alpha}, its logic decision tree can be (post-)pruned.
#'
#' Similar to Breiman et al. (1984), we implement post-pruning by first
#' computing the optimal pruning path and then choosing the tree that is
#' pruned according to the specified complexity penalty.
#'
#' If no validation data is available or if the tree shall be automatically
#' optimally pruned, \code{\link{cv.prune}} should be used instead which
#' employs k-fold cross-validation for finding the best complexity penalty
#' value.
#'
#' @param model A fitted \code{logicDT} model
#' @param alpha A fixed complexity penalty value. This value should be
#'   determined out-of-sample, e.g., performing hyperparameter optimization
#'   on independent validation data.
#' @param simplify Should the pruned model be simplified with regard to the
#'   input terms, i.e., should terms that are no longer in the tree contained
#'   be removed from the model?
#' @return The new \code{logicDT} model containing the pruned tree
#' @export
prune <- function(model, alpha, simplify = TRUE) {
  pruned.trees <- prune.path(model$pet, model$y, model$Z)

  current.alphas <- pruned.trees[[1]]
  ind <- length(current.alphas[current.alphas < alpha])
  if(ind == 0) ind <- 1
  best.pet <- pruned.trees[[2]][[ind]]
  model$pet <- best.pet

  if(simplify) {
    old_n_conj <- nrow(model$disj)
    old_disj <- model$disj
    model <- simple.simplifyDisjunctions(model)
    if(nrow(model$disj) < old_n_conj) {
      model$pet <- simplifyTree(model$pet, old_disj, model$disj)
      if(!is.null(model$ensemble) && length(model$ensemble) > 0) {
        for(i in 1:length(model$ensemble))
          model$ensemble[[i]] <- simplifyTree(model$ensemble[[i]], old_disj, model$disj)
      }
    }
  }
  return(model)
}

#' Optimal pruning via cross-validation
#'
#' Using a fitted \code{\link{logicDT}} model, its logic decision tree can be
#' optimally (post-)pruned utilizing k-fold cross-validation.
#'
#' Similar to Breiman et al. (1984), we implement post-pruning by first
#' computing the optimal pruning path and then using cross-validation for
#' identifying the best generalizing model.
#'
#' In order to handle continuous covariables with fitted regression models in
#' each leaf, similar to the likelihood-ratio splitting criterion in
#' \code{\link{logicDT}}, we propose using the log-likelihood as the impurity
#' criterion in this case for computing the pruning path.
#' In particular, for each node \eqn{t}, the weighted node impurity
#' \eqn{p(t)i(t)} has to be calculated and the inequality
#' \deqn{\Delta i(s,t) := i(t) - p(t_L | t)i(t_L) - p(t_R | t)i(t_R) \geq 0}
#' has to be fulfilled for each possible split \eqn{s} splitting \eqn{t} into
#' two subnodes \eqn{t_L} and \eqn{t_R}. Here, \eqn{i(t)} describes the
#' impurity of a node \eqn{t}, \eqn{p(t)} the proportion of data points falling
#' into \eqn{t}, and \eqn{p(t' | t)} the proportion of data points falling
#' from \eqn{t} into \eqn{t'}.
#' Since the regression models are fitted using maximum likelihood, the
#' maximum likelihood criterion fulfills this property and can also be seen as
#' an extension of the entropy impurity criterion in the case of classification
#' or an extension of the MSE impurity criterion in the case of regression.
#'
#' The default model selection is done by choosing the most parsimonious model
#' that yields a cross-validation error in the range of
#' \eqn{\mathrm{CV}_{\min} + \mathrm{SE}_{\min}}
#' for the minimal cross-validation error \eqn{\mathrm{CV}_{\min}} and its
#' corresponding standard error \eqn{\mathrm{SE}_{\min}}.
#' For a more robust standard error estimation, the scores are calculated per
#' training observation such that the AUC is no longer an appropriate choice
#' and the deviance or the Brier score should be used in the case of
#' classification.
#'
#' @param model A fitted \code{logicDT} model
#' @param nfolds Number of cross-validation folds
#' @param scoring_rule The scoring rule for evaluating the cross-validation
#'   error and its standard error. For classification tasks, \code{"deviance"}
#'   or \code{"Brier"} should be used.
#' @param choose Model selection scheme. If the model that minimizes the
#'   cross-validation error should be chosen, \code{choose = "min"} should be
#'   set. Otherwise, \code{choose = "1se"} leads to simplest model in the range
#'   of one standard error of the minimizing model.
#' @param simplify Should the pruned model be simplified with regard to the
#'   input terms, i.e., should terms that are no longer in the tree contained
#'   be removed from the model?
#' @return A list containing
#'   \item{\code{model}}{The new \code{logicDT} model containing the optimally
#'   pruned tree}
#'   \item{\code{cv.res}}{A data frame containing the penalties, the
#'   cross-validation scores and the corresponding standard errors}
#'   \item{\code{best.beta}}{The ideal penalty value}
#' @references
#' \itemize{
#'   \item Breiman, L., Friedman, J., Stone, C. J. & Olshen, R. A. (1984).
#'   Classification and Regression Trees. CRC Press.
#'   \doi{https://doi.org/10.1201/9781315139470}
#' }
#' @importFrom stats sd
#' @export
cv.prune <- function(model, nfolds = 10, scoring_rule = "deviance", choose = "1se", simplify = TRUE) {
  X <- model$X; y <- model$y; Z <- model$Z
  use_Z <- !is.null(Z)

  X_train <- list(); y_train <- list()
  X_val <- list(); y_val <- list()
  if(use_Z) {
    Z_train <- list(); Z_val <- list()
  } else {
    Z_train <- NULL; Z_val <- NULL
  }

  folds <- cut(seq(1, nrow(X)), breaks=nfolds, labels=FALSE)
  shuffle <- sample(nrow(X))
  X_shuffle <- X[shuffle,]
  y_shuffle <- y[shuffle]
  if(use_Z) Z_shuffle <- Z[shuffle,,drop=FALSE]

  for (i in 1:nfolds) {
    test_ind <- which(folds == i, arr.ind = TRUE)
    train_ind <- -test_ind
    X_train[[i]] <- X_shuffle[train_ind,,drop=FALSE]
    y_train[[i]] <- y_shuffle[train_ind]
    X_val[[i]] <- X_shuffle[test_ind,,drop=FALSE]
    y_val[[i]] <- y_shuffle[test_ind]
    if(use_Z) {
      Z_train[[i]] <- Z_shuffle[train_ind,,drop=FALSE]
      Z_val[[i]] <- Z_shuffle[test_ind,,drop=FALSE]
    }
  }

  if(!model$y_bin) scoring_rule <- "mse"

  pets <- fitPETs(X_train, y_train, X_train, y_train, Z_train, Z_train, FALSE, model$y_bin, model$tree_control$nodesize, model$tree_control$split_criterion, model$tree_control$alpha, model$tree_control$cp, model$tree_control$smoothing, model$tree_control$mtry, model$tree_control$covariable_final, model$disj, sum(rowSums(!is.na(model$disj)) > 0), getScoreRule(scoring_rule), model$gamma, return_full_model = TRUE)

  pruned.trees <- prune.path(model$pet, model$y, model$Z)
  alphas <- unique(pruned.trees[[1]])
  betas <- alphas
  alphas[length(alphas) + 1] <- Inf
  for(i in 1:(length(alphas) - 1)) {
    betas[i] <- sqrt(alphas[i] * alphas[i+1])
  }

  scores <- list()
  for(i in 1:nfolds) {
    current.pruned.trees <- prune.path(pets[[i]], y_train[[i]], Z_train[[i]])
    current.alphas <- current.pruned.trees[[1]]
    for(j in 1:length(betas)) {
      if(i == 1) scores[[j]] <- list()
      ind <- length(current.alphas[current.alphas < betas[j]])
      if(ind == 0) ind <- 1
      current.pet <- current.pruned.trees[[2]][[ind]]
      X.tmp <- as.matrix(X_val[[i]]); mode(X.tmp) <- "integer"
      dm <- getDesignMatrix(X.tmp, model$disj)
      preds <- predict_pet(current.pet, dm, Z = Z_val[[i]])
      scores[[j]][[i]] <- calcScorePerObservation(preds, y_val[[i]], scoring_rule)
    }
  }

  cv.res <- data.frame()
  scores <- lapply(scores, function(x) unlist(x))
  for(j in 1:length(scores)) {
    m <- mean(scores[[j]])
    se <- sd(scores[[j]])/sqrt(length(scores[[j]]))
    cv.res <- rbind(cv.res, data.frame(beta = betas[j], score = m, se = se, score.plus.1se = m + se))
  }

  min.ind <- max(which(cv.res$score == min(cv.res$score)))
  if(choose != "min") {
    max.val <- cv.res$score.plus.1se[min.ind]
    min.ind <- max(which(cv.res$score <= max.val))
  }
  best.beta <- betas[min.ind]

  model <- prune(model, best.beta, simplify = simplify)

  return(list(model = model, cv.res = cv.res, best.beta = best.beta))
}

simple.simplifyDisjunctions <- function(model) {
  splits <- model$pet[[1]]
  n_conj <- sum(rowSums(!is.na(model$disj)) > 0)
  to.remove <- setdiff(1:n_conj, unique(splits))
  # Keep at least one term:
  if(length(to.remove) == n_conj) to.remove <- to.remove[-1]
  if(length(to.remove) == 0) return(model)
  model$disj <- model$disj[-to.remove,,drop=FALSE]
  model$disj <- simplifyMatrix(model$disj)
  model$disj <- dontNegateSinglePredictors(model$disj)
  model$real_disj <- translateLogicPET(model$disj, model$X)
  return(model)
}

calcScorePerObservation <- function(preds, y, scoring_rule) {
  scores <- vector(mode = "numeric", length = length(preds))
  for(i in 1:length(preds)) scores[i] <- calcScore(preds[i], y[i], scoring_rule)
  scores
}

simplifyTree <- function(pet, old_disj, new_disj) {
  trafo <- rep(0, nrow(old_disj))
  for(i in 1:nrow(new_disj)) {
    old_ind <- which(apply(old_disj, 1, function(x) all.equal(x, new_disj[i,])) == "TRUE")
    trafo[old_ind] <- i - 1
  }
  .Call(simplifyTree_, pet, as.integer(trafo))
}

