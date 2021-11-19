#' @importFrom utils combn
#' @importFrom stats fisher.test t.test
vim_adjusted <- function(model, scoring_rule, vim_type, interaction_order, nodesize, alpha,
                         X_oob, y_oob, Z_oob, leaves = "4pl", ...) {
  X <- model$X
  N <- nrow(X)
  disj <- model$disj
  n_conj <- sum(rowSums(!is.na(disj)) > 0)
  n_vars <- rowSums(!is.na(disj[1:n_conj,,drop=FALSE]))
  if(interaction_order < 2 || n_conj < 2) {
    return(NULL)
  }

  # Get all reasonable conjunctions
  new.conjs <- list()
  new.conj.vecs <- NULL
  for(i in 2:min(interaction_order, n_conj)) {
    tmp.conjs <- combn(n_conj, i, simplify = FALSE)
    # Check if new conjunctions are eligible
    # First, prune off interactions whose order is too high
    rem.conjs <- integer()
    tmp.disj <- NULL
    for(j in 1:length(tmp.conjs)) {
      conj.vec <- as.integer(disj[tmp.conjs[[j]],,drop=FALSE])
      conj.vec <- conj.vec[!is.na(conj.vec)]
      if(length(conj.vec) > interaction_order) {
        rem.conjs <- c(rem.conjs, j)
      } else {
        conj.vec  <- c(conj.vec, rep(NA_integer_, interaction_order - length(conj.vec)))
        tmp.disj <- rbind(tmp.disj, matrix(conj.vec, ncol = interaction_order, nrow = 1))
      }
    }
    if(length(rem.conjs) > 0) tmp.conjs <- tmp.conjs[-rem.conjs]

    if(is.null(tmp.disj)) {
      next
    }

    # Now consider the minimum node size
    mode(tmp.disj) <- "integer"
    dm <- getDesignMatrix(X, tmp.disj)
    nodesizes <- colSums(dm)
    rem.conjs <- which((nodesizes < nodesize) | (nodesizes > N - nodesize))
    if(length(rem.conjs) > 0) {
      tmp.conjs <- tmp.conjs[-rem.conjs]
      tmp.disj <- tmp.disj[-rem.conjs,,drop=FALSE]
    }

    new.conjs <- c(new.conjs, tmp.conjs)
    new.conj.vecs <- rbind(new.conj.vecs, tmp.disj)
  }

  vims <- data.frame(matrix(nrow = length(new.conjs), ncol = 2))
  colnames(vims) <- c("var", "vim")
  if(nrow(vims) == 0) return(vims)
  tmp.real.disj <- translateLogicPET(new.conj.vecs, X)
  vims$var <- getPredictorNames(tmp.real.disj, sort_conj = TRUE)
  vims$vim <- 0

  rem.conjs <- integer()
  for(i in 1:length(new.conjs)) {
    vim.vars <- list()
    rem.vars <- new.conjs[[i]]
    for(j in length(new.conjs[[i]]):1) {
      vim.vars <- c(vim.vars, combn(new.conjs[[i]], j, simplify = FALSE))
    }
    tmp.vims <- vim_single(model, scoring_rule, vim_type, vim.vars = vim.vars, rem.vars = rem.vars,
                           X_oob = X_oob, y_oob = y_oob, Z_oob = Z_oob, leaves = leaves, ...)
    for(j in 1:length(vim.vars)) {
      coefficient <- ((length(new.conjs[[i]]) - length(vim.vars[[j]])) %% 2) * (-2) + 1
      vims$vim[i] <- vims$vim[i] + coefficient * tmp.vims[j]
    }
  }

  # Post-hoc adjustment for identifying concrete conjunctions
  # (e.g., X_1^c & X_2 instead of X_1 & X_2)
  # Due to the issue of inverting conjunctions such as
  # (X_1 & X_2)^c & X_3, they can also be inverted as a whole
  keep.ind <- is.finite(vims$vim)
  vims <- vims[keep.ind,]

  if(alpha > 0 && nrow(vims) > 0) {
    new.conjs <- new.conjs[keep.ind]
    new.conj.vecs <- new.conj.vecs[keep.ind,,drop=FALSE]
    dm <- getDesignMatrix(X, model$disj)
    y <- model$y
    y.sum <- sum(y)
    name.pool <- getPredictorNames(model$real_disj, sort_conj = TRUE)

    for(i in 1:length(new.conjs)) {
      current.conj <- new.conjs[[i]]
      combs <- expand.grid(rep(list(0:1), length(current.conj)))
      dm.tmp <- dm[,current.conj,drop=FALSE]
      p.values <- rep(1.0, nrow(combs))
      for(j in 1:nrow(combs)) {
        # comb.rows <- unname(apply(dm.tmp, 1, function(x) all(x == combs[j,])))
        # Way faster:
        stacked <- matrix(rep(combs[j,], N), byrow = TRUE, nrow = N)
        comb.rows <- unname(rowMeans(dm.tmp == stacked)) == 1
        if(model$y_bin) {
          y.comb.sum <- sum(y[comb.rows])
          N.comb <- sum(comb.rows)
          # p.values[j] <- stats::prop.test(c(y.comb.sum, y.sum - y.comb.sum), c(N.comb, N - N.comb),
          #                                 alternative = "two.sided", correct = FALSE)$p.value
          # Fisher's exact test (allowing for small sample sizes/expected values < 5):
          p.values[j] <- fisher.test(matrix(c(y.comb.sum, y.sum - y.comb.sum,
                                              N.comb - y.comb.sum, N - N.comb - (y.sum - y.comb.sum)), ncol = 2),
                                     alternative = "two.sided", conf.int = FALSE)$p.value
        } else {
          p.values[j] <- t.test(y[comb.rows], y[!comb.rows], alternative = "two.sided",
                                paired = FALSE, var.equal = FALSE)$p.value
        }
      }
      if(min(p.values) < alpha) {
        j <- which.min(p.values)
        comb <- as.integer(combs[j,])
        new.vars <- character()
        for(k in 1:length(comb)) {
          current.name <- name.pool[current.conj[k]]
          # Negate term
          if(comb[k] == 0) {
            if(n_vars[current.conj[k]] > 1) {
              new.vars <- c(new.vars, paste("-(", current.name, ")", sep = ""))
            } else {
              new.vars <- c(new.vars, ifelse(startsWith(current.name, "-"),
                                             substr(current.name, 2, nchar(current.name)),
                                             paste("-", current.name, sep = "")))
            }
          } else {
            # Split non negated terms for proper sorting
            if(n_vars[current.conj[k]] > 1) {
              single.vars <- strsplit(current.name, "\\^")[[1]]
              new.vars <- c(new.vars, single.vars)
            } else {
              new.vars <- c(new.vars, current.name)
            }
          }
        }
        vims$var[i] <- paste(sort(new.vars), collapse = "^")
      }
    }
  }

  return(vims)
}

vim_single <- function(model, scoring_rule, vim_type, vim.vars, rem.vars,
                       X_oob, y_oob, Z_oob, leaves = "4pl", ...) {
  FUN <- getPerfFUN(vim_type)
  base.perf <- FUN(model, scoring_rule, rem.vars,
                   X_oob = X_oob, y_oob = y_oob, Z_oob = Z_oob, leaves = leaves, ...)
  if(is.list(vim.vars)) {
    vim <- numeric(length(vim.vars))
    for(i in 1:length(vim.vars)) {
      adj.rem.vars <- setdiff(rem.vars, vim.vars[[i]])
      if(length(adj.rem.vars) > 0) {
        vim[i] <- base.perf - FUN(model, scoring_rule, adj.rem.vars,
                                  X_oob = X_oob, y_oob = y_oob, Z_oob = Z_oob, leaves = leaves, ...)
      } else {
        vim[i] <- base.perf - perf(model, scoring_rule,
                                   X_oob = X_oob, y_oob = y_oob, Z_oob = Z_oob, leaves = leaves)
      }
    }
  } else {
    vim <- base.perf - perf(model, scoring_rule,
                            X_oob = X_oob, y_oob = y_oob, Z_oob = Z_oob, leaves = leaves)
  }
  return(vim)
}

perf <- function(model, scoring_rule,
                 X_oob, y_oob, Z_oob, leaves = "4pl") {
  orig_score <- 0
  n_folds <- length(model$X_train)
  for(j in 1:n_folds) {
    dm_orig <- getDesignMatrix(X_oob[[j]], model$disj)
    Z_temp <- NULL
    if(!is.null(Z_oob))
      Z_temp <- Z_oob[[j]]
    probs_orig <- predict_pet(model$ensemble[[j]], dm_orig, Z = Z_temp, "prob", leaves)

    orig_score <- orig_score + calcScore(probs_orig, y_oob[[j]], scoring_rule)
  }
  orig_score/n_folds
}

getPerfFUN <- function(vim_type) {
  if(vim_type == "logic")
    FUN <- perf.logic
  else if(vim_type == "remove")
    FUN <- perf.remove
  else if(vim_type == "permutation")
    FUN <- perf.permutation
  FUN
}

perf.logic <- function(model, scoring_rule, rem.vars,
                       X_oob, y_oob, Z_oob, leaves = "4pl",
                       average = "before") {
  n_folds <- length(model$X_train)
  score <- 0
  for(j in 1:n_folds) {
    Z_temp <- NULL
    if(!is.null(Z_oob))
      Z_temp <- Z_oob[[j]]
    dm_orig <- getDesignMatrix(X_oob[[j]], model$disj)
    combs <- expand.grid(rep(list(0:1), length(rem.vars)))
    probs <- 0
    tmp.score <- 0
    for(k in 1:nrow(combs)) {
      dm_inv <- dm_orig
      dm_inv[,rem.vars] <- rep(as.integer(combs[k,]), each = nrow(dm_orig))
      if(average == "before")
        probs <- probs + predict_pet(model$ensemble[[j]], dm_inv, Z = Z_temp, "prob", leaves)
      else {
        probs <- predict_pet(model$ensemble[[j]], dm_inv, Z = Z_temp, "prob", leaves)
        tmp.score <- tmp.score + calcScore(probs, y_oob[[j]], scoring_rule)
      }
    }
    if(average == "before") {
      probs <- probs/nrow(combs)
      score <- score + calcScore(probs, y_oob[[j]], scoring_rule)
    } else
      score <- score + tmp.score/nrow(combs)
  }
  score/n_folds
}

perf.remove <- function(model, scoring_rule, rem.vars,
                        X_oob, y_oob, Z_oob, leaves = "4pl",
                        empty.model = "mean") {
  n_folds <- length(model$X_train)
  disj <- model$disj
  disj <- disj[-rem.vars,,drop=FALSE]
  n_conj <- sum(rowSums(!is.na(disj)) > 0)
  if(n_conj == 0) {
    if(empty.model == "none")
      return(-Inf)
    else {
      score <- 0
      for(j in 1:n_folds) {
        score <- score + calcScore(rep(mean(model$y_train[[j]]), length(y_oob[[j]])), y_oob[[j]], scoring_rule)
      }
      return(score/n_folds)
    }
  }
  ensemble <- fitPETs(model$X_train, model$y_train, model$X_val, model$y_val, model$Z_train, model$Z_val, model$use_validation, model$y_bin,
                      model$tree_control$nodesize, model$tree_control$cp, model$tree_control$smoothing, model$tree_control$mtry, model$tree_control$covariable_final,
                      disj, n_conj, getScoreRule(scoring_rule), TRUE)
  score <- 0
  for(j in 1:n_folds) {
    Z_temp <- NULL
    if(!is.null(Z_oob))
      Z_temp <- Z_oob[[j]]
    dm <- getDesignMatrix(X_oob[[j]], disj)
    probs <- predict_pet(ensemble[[j]], dm, Z = Z_temp, "prob", leaves)
    score <- score + calcScore(probs, y_oob[[j]], scoring_rule)
  }
  score/n_folds
}

perf.permutation <- function(model, scoring_rule, rem.vars,
                             X_oob, y_oob, Z_oob, leaves = "4pl",
                             n.perm = 100) {
  n_folds <- length(model$X_train)
  disj <- model$disj
  score <- 0
  for(i in 1:n.perm) {
    for(j in 1:n_folds) {
      Z_temp <- NULL
      if(!is.null(Z_oob))
        Z_temp <- Z_oob[[j]]
      dm <- getDesignMatrix(X_oob[[j]], disj)
      for(k in 1:length(rem.vars)) {
        perm <- sample(nrow(X_oob[[j]]))
        dm[,rem.vars[k]] <- dm[perm, rem.vars[k]]
      }
      probs <- predict_pet(model$ensemble[[j]], dm, Z = Z_temp, "prob", leaves)
      score <- score + calcScore(probs, y_oob[[j]], scoring_rule)
    }
  }
  score/(n_folds * n.perm)
}

calcScore <- function(preds, y, scoring_rule) {
  if(scoring_rule == "deviance") {
    score <- calcDev(preds, y)
  } else if(scoring_rule == "brier") {
    score <- calcBrier(preds, y)
  } else if(scoring_rule == "mis") {
    score <- calcMis(preds, y)
  } else if(scoring_rule == "auc") {
    score <- -calcAUCFast(preds, y)
  } else if(scoring_rule == "nce") {
    score <- calcNCE(preds, y)
  } else if(scoring_rule == "mse") {
    score <- calcMSE(preds, y)
  } else if(scoring_rule == "nrmse") {
    score <- calcNRMSE(preds, y)
  }
  score
}

#' Variable Importance Measures (VIMs)
#'
#' Calculate variable importance measures (VIMs) based on different
#' approaches.
#'
#' Three different VIM methods are implemented:
#' \itemize{
#'   \item Permutation VIMs: Random permutations of the respective
#'     identified logic terms
#'   \item Removal VIMs: Removing single logic terms
#'   \item Logic VIMs: Prediction with both possible outcomes
#'     of a logic term
#' }
#' Details on the calculation of these VIMs are given below.
#'
#' By variable importance, importance of identified logic terms
#' is meant. These terms can also be single predictors but also
#' conjunctions in the spirit of this software package.
#'
#' @section Permutation VIMs:
#' Permutation VIMs are computed by comparing the the model's
#' performance using the original data and data with random
#' permutations of single terms. This approach was originally
#' proposed by Breiman & Cutler (2003).
#'
#' @section Removal VIMs:
#' Removal VIMs are constructed removing specific logic
#' term from the set of predictors, refitting the decision
#' tree and comparing the performance to the original model.
#' Thus, this approach requires that at least two terms were
#' found by the algorithm. Therefore, no VIM will be
#' calculated if \code{empty.model = "none"} was specified.
#' Alternatively, \code{empty.model = "mean"} can be set to
#' use the constant mean response model for approximating
#' the empty model.
#'
#' @section Logic VIMs:
#' Logic VIMs use the fact that Boolean conjunctions are
#' Boolean variables themselves and therefore are equal to
#' 0 or 1. To compute the VIM for a specific term,
#' predictions are performed once for this term fixed to
#' 0 and once for this term fixed to 1. Then, the arithmetic
#' mean of these two (risk or regression) predictions is
#' is used for calculating the performance. This performance
#' is then compared to the original one as in the other
#' VIM approaches (average = "before"). Alternatively,
#' predictions for each fixed 0-1 scenario of the considered
#' term can be performed leading to individual performances
#' which then are averaged and compared to the original
#' performance (average = "after").
#'
#' @section Validation:
#' Validation data sets which
#' were not used in the fitting of the model are prefered
#' preventing an overfitting of the VIMs themselves.
#' These should be specified by the \code{_oob} arguments,
#' if neither bagging nor inner validation was used for fitting
#' the model.
#'
#' @section Bagging:
#' For the bagging version, out of bag (OOB) data are naturally
#' used for the calculation of VIMs.
#'
#' @section VIM Adjustment for Interactions:
#' Since decision trees can naturally include interactions
#' between single predictors (especially when strong marginal
#' effects are present as well), logicDT models might, e.g.,
#' include the single input variables \eqn{X_1} and \eqn{X_2} but
#' not their interaction \eqn{X_1 \land X_2} although an interaction
#' effect is present. We, therefore, developed and implemented an
#' adjustment approach for calculating VIMs for such
#' unidentified interactions nonetheless.
#' For predictors \eqn{X_{i_1}, \ldots, X_{i_k} =: Z}, this interaction
#' importance is given by
#' \deqn{\mathrm{VIM}(X_{i_1} \land \ldots \land X_{i_k}) =
#' \mathrm{VIM}(X_{i_1}, \ldots, X_{i_k} \mid X \setminus Z) -
#' \sum_{\lbrace j_1, \ldots, j_l \rbrace {\subset \atop \neq}
#' \lbrace i_1, \ldots, i_k \rbrace}
#' \mathrm{VIM}(X_{j_1} \land \ldots \land X_{j_l} \mid X \setminus Z)}
#' and can basically be applied to all black-box models.
#' By \eqn{\mathrm{VIM}(A \mid X \setminus Z)}, the VIM of \eqn{A}
#' considering the predictor set excluding the variables in \eqn{Z}
#' is meant, i.e., the improvement of additionally considering \eqn{A}
#' while regarding only the predictors in \eqn{X \setminus Z}.
#' The proposed interaction VIM can be recursively calculated through
#' \deqn{\mathrm{VIM}(X_{i_1} \land X_{i_2}) =
#' \mathrm{VIM}(X_{i_1}, X_{i_2} \mid X \setminus Z) -
#' \mathrm{VIM}(X_{i_1} \mid X \setminus Z) -
#' \mathrm{VIM}(X_{i_2} \mid X \setminus Z)}
#' for \eqn{Z = X_{i_1}, X_{i_2}}.
#' This leads to the relationship
#' \deqn{\mathrm{VIM}(X_{i_1} \land \ldots \land X_{i_k}) =
#' \sum_{\lbrace j_1, \ldots, j_l \rbrace \subseteq \lbrace i_1, \ldots, i_k \rbrace}
#' (-1)^{k-l} \cdot \mathrm{VIM}(X_{j_1}, \ldots, X_{j_l} \mid X \setminus Z).}
#'
#' @section Identification of Concrete Conjunctions:
#' The aforementioned VIM adjustment approach only captures the importance
#' of a general definition of interactions, i.e., it just considers
#' the question whether some variables do interact in any way.
#' Since logicDT is aimed at identifying specific conjunctions (and also assigns
#' them VIMs if they were identified by \code{\link{logicDT}}), a further
#' adjustment approach is implemented which tries to identify the specific
#' conjunction leading to an interaction effect.
#' The idea of this method is to consider the response for each possible
#' scenario of the interacting variables, e.g., for \eqn{X_1 \land (X_2^c \land X_3)}
#' where the second term \eqn{X_2^c \land X_3} was identified by \code{\link{logicDT}}
#' and, thus, two interacting terms are regarded,
#' the \eqn{2^2 = 4} possible scenarios
#' \eqn{\lbrace (i, j) \mid i, j \in \lbrace 0, 1 \rbrace \rbrace}
#' are considered. For each setting, the corresponding response is compared with
#' outcome values of the complementary set. For continuous outcomes, a two sample
#' t-test (with Welch correction for potentially unequal variances) is performed
#' comparing the means between these two groups. For binary outcomes, Fisher's exact
#' test is performed testing different underlying case probabilities.
#' If at least one test rejects the null hypothesis of equal outcomes (without adjusting
#' for multiple testing), the combination with the lowest p-value is chosen as the
#' explanatory term for the interaction effect. For example, if the most significant
#' deviation results from \eqn{X_1 = 0} and \eqn{(X_2^c \land X_3) = 1} from the example
#' above, the term \eqn{X_1^c \land (X_2^c \land X_3)} is chosen.
#'
#' @param model The fitted \code{logicDT} or \code{logic.bagged}
#'   model
#' @param scoring_rule The scoring rule for assessing the model
#'   performance. As in \code{\link{logicDT}}, "auc", "nce",
#'   "deviance" and "brier" are possible for binary outcomes.
#'   For regression, the mean squared error is used.
#' @param vim_type The type of VIM to be calculated. This can
#'   either be \code{"logic"}, \code{"remove"} or
#'   \code{"permutation"}. See below for details.
#' @param adjust Shall adjusted interaction VIMs be additionally
#'   (to the VIMs of identified terms) computed? See below for
#'   details.
#' @param interaction_order If \code{adjust = TRUE}, up to which
#'   interaction order shall adjusted interaction VIMs be
#'   computed?
#' @param nodesize If \code{adjust = TRUE}, how many observations
#'   need to be discriminated by an interaction in order to being
#'   considered? Similar to \code{conjsize} in \code{\link{logicDT}}
#'   and \code{nodesize} in \code{\link{tree.control}}.
#' @param alpha If \code{adjust = TRUE}, a further adjustment can be
#'   performed trying to identify the concrete conjunctions responsible
#'   for the interaction of the considered binary predictors.
#'   \code{alpha} specifies the significance level for statistical tests
#'   testing the alternative of a difference in the response for specific
#'   conjunctions. \code{alpha = 0} leads to no further adjustment.
#'   See below for details.
#' @param X_oob The predictor data which should be used for
#'   calculating the VIMs.
#'   Preferably some type of validation
#'   data independent of the training data.
#' @param y_oob The outcome data for computing the VIMs.
#'   Preferably some type of validation
#'   data independent of the training data.
#' @param Z_oob The optional covariable data for computing the
#'   VIMs.
#'   Preferably some type of validation
#'   data independent of the training data.
#' @param leaves The prediction mode if 4pL models were fitted
#'   in the leaves. As in \code{\link{predict.logicDT}},
#'   "4pl" and "constant" are the possible settings.
#' @param ... Parameters passed to the different VIM type functions.
#'   For \code{vim_type = "logic"}, the argument \code{average} can
#'   be specified as \code{"before"} or \code{"after"}. For
#'   \code{vim_type = "permutation"}, \code{n.perm} can be set to
#'   the number of random permutations. See below for details.
#'   For \code{vim_type = "remove"}, \code{empty.model} can be specified
#'   as either \code{"none"} ignoring empty models with all predictive
#'   terms removed or \code{"mean"} using the response mean as prediction
#'   in the case of an empty model.
#' @return A data frame with two columns:
#'   \item{\code{var}}{Short descriptions of the terms for which the
#'     importance was measured. For example \code{-X1^X2} for
#'     \eqn{X_1^c \land X_2}.}
#'   \item{\code{vim}}{The actual calculated VIM values.}
#'   The rows of such a data frame are sorted decreasingly by the VIM values.
#'
#' @references
#' \itemize{
#'   \item Breiman, L. (2001). Random Forests. Machine Learning 45(1):5-32.
#'     \doi{https://doi.org/10.1023/A:1010933404324}
#'   \item Breiman, L. & Cutler, A. (2003). Manual on Setting Up, Using,
#'     and Understanding Random Forests V4.0. University of California,
#'     Berkeley, Department of Statistics.
#'     \url{https://www.stat.berkeley.edu/~breiman/Using_random_forests_v4.0.pdf}
#' }
#'
#' @export
vim <- function(model, scoring_rule = "auc", vim_type = "logic",
                adjust = TRUE, interaction_order = 3, nodesize = NULL, alpha = 0.05,
                X_oob = NULL, y_oob = NULL, Z_oob = NULL, leaves = "4pl", ...) {
  type <- class(model)
  if(type == "logicDT") {
    if(is.null(X_oob)) {
      X_oob <- model$X_val
      y_oob <- model$y_val
      Z_oob <- model$Z_val
    } else {
      XyZ <- prepareXyZ(X_oob, y_oob, Z_oob, model$y_bin, make.list = TRUE)
      X_oob <- XyZ$X; y_oob <- XyZ$y; Z_oob <- XyZ$Z
    }

    if(!model$y_bin && scoring_rule == "auc") scoring_rule <- "nrmse"

    if(is.null(nodesize)) nodesize <- model$conjsize
    disj <- model$disj
    n_conj <- sum(rowSums(!is.na(disj)) > 0)
    vims <- data.frame(matrix(nrow = n_conj, ncol = 2))
    colnames(vims) <- c("var", "vim")
    vims$var <- getPredictorNames(model$real_disj, sort_conj = TRUE)
    vims$vim <- 0
    for(i in 1:n_conj) {
      vims$vim[i] <- vim_single(model, scoring_rule, vim_type, NULL, i,
                                X_oob = X_oob, y_oob = y_oob, Z_oob = Z_oob, leaves = leaves, ...)
    }
    if(adjust) {
      vims <- rbind(vims, vim_adjusted(model, scoring_rule = scoring_rule, vim_type = vim_type,
                                       interaction_order = interaction_order, nodesize = nodesize, alpha = alpha,
                                       X_oob = X_oob, y_oob = y_oob, Z_oob = Z_oob, leaves = leaves, ...))
    }
  } else if(type == "logic.bagged") {
    models <- model$models
    bags <- model$bags

    XyZ <- prepareXyZ(model$X, model$y, model$Z, models[[1]]$y_bin, make.list = FALSE)
    X <- XyZ$X; y <- XyZ$y; Z <- XyZ$Z

    N <- nrow(X)
    bagging.iter <- length(models)

    vims <- data.frame(matrix(nrow = 0, ncol = 2))
    colnames(vims) <- c("var", "vim")

    if(is.null(nodesize)) nodesize <- models[[1]]$conjsize

    if(!models[[1]]$y_bin && scoring_rule == "auc") scoring_rule <- "nrmse"

    Z_temp <- NULL
    vim.iter <- bagging.iter

    for(i in 1:bagging.iter) {
      oob <- setdiff(1:N, bags[[i]])

      n_folds <- length(models[[i]]$ensemble)
      if(!is.null(Z))
        Z_temp <- rep(list(Z[oob,]), n_folds)
      ret <- vim(models[[i]], scoring_rule, vim_type, adjust = adjust,
                 interaction_order = interaction_order, nodesize = nodesize, alpha = alpha,
                 X_oob = rep(list(X[oob,]), n_folds), y_oob = rep(list(y[oob]), n_folds),
                 Z_oob = Z_temp, leaves = leaves, ...)$vims
      predictor_names <- ret$var
      current_vims <- ret$vim

      if(any(!is.finite(current_vims))) {
        vim.iter <- vim.iter - 1
        next
      }

      for(j in 1:length(predictor_names)) {
        if(predictor_names[j] %in% vims$var) {
          vims$vim[vims$var == predictor_names[j]] <- vims$vim[vims$var == predictor_names[j]] + current_vims[j]
        } else {
          vims <- rbind(vims, data.frame(var = predictor_names[j], vim = current_vims[j]))
        }
      }
    }
    vims$vim <- vims$vim/vim.iter
  } else {
    stop("VIMs can only be calculated for models of type logicDT or logic.bagged!")
  }
  vims <- vims[order(vims$vim, decreasing = TRUE),]
  rownames(vims) <- 1:nrow(vims)
  ret <- list(vims = vims, vim_type = vim_type, scoring_rule = scoring_rule)
  class(ret) <- "vim"
  return(ret)
}

#' @importFrom stats sd
permutation_test <- function(x1, x2, n_perm_t = 10000) {
  # x2 - x1
  n <- length(x1)
  t_perm <- vector()

  for(i in 0:n_perm_t) {
    if(i == 0)
      perm <- rep(FALSE, n)
    else
      perm <- sample(c(FALSE, TRUE), n, replace = TRUE)

    x1p <- c(x1[!perm], x2[perm])
    x2p <- c(x2[!perm], x1[perm])

    inner_rank_x1 <- rank(x1p)
    inner_rank_x2 <- rank(x2p)
    full_rank <- rank(c(x1p,x2p))
    full_rank_x1 <- full_rank[1:n]
    full_rank_x2 <- full_rank[(n+1):(2*n)]

    p_est <- 1/(2*n) * mean(full_rank_x2 - full_rank_x1) + 0.5
    Z_k <- 1/n * (full_rank_x2 - inner_rank_x2 - full_rank_x1 + inner_rank_x1)
    sd_est <- sd(Z_k)
    t_perm <- c(t_perm, sqrt(n) * (p_est - 0.5)/sd_est)
  }
  return(mean(t_perm[1] <= t_perm[-1]))
}

prepareXyZ <- function(X, y, Z, y_bin, make.list = TRUE) {
  if(!is.list(y)) {
    X <- as.matrix(X)
    mode(X) <- "integer"
    if(!y_bin) {
      y <- as.numeric(y)
    } else {
      y <- as.integer(y)
    }
    if(make.list) {
      X <- list(X)
      y <- list(y)
    }
    if(!is.null(Z)) {
      Z <- as.matrix(Z)
      mode(Z) <- "double"
      if(make.list) {
        Z <- list(Z)
      }
    }
  } else {
    n_folds <- length(X)
    for(i in 1:n_folds) {
      X[[i]] <- as.matrix(X[[i]])
      mode(X[[i]]) <- "integer"
      if(!y_bin) {
        y[[i]] <- as.numeric(y[[i]])
      } else {
        y[[i]] <- as.integer(y[[i]])
      }
      if(!is.null(Z)) {
        Z[[i]] <- as.matrix(Z[[i]])
        mode(Z[[i]]) <- "double"
      }
    }
  }
  return(list(X = X, y = y, Z = Z))
}



