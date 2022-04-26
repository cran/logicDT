#' @useDynLib logicDT fitPETs_
#' @useDynLib logicDT fitPET_
#' @useDynLib logicDT predict_
#' @useDynLib logicDT predictEnsemble_
#' @useDynLib logicDT vim_permutation_
#' @useDynLib logicDT getDesignMatrix_
#' @useDynLib logicDT geneticProgramming_
#' @useDynLib logicDT predictGP_
#' @useDynLib logicDT simulatedAnnealing_
#' @useDynLib logicDT greedySearch_
#' @useDynLib logicDT calcAUC_
#' @useDynLib logicDT fit4plModel_
#' @useDynLib logicDT fit4plModelWithGroups_
#' @useDynLib logicDT fitLinearModel_

fitPETs <- function(X_train, y_train, X_val, y_val, Z_train, Z_val, use_validation, y_bin, nodesize, split_criterion, alpha, cp, smoothing, mtry, covariable_mode, disj, real_n_conj, scoring_rule, gamma, return_full_model = FALSE) {
  return(.Call(fitPETs_, X_train, y_train, X_val, y_val, Z_train, Z_val, use_validation, y_bin, nodesize, split_criterion, alpha, cp, smoothing, mtry, covariable_mode, disj, real_n_conj, as.integer(scoring_rule), gamma, return_full_model))
}

fitPET <- function(X, y, Z = NULL, nodesize, split_criterion, alpha, cp, smoothing, mtry, covariable_mode) {
  pet <- .Call(fitPET_, X, y, Z, nodesize, split_criterion, alpha, cp, smoothing, mtry, covariable_mode)
  return(pet)
}

predict_pet <- function(pet, X, Z = NULL, type = "prob", leaves = "4pl") {
  if (type == "prob")
    typus <- FALSE
  else
    typus <- TRUE
  if (leaves == "4pl")
    leaves <- 1
  else
    leaves <- 0
  if(is.null(Z) & any(pet[[2]] == 1))
    stop("The covariables are missing!")
  preds <- .Call(predict_, pet, X, Z, typus, leaves)
  return(preds)
}

predict_ensemble <- function(ensemble, X, Z = NULL, type = "prob", leaves = "4pl") {
  if (type == "prob")
    typus <- FALSE
  else
    typus <- TRUE
  if (leaves == "4pl")
    leaves <- 1
  else
    leaves <- 0
  preds <- .Call(predictEnsemble_, ensemble, X, Z, typus, leaves)
  return(preds)
}

#' @export
logicDT <- function(X, ...) UseMethod("logicDT")

#' Fitting logic decision trees
#'
#' Main function for fitting logicDT models.
#'
#' logicDT is a method for finding response-associated interactions between
#' binary predictors. A global search for the best set of predictors and
#' interactions between predictors is performed trying to find the global
#' optimal decision trees. On the one hand, this can be seen as a variable
#' selection. On the other hand, Boolean conjunctions between binary predictors
#' can be identified as impactful which is particularly useful if the
#' corresponding marginal effects are negligible due to the greedy fashion of
#' choosing splits in decision trees.
#'
#' Three search algorithms are implemented:
#' \itemize{
#'   \item Simulated annealing. An exhaustive stochastic optimization procedure.
#'         Recommended for single models (without [outer] bagging or boosting).
#'   \item Greedy search. A very fast search always looking for the best
#'         possible improvement. Recommended for ensemble models.
#'   \item Genetic programming. A more or less intensive search holding
#'         several competetive models at each generation. Niche method
#'         which is only recommended if multiple (simple) models do
#'         explain the variation in the response.
#' }
#'
#' Furthermore, the option of a so-called "inner validation" is available.
#' Here, the search is guided using several train-validation-splits and
#' the average of the validation performance. This approach is
#' computationally expensive but can lead to more robust single models.
#'
#' For minimizing the computation time, two-dimensional hash tables
#' are used saving evaluated models. This is irrelevant for the greedy
#' search but can heavily improve the fitting times when employing
#' a search with simulated annealing or genetic programming, especially
#' when choosing an inner validation.
#'
#' @section Saving and Loading:
#' logicDT models can be saved and loaded using \code{save(...)} and
#'   \code{load(...)}. The internal C structures will not be saved
#'   but rebuilt from the R representations if necessary.
#'
#' @param X Matrix or data frame of binary predictors coded as 0 or 1.
#' @param y Response vector. 0-1 coding for binary responses.
#'  Otherwise, a regression task is assumed.
#' @param max_vars Maximum number of predictors in the set of predictors.
#'   For the set \eqn{[X_1 \land X_2^c, X_1 \land X_3]}, this parameter
#'   is equal to 4.
#' @param max_conj Maximum number of input variables for the decision
#'   trees. For the set \eqn{[X_1 \land X_2^c, X_1 \land X_3]}, this
#'   parameter is equal to 2.
#' @param Z Optional matrix or data frame of quantitative/continuous
#'   covariables. Multiple covariables allowed for splitting the trees.
#'   If four parameter logistic models shall be fitted in the leaves,
#'   only the first given covariable is used.
#' @param search_algo Search algorithm for guiding the global search.
#'   This can either be "sa" for simulated annealing, "greedy" for a
#'   greedy search or "gp" for genetic programming.
#' @param cooling_schedule Cooling schedule parameters if simulated
#'   annealing is used. The required object should be created via
#'   the function \code{\link{cooling.schedule}}.
#' @param scoring_rule Scoring rule for guiding the global search.
#'   This can either be "auc" for the area under the receiver
#'   operating characteristic curve (default for binary reponses),
#'   "deviance" for the deviance, "nce" for the normalized cross entropy
#'   or "brier" for the Brier score.
#'   For regression purposes, the MSE (mean squared error) is
#'   automatically chosen.
#' @param tree_control Parameters controlling the fitting of
#'   decision trees. This should be configured via the
#'   function \code{\link{tree.control}}.
#' @param gamma Complexity penalty added to the score.
#'   If \eqn{\texttt{gamma} > 0} is given, \eqn{\texttt{gamma} \cdot ||m||_0}
#'   is added to the score with \eqn{||m||_0} being the total number of
#'   variables contained in the current model \eqn{m}.
#'   The main purpose of this penalty is for fitting logicDT stumps
#'   in conjunction with boosting. For regular logicDT models or bagged
#'   logicDT models, instead, the model complexity parameters \code{max_vars}
#'   and \code{max_conj} should be tuned.
#' @param simplify Should the final fitted model be simplified?
#'   This means, that unnecessary terms as a whole ("conj") will
#'   be removed if they cannot improve the score.
#'   simplify = "vars" additionally tries to prune individual
#'   conjunctions by removing unnecessary variables in those.
#'   simplify = "none" will not modify the final model.
#' @param val_method Inner validation method. "rv" leads to a
#'   repeated validation where val_reps times the original data
#'   set is divided into \eqn{\texttt{val\_frac} \cdot 100\%} validation
#'   data and \eqn{(1-\texttt{val\_frac}) \cdot 100\%} training data.
#'   "bootstrap" draws bootstrap samples and uses the out-of-bag
#'   data as validation data. "cv" employs cross-validation with
#'   val_reps folds.
#' @param val_frac Only used if val_method = "rv". See description
#'   of val_method.
#' @param val_reps Number of inner validation partitionings.
#' @param allow_conj_removal Should it be allowed to remove
#'   complete terms/conjunctions in the search?
#'   If exact numbers of terms are aimed at, this should be set
#'   to FALSE. If extensive hyperparameter optimizations are
#'   feasible, allow_conj_removal = FALSE with a proper search
#'   over max_vars and max_conj is advised for fitting single
#'   models. For bagging or boosting with a greedy search,
#'   allow_conj_removal = TRUE together with a small number for
#'   max_vars = max_conj is recommended, e.g., 2 or 3.
#' @param conjsize The minimum of training samples that have to
#'   belong to a conjunction. This parameters prevents including
#'   unnecessarily complex conjunctions that rarely occur.
#' @param randomize_greedy Should the greedy search be randomized
#'   by only considering \eqn{\sqrt{\mathrm{Neighbour\ states}}}
#'   neighbors at each iteration, similar to random forests.
#'   Speeds up the greedy search but can lead to inferior results.
#' @param greedy_mod Should modifications of conjunctions be
#'   considered in a greedy search?
#'   Speeds up the greedy search but can lead to inferior results.
#' @param greedy_rem Should the removal of conjunctions be
#'   considered in a greedy search?
#'   Speeds up the greedy search but can lead to inferior results.
#' @param max_gen Maximum number of generations for genetic
#'   programming.
#' @param gp_sigma Parameter \eqn{\sigma} for fitness sharing in
#'   genetic programming. Very small values (e.g., 0.001) are
#'   recommended leading to only penalizing models with yield
#'   the exact same score.
#' @param gp_fs_interval Interval for fitness sharing in
#'   genetic programming. The fitness calculation can be
#'   computationally expensive if many models exist in one
#'   generation. gp_fs_interval = 10 leads to performing
#'   fitness sharing only every 10th generation.
#' @return An object of class \code{logicDT}. This is a list
#'   containing
#'   \item{\code{disj}}{A matrix of identified set of predictors
#'     and conjunctions of predictors. Each entry corresponds to
#'     the column index in X. Negative values indicate negations.
#'     Missing values mean that the term does not contain any
#'     more variables.}
#'   \item{\code{real_disj}}{Human readable form of \code{disj}.
#'     Here, variable names are directly depicted.}
#'   \item{\code{score}}{Score of the best model. Smaller values
#'     are prefered.}
#'   \item{\code{pet}}{Decision tree fitted on the best set of
#'     input terms. This is a list containing the pointer to the
#'     C representation of the tree and R representations of the
#'     tree structure such as the splits and predictions.}
#'   \item{\code{ensemble}}{List of decision trees. Only relevant
#'     if inner validation was used.}
#'   \item{\code{total_iter}}{The total number of search
#'     iterations, i.e., tested configurations by fitting a tree
#'     (ensemble) and evaluating it.}
#'   \item{\code{prevented_evals}}{The number of prevented tree
#'     fitting by using the 2-dimensional hash table.}
#'   \item{\code{...}}{Supplied parameters of the functional call
#'     to \code{\link{logicDT}}.}
#' @example examples/toy.R
#' @references
#' \itemize{
#'   \item Lau, M., Schikowski, T. & Schwender, H. (2021).
#'   logicDT: A Procedure for Identifying Response-Associated
#'   Interactions Between Binary Predictors. To be submitted.
#'   \item Breiman, L., Friedman, J., Stone, C. J. & Olshen, R. A. (1984).
#'   Classification and Regression Trees. CRC Press.
#'   \doi{https://doi.org/10.1201/9781315139470}
#'   \item Kirkpatrick, S., Gelatt C. D. & Vecchi M. P. (1983).
#'   Optimization by Simulated Annealing. Science 220(4598):671–680.
#'   \doi{https://doi.org/10.1126/science.220.4598.671}
#' }
#'
#' @name logicDT
#' @method logicDT default
#' @importFrom stats var quantile
#' @export
logicDT.default <- function(X, y, max_vars = 3, max_conj = 3, Z = NULL,
                            search_algo = "sa",
                            cooling_schedule = cooling.schedule(),
                            scoring_rule = "auc", tree_control = tree.control(),
                            gamma = 0,
                            simplify = "vars",
                            val_method = "none", val_frac = 0.5, val_reps = 10,
                            allow_conj_removal = TRUE, conjsize = 1,
                            randomize_greedy = FALSE, greedy_mod = TRUE, greedy_rem = FALSE,
                            max_gen = 10000, gp_sigma = 0.15, gp_fs_interval = 1, ...) {

  total_iter <- 0
  prevented_evals <- 0

  p <- ncol(X)

  X <- as.matrix(X)
  mode(X) <- "integer"
  # if(setequal(unique(y), 0:1))
  if(any(!(y %in% c(0, 1)))) {
    y <- as.numeric(y)
    y_bin <- FALSE

    if(tree_control$cp_orig > 0) {
      tree_control$cp <- tree_control$cp_orig * var(y)
    }
  } else {
    y <- as.integer(y)
    y_bin <- TRUE
  }

  X_train <- list()
  y_train <- list()

  X_val <- list()
  y_val <- list()

  if(!is.null(Z)) {
    Z <- as.matrix(Z)
    pZ <- ncol(Z)
    mode(Z) <- "double"
    Z_train <- list()
    Z_val <- list()
    use_Z <- TRUE
  } else {
    pZ <- 0
    Z_train <- NULL
    Z_val <- NULL
    use_Z <- FALSE
  }

  if (val_method %in% c("bootstrap", "rv", "cv") & val_reps > 0) {
    if (val_method == "cv")
      folds <- cut(seq(1, nrow(X)), breaks=val_reps, labels=FALSE)

    for (i in 1:val_reps) {
      if (val_method == "bootstrap") {
        train_ind <- sample(1:nrow(X), nrow(X), replace = TRUE)
        test_ind <- setdiff(1:nrow(X), train_ind)
        X_train[[i]] <- X[train_ind,,drop=FALSE]
        y_train[[i]] <- y[train_ind]
        X_val[[i]] <- X[test_ind,,drop=FALSE]
        y_val[[i]] <- y[test_ind]
        if(use_Z) {
          Z_train[[i]] <- Z[train_ind,,drop=FALSE]
          Z_val[[i]] <- Z[test_ind,,drop=FALSE]
        }
      } else if (val_method == "rv") {
        N <- floor((1-val_frac)*nrow(X))
        train_ind <- sample(1:nrow(X), N)
        test_ind <- -train_ind
        X_train[[i]] <- X[train_ind,,drop=FALSE]
        y_train[[i]] <- y[train_ind]
        X_val[[i]] <- X[test_ind,,drop=FALSE]
        y_val[[i]] <- y[test_ind]
        if(use_Z) {
          Z_train[[i]] <- Z[train_ind,,drop=FALSE]
          Z_val[[i]] <- Z[test_ind,,drop=FALSE]
        }
      } else {
        shuffle <- sample(nrow(X))
        X_shuffle <- X[shuffle,]
        y_shuffle <- y[shuffle]
        test_ind <- which(folds == i, arr.ind = TRUE)
        train_ind <- -test_ind
        X_train[[i]] <- X_shuffle[train_ind,,drop=FALSE]
        y_train[[i]] <- y_shuffle[train_ind]
        X_val[[i]] <- X_shuffle[test_ind,,drop=FALSE]
        y_val[[i]] <- y_shuffle[test_ind]
        if(use_Z) {
          Z_shuffle <- Z[shuffle,,drop=FALSE]
          Z_train[[i]] <- Z_shuffle[train_ind,,drop=FALSE]
          Z_val[[i]] <- Z_shuffle[test_ind,,drop=FALSE]
        }
      }

      if(scoring_rule == "auc") {
        train_order <- order(y_train[[i]])
        X_train[[i]] <- X_train[[i]][train_order,,drop=FALSE]
        y_train[[i]] <- y_train[[i]][train_order]
        val_order <- order(y_val[[i]])
        X_val[[i]] <- X_val[[i]][val_order,,drop=FALSE]
        y_val[[i]] <- y_val[[i]][val_order]
        if(use_Z) {
          Z_train[[i]] <- Z_train[[i]][train_order,,drop=FALSE]
          Z_val[[i]] <- Z_val[[i]][val_order,,drop=FALSE]
        }
      }
    }
    use_validation <- TRUE
  } else {
    X_train[[1]] <- X
    y_train[[1]] <- y
    X_val[[1]] <- X
    y_val[[1]] <- y
    if(use_Z) {
      Z_train[[1]] <- Z
      Z_val[[1]] <- Z
    }

    if(scoring_rule == "auc") {
      train_order <- order(y)
      X_train[[1]] <- X[train_order,,drop=FALSE]
      y_train[[1]] <- y[train_order]
      X_val[[1]] <- X[train_order,,drop=FALSE]
      y_val[[1]] <- y[train_order]
      if(use_Z) {
        Z_train[[1]] <- Z[train_order,,drop=FALSE]
        Z_val[[1]] <- Z[train_order,,drop=FALSE]
      }
    }

    val_reps <- 1
    use_validation <- FALSE
  }

  disj <- matrix(nrow = max_conj, ncol = max_vars)
  mode(disj) <- "integer"

  evaluated_models <- list()

  if (allow_conj_removal) {
    min_score <- 1e35
  } else {
    disj[,1] <- sample(p, max_conj, replace = FALSE)
    min_score <- buildModel(X_train, y_train, Z_train, Z_val, disj, list(), tree_control, scoring_rule, gamma,
                            X_val, y_val, use_validation, y_bin)$new_score
  }

  min_conj <- disj
  score <- min_score

  score_rule <- getScoreRule(scoring_rule)

  max_vars <- as.integer(max_vars)
  max_conj <- as.integer(max_conj)
  conjsize <- as.integer(conjsize)

  if(search_algo == "greedy") {
    eval <- .Call(greedySearch_, X_train, y_train, max_vars, max_conj, Z_train, Z_val, disj, score, randomize_greedy, greedy_mod, greedy_rem,
                  tree_control$nodesize, tree_control$split_criterion, tree_control$alpha, tree_control$cp, tree_control$smoothing, tree_control$mtry, tree_control$covariable_search, score_rule, gamma,
                  X_val, y_val, use_validation, y_bin, allow_conj_removal, conjsize, X)
    min_conj <- eval[[1]]
    min_score <- eval[[2]]
    total_iter <- eval[[3]]
    move_counts <- eval[[4]]
  } else if (search_algo == "sa") {

    cs <- cooling_schedule

    # L <- max_conj + p + max_conj * (2 * (p - 1) + max_vars + max_vars * (2 * (p - 1) + 1))
    # (Van Laarhoven and Aarts)

    # L <- 1000
    # L_print_length <- nchar(toString(L))

    if(cs$auto_start_temp) {
      # Configuring start temperature
      start_temp_scores <- vector()
      disj_buffer <- disj
      score_buffer <- Inf
      for(i in 1:cs$start_temp_steps) {
        eval <- simulatedAnnealingStep(X_train, y_train, max_vars, max_conj, Z_train, Z_val, disj_buffer, Inf, 0, evaluated_models, tree_control, scoring_rule, gamma, X_val, y_val, use_validation, y_bin, allow_conj_removal, conjsize, X)
        disj_buffer <- eval$disj

        if(eval$score > score_buffer) {
        #if(i > 1) {
          start_temp_scores <- c(start_temp_scores, eval$score - score_buffer)
        }
        score_buffer <- eval$score
      }
      if(!cs$acc_type2)
        t <- -mean(start_temp_scores)/log(cs$start_acc_ratio)
      else
        t <- quantile(start_temp_scores, cs$start_acc_ratio, names = FALSE)
      if(is.na(t))
        t <- cs$real_start_temp
    } else {
      t <- cs$real_start_temp
    }

    current_scores <- vector()
    current_acc <- vector()

    frozen <- 0

    eval <- .Call(simulatedAnnealing_, X_train, y_train, max_vars, max_conj, Z_train, Z_val, disj, t, score,
                  tree_control$nodesize, tree_control$split_criterion, tree_control$alpha, tree_control$cp, tree_control$smoothing, tree_control$mtry, tree_control$covariable_search, score_rule, gamma,
                  X_val, y_val, use_validation, y_bin, allow_conj_removal, conjsize, X, cs)
    min_conj <- eval[[1]]
    min_score <- eval[[2]]
    total_iter <- eval[[3]]
    prevented_evals <- eval[[4]]
  } else if (search_algo == "genetic") {
    eval <- .Call(geneticProgramming_, X_train, y_train, max_vars, max_conj,
                  as.integer(max_gen), as.numeric(gp_sigma), as.integer(gp_fs_interval),
                  Z_train, Z_val, tree_control$nodesize, tree_control$split_criterion, tree_control$alpha, tree_control$cp,
                  tree_control$smoothing, tree_control$mtry, tree_control$covariable_search, score_rule, gamma, X_val, y_val, use_validation, y_bin,
                  allow_conj_removal, conjsize, X)

    min_conj <- eval[[1]]
    min_score <- eval[[2]]
    total_iter <- eval[[3]]
    prevented_evals <- eval[[4]]
    best_score <- eval[[5]]
  }

  model <- list(disj = min_conj, score = min_score, evaluated_models = evaluated_models,
                X_train = X_train, y_train = y_train, Z_train = Z_train, Z_val = Z_val,
                tree_control = tree_control, scoring_rule = scoring_rule, gamma = gamma, val_method = val_method,
                X_val = X_val, y_val = y_val, use_validation = use_validation, y_bin = y_bin, X = X, y = y, Z = Z,
                total_iter = total_iter, prevented_evals = prevented_evals, conjsize = conjsize)
  class(model) <- "logicDT"

  if(search_algo == "genetic") {
    class(model) <- "geneticLogicPET"
    model$best_score <- best_score

    n_ind <- length(model$disj)
    model$ensemble <- vector("list", length = n_ind)
    model$real_disj <- vector("list", length = n_ind)

    for(i in 1:n_ind) {
      model2 <- model
      model2$disj <- min_conj[[i]]
      model2$score <- min_score[i]
      if (simplify %in% c("vars", "conj")) {
        model2 <- simplifyDisjunctions(model2)
      }
      if (simplify == "vars") {
        model2 <- simplifyConjunctions(model2)
      }
      model$disj[[i]] <- model2$disj

      model$disj[[i]] <- simplifyMatrix(model$disj[[i]])
      model$disj[[i]] <- dontNegateSinglePredictors(model$disj[[i]])
      model$real_disj[[i]] <- translateLogicPET(model$disj[[i]], X)

      ensemble <- fitPETs(X_train, y_train, X_val, y_val, Z_train, Z_val, use_validation, y_bin,
                          tree_control$nodesize, tree_control$split_criterion, tree_control$alpha, tree_control$cp, tree_control$smoothing, tree_control$mtry, tree_control$covariable_final,
                          model$disj[[i]], sum(rowSums(!is.na(model$disj[[i]])) > 0), score_rule, gamma, TRUE)
      model$ensemble[[i]] <- ensemble
    }
    return(model)
  }

  if (simplify %in% c("vars", "conj")) {
    model <- simplifyDisjunctions(model)
  }
  if (simplify == "vars") {
    model <- simplifyConjunctions(model)
  }

  model$disj <- simplifyMatrix(model$disj)
  model$disj <- dontNegateSinglePredictors(model$disj)
  model$real_disj <- translateLogicPET(model$disj, X)

  X2 <- getDesignMatrix(X, model$disj)
  pet <- fitPET(X2, y, Z, tree_control$nodesize, tree_control$split_criterion, tree_control$alpha, tree_control$cp, tree_control$smoothing, tree_control$mtry, tree_control$covariable_final)
  ensemble <- fitPETs(X_train, y_train, X_val, y_val, Z_train, Z_val, use_validation, y_bin,
                      tree_control$nodesize, tree_control$split_criterion, tree_control$alpha, tree_control$cp, tree_control$smoothing, tree_control$mtry, tree_control$covariable_final,
                      model$disj, sum(rowSums(!is.na(model$disj)) > 0), score_rule, gamma, TRUE)
  model$pet <- pet
  model$ensemble <- ensemble

  if(search_algo == "greedy") {
    model$move_counts <- move_counts
  }

  return(model)
}

#' @param formula An object of type \code{formula} describing the
#'   model to be fitted.
#' @param data A data frame containing the data for the corresponding
#'   \code{formula} object. Must also contain quantitative covariables
#'   if they should be included as well.
#' @param ... Arguments passed to \code{\link{logicDT.default}}
#' @rdname logicDT
#' @name logicDT
#' @method logicDT formula
#' @importFrom stats model.frame model.response model.matrix
#' @export
logicDT.formula <- function(formula, data, ...) {
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  predictors <- model.matrix(formula, data = mf)[,-1]
  quant.preds <- apply(predictors, 2, function(x) any(!(x %in% c(0, 1))))
  X <- predictors[,!quant.preds,drop=FALSE]
  Z <- predictors[,quant.preds,drop=FALSE]
  if(ncol(Z) == 0) Z <- NULL
  logicDT(X, y, Z = Z, ...)
}

greedyStep <- function(X_train, y_train, max_vars, max_conj, disj, min_score, evaluated_models, tree_control, scoring_rule, gamma,
                       X_val, y_val, use_validation, allow_conj_removal, conjsize, X) {
  p <- ncol(X_train[[1]])
  min_conj <- disj

  n_conj <- sum(rowSums(!is.na(disj)) > 0)
  n_vars <- rowSums(!is.na(disj[1:n_conj,,drop=FALSE]))
  n_vars_total <- sum(n_vars)

  # Add
  if (n_conj < max_conj & n_vars_total < max_vars) {
    for(j in 1:p) {
      disj2 <- disj
      disj2[n_conj + 1, 1] <- j
      score_and_models <- buildModel(X_train, y_train, disj2, evaluated_models, tree_control, scoring_rule, gamma, X_val, y_val, use_validation)
      if(score_and_models$new_score < min_score) {
        min_score <- score_and_models$new_score
        min_conj <- disj2
      }
    }
  }

  # Remove
  if (n_conj > 1 & allow_conj_removal) {
    for(j in 1:n_conj) {
      disj2 <- disj
      disj2[j,] <- disj2[n_conj,]
      disj2[n_conj,] <- NA
      score_and_models <- buildModel(X_train, y_train, disj2, evaluated_models, tree_control, scoring_rule, gamma, X_val, y_val, use_validation)
      if(score_and_models$new_score < min_score) {
        min_score <- score_and_models$new_score
        min_conj <- disj2
      }
    }
  }

  # Modify
  for(i in 1:n_conj) {
    n_vars_here <- n_vars[i]
    current_conj <- disj[i,]
    unused_vars <- setdiff(1:p, abs(current_conj))
    unused_vars <- c(unused_vars, -unused_vars)

    # Add
    if (n_vars_total < max_vars) {
      for(j in unused_vars) {
        disj2 <- disj
        disj2[i, n_vars_here + 1] <- j

        start2 <- Sys.time()
        conjsum <- sum(getDesignMatrix(X, disj2[i,,drop=FALSE]))
        if (conjsum < conjsize | conjsum > nrow(X) - conjsize) {
          next
        }

        score_and_models <- buildModel(X_train, y_train, disj2, evaluated_models, tree_control, scoring_rule, gamma, X_val, y_val, use_validation)
        if(score_and_models$new_score < min_score) {
          min_score <- score_and_models$new_score
          min_conj <- disj2
        }
      }
    }

    # Remove
    if (n_vars_here > 1) {
      for(j in 1:n_vars_here) {
        disj2 <- disj
        disj2[i, j] <- disj2[i, n_vars_here]
        disj2[i, n_vars_here] <- NA
        score_and_models <- buildModel(X_train, y_train, disj2, evaluated_models, tree_control, scoring_rule, gamma, X_val, y_val, use_validation)
        if(score_and_models$new_score < min_score) {
          min_score <- score_and_models$new_score
          min_conj <- disj2
        }
      }
    }

    # Modify
    for(j in 1:n_vars_here) {
      available_vars <- c(unused_vars, -current_conj[j])
      for(k in available_vars) {
        disj2 <- disj
        disj2[i, j] <- k

        start2 <- Sys.time()
        conjsum <- sum(getDesignMatrix(X, disj2[i,,drop=FALSE]))
        if (conjsum < conjsize | conjsum > nrow(X) - conjsize) {
          next
        }

        score_and_models <- buildModel(X_train, y_train, disj2, evaluated_models, tree_control, scoring_rule, gamma, X_val, y_val, use_validation)
        if(score_and_models$new_score < min_score) {
          min_score <- score_and_models$new_score
          min_conj <- disj2
        }
      }
    }
  }

  ret <- list(score = min_score, disj = min_conj)
  return(ret)
}

simulatedAnnealingStep <- function(X_train, y_train, max_vars, max_conj, Z_train, Z_val, disj, t, score, evaluated_models, tree_control, scoring_rule, gamma,
                                   X_val, y_val, use_validation, y_bin, allow_conj_removal, conjsize, X) {
  p <- ncol(X_train[[1]])
  disj2 <- disj

  main_moves <- c("Modify", "Add", "Remove")

  n_conj <- sum(rowSums(!is.na(disj2)) > 0)
  n_vars <- rowSums(!is.na(disj2[1:n_conj,,drop=FALSE]))
  n_vars_total <- sum(n_vars)

  # Uniform distribution over the neighbours of a specific state
  # as suggested by Van Laarhoven and Aarts

  poss_add_moves <- 2 * p * (n_conj < max_conj) * (n_vars_total < max_vars)
  poss_rem_moves <- n_conj * (n_conj > 1) * allow_conj_removal

  poss_mod_moves <- 2 * (p - n_vars) * (n_vars_total < max_vars) * (n_conj > 0) + n_vars * (n_vars > 1) + n_vars * (2 * (p - n_vars) + 1)
  poss_mod_moves_count <- sum(poss_mod_moves)

  main_move_prob <- c(poss_mod_moves_count, poss_add_moves, poss_rem_moves)/(poss_mod_moves_count + poss_add_moves + poss_rem_moves)

  main_move <- sample(main_moves, 1, prob = main_move_prob)

  if (main_move == "Modify") {
    which_conj <- sample(1:n_conj, 1, prob = poss_mod_moves/poss_mod_moves_count)

    current_conj <- disj2[which_conj,]

    n_vars_here <- n_vars[which_conj]

    poss_mod_add_moves <- 2 * (p - n_vars_here) * (n_vars_total < max_vars)
    poss_mod_rem_moves <- n_vars_here * (n_vars_here > 1)
    poss_mod_mod_moves <- n_vars_here * (2 * (p - n_vars_here) + 1)

    mod_move_prob <- c(poss_mod_mod_moves, poss_mod_add_moves, poss_mod_rem_moves)/(poss_mod_mod_moves + poss_mod_add_moves + poss_mod_rem_moves)

    mod_move <- sample(main_moves, 1, prob = mod_move_prob)
    unused_vars <- setdiff(1:p, abs(current_conj))

    if (mod_move == "Modify") {
      which_var <- sample(1:n_vars_here, 1)
      available_vars <- c(unused_vars, -unused_vars, -current_conj[which_var])

      disj2[which_conj, which_var] <- sample(available_vars, 1)
    } else if (mod_move == "Add") {
      available_vars <- c(unused_vars, -unused_vars)

      disj2[which_conj, n_vars_here + 1] <- sample(available_vars, 1)
    } else if (mod_move == "Remove") {
      which_var <- sample(1:n_vars_here, 1)

      disj2[which_conj, which_var] <- disj2[which_conj, n_vars_here]
      disj2[which_conj, n_vars_here] <- NA
    }

    if (mod_move %in% c("Modify", "Add")) {
      start2 <- Sys.time()
      conjsum <- sum(getDesignMatrix(X, disj2[which_conj,,drop=FALSE]))
      if (conjsum < conjsize | conjsum > nrow(X) - conjsize) {
        return(simulatedAnnealingStep(X_train, y_train, max_vars, max_conj, Z_train, Z_val, disj, t, score, evaluated_models, tree_control, scoring_rule, gamma,
                                      X_val, y_val, use_validation, y_bin, allow_conj_removal, conjsize, X))
      }
    }

  } else if (main_move == "Add") {
    # disj2[n_conj + 1, 1] <- sample(1:p, 1)
    disj2[n_conj + 1, 1] <- sample(c(1:p, -(1:p)), 1)
  } else if (main_move == "Remove") {
    which_rem <- sample(1:n_conj, 1)

    disj2[which_rem,] <- disj2[n_conj,]
    disj2[n_conj,] <- NA
  }

  return(evaluateModel(X_train, y_train, Z_train, Z_val, t, disj, disj2, score, evaluated_models, tree_control, scoring_rule, gamma, X_val, y_val, use_validation, y_bin))
}

#' @importFrom stats runif
evaluateModel <- function(X_train, y_train, Z_train, Z_val, t, disj, disj2, old_score, evaluated_models, tree_control, scoring_rule, gamma,
                          X_val, y_val, use_validation, y_bin) {
  score_and_models <- buildModel(X_train, y_train, Z_train, Z_val, disj2, evaluated_models, tree_control, scoring_rule, gamma, X_val, y_val, use_validation, y_bin)
  new_score <- score_and_models$new_score
  evaluated_models <- score_and_models$evaluated_models
  disj2 <- score_and_models$disj2

  acc <- min(exp((old_score - new_score)/t), 1)
  rnd <- runif(1, 0, 1)

  # Not needed here in R, since this function is only called with t = Inf
  # acc <- old_score - new_score
  # rnd <- -t

  if (acc > rnd) {
    return(list(disj = disj2, score = new_score, evaluated_models = evaluated_models, acc = TRUE))
  } else {
    return(list(disj = disj, score = old_score, evaluated_models = evaluated_models, acc = FALSE))
  }
}

buildModel <- function(X_train, y_train, Z_train, Z_val, disj2, evaluated_models, tree_control, scoring_rule, gamma,
                       X_val, y_val, use_validation, y_bin, submodel = "tree") {
  score_rule <- getScoreRule(scoring_rule)
  new_score <- fitPETs(X_train, y_train, X_val, y_val, Z_train, Z_val, use_validation, y_bin,
                       tree_control$nodesize, tree_control$split_criterion, tree_control$alpha, tree_control$cp, tree_control$smoothing, tree_control$mtry, tree_control$covariable_search,
                       disj2, sum(rowSums(!is.na(disj2)) > 0), score_rule, gamma)

  return(list(new_score = new_score, disj2 = disj2, evaluated_models = evaluated_models))
}

sortMatrix <- function(m) {
  return(m[do.call(order, split(m, col(m))),,drop=FALSE])
}

sortModel <- function(disj) {
  sorted_conj <- apply(disj, 1, sort, na.last = TRUE)
  if (ncol(disj) == 1) {
    sorted_conj <- matrix(sorted_conj, nrow = 1)
  }
  return(sortMatrix(t(sorted_conj)))
}

getPredictorNames <- function(real_disj, sort_conj = FALSE) {
  n_conj <- sum(rowSums(!is.na(real_disj)) > 0)
  disj2 <- split(real_disj[1:n_conj,,drop=FALSE], 1:n_conj)
  if(sort_conj)
    disj2 <- lapply(disj2, sort)
  disj2 <- lapply(disj2, function(x) x[!is.na(x)])
  disj2 <- lapply(disj2, paste, collapse="^")
  disj2 <- unlist(disj2, use.names = FALSE)
  return(disj2)
}

encodeModel <- function(disj) {
  disj2 <- apply(disj, 1, paste, collapse="^")
  return(paste(disj2, collapse = "v"))
}

# Without NA's
decodeModel <- function(code) {
  n_conj <- length(code)
  # n_vars <- lengths(regmatches(code, gregexpr("\\^", code))) + 1
  conj <- strsplit(code, "\\^")
  n_vars <- lengths(conj)
  max_len <- max(n_vars)
  conj <- lapply(conj, function (x) {c(as.integer(x), rep(NA, max_len - length(x)))})
  disj <- do.call(rbind, conj)
  return(sortMatrix(disj))
}

#' Design matrix for the set of conjunctions
#'
#' Transform the original predictor matrix X into the conjunction design matrix
#' which contains for each conjunction a corresponding column.
#'
#' @param X The original (binary) predictor matrix. This has to be of type
#'   integer.
#' @param disj The conjunction matrix which can, e.g., be extracted from a
#'   fitted \code{logicDT} model via $disj.
#' @return The transformed design matrix.
#'
#' @export
getDesignMatrix <- function(X, disj) {
  dm <- .Call(getDesignMatrix_, X, disj, sum(rowSums(!is.na(disj)) > 0))
  return(dm)
}

simplifyDisjunctions <- function(model) {
  X_train <- model$X_train
  y_train <- model$y_train
  Z_train <- model$Z_train
  Z_val <- model$Z_val
  disj <- model$disj
  evaluated_models <- model$evaluated_models
  tree_control <- model$tree_control
  scoring_rule <- model$scoring_rule
  gamma <- model$gamma
  X_val <- model$X_val
  y_val <- model$y_val
  use_validation <- model$use_validation
  y_bin <- model$y_bin

  score <- model$score
  n_conj <- sum(rowSums(!is.na(disj)) > 0)

  if(n_conj < 2)
    return(model)

  for(j in 1:n_conj) {
    disj2 <- disj
    disj2[j,] <- disj2[n_conj,]
    disj2[n_conj,] <- NA
    score_and_models <- buildModel(X_train, y_train, Z_train, Z_val, disj2, evaluated_models, tree_control, scoring_rule, gamma,
                                   X_val, y_val, use_validation, y_bin)
    if(score_and_models$new_score <= score) {
      model$disj <- disj2
      return(simplifyDisjunctions(model))
    }
  }
  return(model)
}

simplifyConjunctions <- function(model) {
  disj <- model$disj
  score <- model$score

  n_conj <- sum(rowSums(!is.na(disj)) > 0)

  n_vars <- rowSums(!is.na(disj[1:n_conj,,drop=FALSE]))
  X_train <- model$X_train
  y_train <- model$y_train
  Z_train <- model$Z_train
  Z_val <- model$Z_val
  evaluated_models <- model$evaluated_models
  tree_control <- model$tree_control
  scoring_rule <- model$scoring_rule
  gamma <- model$gamma
  X_val <- model$X_val
  y_val <- model$y_val
  use_validation <- model$use_validation
  y_bin <- model$y_bin

  # Hm
  disj <- sortModel(disj)

  for (i in 1:n_conj) {
    current_conj <- disj[i,]
    current_n_vars <- n_vars[i]
    min_n_vars <- current_n_vars

    results <- simplifyConjunction(disj, i, current_n_vars, score, list(),
                                   X_train, y_train, Z_train, Z_val, evaluated_models, tree_control, scoring_rule, gamma,
                                   X_val, y_val, use_validation, y_bin)
    evaluated_models <- results$evaluated_models
    variations <- results$variations

    if (length(variations) > 0) {
      variation <- names(sort(unlist(variations))[1])

      conj <- strsplit(variation, "\\^")[[1]]
      conj <- c(as.integer(conj), rep(NA, ncol(disj) - length(conj))) # HERE
      disj[i,] <- conj
    }
  }

  disj <- sortModel(disj)
  model$disj <- disj
  model$evaluated_models <- evaluated_models

  return(model)
}

simplifyConjunction <- function(disj, disj_ind, current_n_vars, score, variations,
                                X_train, y_train, Z_train, Z_val, evaluated_models, tree_control, scoring_rule, gamma,
                                X_val, y_val, use_validation, y_bin) {
  if (current_n_vars < 2) {
    return(list(variations = variations, evaluated_models = evaluated_models))
  }

  for (j in 1:current_n_vars) {
    disj2 <- disj
    disj2[disj_ind, j] <- disj2[disj_ind, current_n_vars]
    disj2[disj_ind, current_n_vars] <- NA

    score_and_models <- buildModel(X_train, y_train, Z_train, Z_val, disj2, evaluated_models, tree_control, scoring_rule, gamma,
                                   X_val, y_val, use_validation, y_bin)
    evaluated_models <- score_and_models$evaluated_models
    new_score <- score_and_models$new_score

    if (new_score <= score) {
      variations[[paste(disj2[disj_ind, 1:(current_n_vars - 1)], collapse = "^")]] <- current_n_vars - 1
      rec <- simplifyConjunction(disj2, disj_ind, current_n_vars - 1, score, variations,
                                 X_train, y_train, Z_train, Z_val, evaluated_models, tree_control, scoring_rule, gamma,
                                 X_val, y_val, use_validation, y_bin)
      variations <- rec$variations
      evaluated_models <- rec$evaluated_models
    }
  }
  return(list(variations = variations, evaluated_models = evaluated_models))
}

dontNegateSinglePredictors <- function(disj) {
  n_vars <- rowSums(!is.na(disj[,,drop=FALSE]))
  disj[n_vars == 1, 1] <- abs(disj[n_vars == 1, 1])
  disj
}

#' Refit the logic decision trees
#'
#' Newly fit the decision trees in the \code{logicDT} model using
#' the supplied tree control parameters.
#' This is especially useful if, e.g., the model was initially trained
#' without utilizing a continuous covariable or fitting linear models and
#' now 4pL model shall be fitted.
#'
#' @param model A fitted \code{logicDT} model
#' @param tree_control Tree control parameters. This object should be
#'   constructed using the function \code{\link{tree.control}}.
#'   Alternatively, the old \code{tree_control} from \code{model} can be
#'   modified and specified here.
#' @return The \code{logicDT} model with newly fitted trees
#'
#' @export
refitTrees <- function(model, tree_control) {
  model$tree_control <- tree_control

  score_rule <- getScoreRule(model$scoring_rule)
  X2 <- getDesignMatrix(model$X, model$disj)
  pet <- fitPET(X2, model$y, model$Z, model$tree_control$nodesize, model$tree_control$split_criterion,
                model$tree_control$alpha, model$tree_control$cp, model$tree_control$smoothing,
                model$tree_control$mtry, model$tree_control$covariable_final)
  ensemble <- fitPETs(model$X_train, model$y_train, model$X_val, model$y_val, model$Z_train, model$Z_val,
                      model$use_validation, model$y_bin, model$tree_control$nodesize,
                      model$tree_control$split_criterion, model$tree_control$alpha, model$tree_control$cp,
                      model$tree_control$smoothing, model$tree_control$mtry, model$tree_control$covariable_final,
                      model$disj, sum(rowSums(!is.na(model$disj)) > 0), score_rule, model$gamma, TRUE)
  model$pet <- pet
  model$ensemble <- ensemble
  return(model)
}

simplifyMatrix <- function(disj) {
  n_conj <- sum(rowSums(!is.na(disj)) > 0)
  n_vars <- rowSums(!is.na(disj[1:n_conj,,drop=FALSE]))
  return(disj[1:n_conj, 1:max(n_vars), drop=FALSE])
}

calcAUCFast <- function(probs, y, sorted = FALSE) {
  return(.Call(calcAUC_, probs, y, sorted))
}

getScoreRule <- function(scoring_rule) {
  if(scoring_rule == "deviance")
    score_rule <- as.integer(0)
  else if(scoring_rule == "brier")
    score_rule <- as.integer(1)
  else if(scoring_rule == "nce")
    score_rule <- as.integer(5)
  else
    score_rule <- as.integer(4)
  score_rule
}

#' Prediction for logicDT models
#'
#' Supply new input data for predicting the outcome with a fitted
#' logicDT model.
#'
#' @param object Fitted logicDT model. Usually a product of a call
#'   to \code{\link{logicDT}}.
#' @param X Matrix or data frame of binary input data. This
#'   object should correspond to the binary matrix for fitting
#'   the model.
#' @param Z Optional quantitative covariables supplied as a
#'   matrix or data frame. Only used (and required) if the
#'   model was fitted using them.
#' @param type Prediction type. This can either be "prob" for
#'   probability estimates or "class" for
#'   classification in binary responses. Ignored for regression.
#' @param ensemble If the model was fitted using the inner
#'   validation approach, shall the prediction be constructed
#'   using the final validated ensemble (TRUE) or using the
#'   single final tree (FALSE)?
#' @param leaves If four parameter logistic models were fitted
#'   for each leaf, shall they be used for the prediction
#'   ("4pl") or shall the constant leaf means be used
#'   ("constant")?
#' @param models Which models of logicDT model fitted with
#'   genetic programming shall be used for prediction?
#'   "best" leads to the single best model in the final
#'   generation, "all" uses the average over the final
#'   generation and "n_models" uses the n_models best models.
#' @param n_models How many models shall be used if
#'   models = "n_models" and genetic programming was employed?
#' @param ... Parameters supplied to \code{\link{predict.logicDT}}
#' @return A numeric vector of predictions. For binary outcomes,
#'   this is a vector with estimates for \eqn{P(Y=1 \mid X = x)}.
#'
#' @export
predict.logicDT <- function(object, X, Z = NULL, type = "prob", ensemble = FALSE, leaves = "4pl", ...) {
  X <- as.matrix(X)
  mode(X) <- "integer"
  X2 <- getDesignMatrix(X, object$disj)
  if(!is.null(Z)) {
    Z <- as.matrix(Z)
    mode(Z) <- "double"
  } else {
    if(!is.null(object$Z_train)) {
      stop("The covariables are missing!")
    }
  }
  if (ensemble)
    return(predict_ensemble(object$ensemble, X2, Z, type, leaves))
  else
    return(predict_pet(object$pet, X2, Z, type, leaves))
}

#' @rdname predict.logicDT
#' @export
predict.geneticLogicPET <- function(object, X, Z = NULL, models = "best", n_models = 10, ensemble = NULL, leaves = "4pl", ...) {
  X <- as.matrix(X)
  mode(X) <- "integer"
  if(!is.null(Z)) {
    Z <- as.matrix(Z)
    mode(Z) <- "double"
  } else {
    if(!is.null(object$Z_train)) {
      stop("The covariables are missing!")
    }
  }

  if(models == "best")
    type <- 0
  else if(models == "all")
    type <- 1
  else if(models == "n_models")
    type <- 2

  if(leaves == "4pl")
    leaves <- 1
  else
    leaves <- 0

  return(.Call(predictGP_, object, X, Z, type, n_models, leaves))
}

translateLogicPET <- function(disj, X) {
  translated <- t(apply(disj, 1, function(row) ifelse(is.na(row), NA, paste(ifelse(row < 0, "-", ""), colnames(X)[abs(row)], sep=""))))
  if(ncol(disj) == 1)
    translated <- t(translated)
  return(translated)
}

#' Fitting 4pL models
#'
#' Method for fitting four parameter logistic models.
#' In the fashion of this package, only binary and quantitative
#' outcomes are supported.
#'
#' 4pL models are non-linear regression models of the shape
#' \deqn{Y = f(x, b, c, d, e) + \varepsilon =
#'   c + \frac{d-c}{1+\exp(b \cdot (x-e))} + \varepsilon}
#' with \eqn{\varepsilon} being a random error term.
#'
#' @param y Response vector. 0-1 coding for binary outcomes,
#'   otherwise conventional regression is performed.
#' @param Z Numeric vector of (univariate) input samples.
#' @return An object of class "4pl" which contains a numeric
#'   vector of the fitted parameters b, c, d, and e.
#'
#' @export
fit4plModel <- function(y, Z) {
  if(any(!(y %in% c(0, 1)))) {
    y <- as.numeric(y)
  } else {
    y <- as.integer(y)
  }
  .Call(fit4plModel_, y, Z)
}

#' Prediction for 4pL models
#'
#' Use new input data and a fitted four parameter logistic
#' model to predict corresponding outcomes.
#'
#' @param object Fitted 4pl model
#' @param Z Numeric vector of new input samples
#' @param ... Ignored additional parameters
#' @return A numeric vector of predictions. For binary outcomes,
#'   this is a vector with estimates for
#'   \eqn{P(Y=1 \mid X = x)}.
#'
#' @export
predict.4pl <- function(object, Z, ...) {
  bcde <- object[[1]]
  y_bin <- object[[2]]
  ret <- bcde[2] + (bcde[3]-bcde[2])/(1+exp(bcde[1]*(Z-bcde[4])))
  if(y_bin) {
    ret[ret > 1] <- 1
    ret[ret < 0] <- 0
  }
  return(ret)
}

#' Fitting linear models
#'
#' Method for fitting linear models.
#' In the fashion of this package, only binary and quantitative
#' outcomes are supported.
#'
#' For binary outcomes, predictions are cut at 0 or 1 for generating
#' proper probability estimates.
#'
#' @param y Response vector. 0-1 coding for binary outcomes,
#'   otherwise conventional regression is performed.
#' @param Z Numeric vector of (univariate) input samples.
#' @return An object of class "linear" which contains a numeric
#'   vector of the fitted parameters b and c.
#'
#' @export
fitLinearModel <- function(y, Z) {
  if(any(!(y %in% c(0, 1)))) {
    y <- as.numeric(y)
  } else {
    y <- as.integer(y)
  }
  .Call(fitLinearModel_, y, Z)
}

#' Prediction for linear models
#'
#' Use new input data and a fitted linear
#' model to predict corresponding outcomes.
#'
#' For binary outcomes, predictions are cut at 0 or 1 for generating
#' proper probability estimates.
#'
#' @param object Fitted linear model
#' @param Z Numeric vector of new input samples
#' @param ... Ignored additional parameters
#' @return A numeric vector of predictions. For binary outcomes,
#'   this is a vector with estimates for
#'   \eqn{P(Y=1 \mid X = x)}.
#'
#' @export
predict.linear <- function(object, Z, ...) {
  bcde <- object[[1]]
  y_bin <- object[[2]]
  ret <- bcde[1] + bcde[2] * Z
  if(y_bin) {
    # Logistic transformation due to LDA model:
    ret <- 1/(1 + exp(-ret))
  }
  return(ret)
}

.onUnload <- function (libpath) {
  library.dynam.unload("logicDT", libpath)
}
