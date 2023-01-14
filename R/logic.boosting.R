#' @export
#' @rawNamespace S3method(logicDT,boosting)
logicDT.boosting <- function(X, ...) UseMethod("logicDT.boosting")

#' Fitting boosted logicDT models
#'
#' Function for fitting gradient boosted logicDT models.
#'
#' Details on single logicDT models can be found in \code{\link{logicDT}}.
#'
#' @param X Matrix or data frame of binary predictors coded as 0 or 1.
#' @param y Response vector. 0-1 coding for binary responses.
#'  Otherwise, a regression task is assumed.
#' @param Z Optional matrix or data frame of quantitative/continuous
#'   covariables. Multiple covariables allowed for splitting the trees.
#'   If four parameter logistic models shall be fitted in the leaves,
#'   only the first given covariable is used.
#' @param boosting.iter Number of boosting iterations
#' @param learning.rate Learning rate for boosted models.
#'   Values between 0.001 and 0.1 are recommended.
#' @param subsample.frac Subsample fraction for each
#'   boosting iteration. E.g., 0.5 means that are random draw
#'   of 50% of the original training data observations
#'   is used in each iteration.
#' @param replace Should the random draws with subsample.frac
#'   in boosted models be performed with or without
#'   replacement? TRUE or FALSE
#' @param line.search Type of line search for gradient boosting.
#'   "min" performs a real minimization while "binary" performs
#'   a loose binary search for a boosting coefficient that
#'   just reduces the score.
#' @param ... Arguments passed to \code{logicDT}
#' @return An object of class \code{logic.boosted}. This is a list
#'   containing
#'   \item{\code{models}}{A list of fitted logicDT models}
#'   \item{\code{rho}}{A vector of boosting coefficient corresponding
#'   to each model}
#'   \item{\code{initialModel}}{Initial model which is usually the
#'   observed mean}
#'   \item{\code{...}}{Supplied parameters of the functional call
#'     to \code{\link{logicDT.boosting}}.}
#' @references
#' \itemize{
#'   \item Lau, M., Schikowski, T. & Schwender, H. (2021).
#'   logicDT: A Procedure for Identifying Response-Associated
#'   Interactions Between Binary Predictors. To be submitted.
#'   \item Friedman, J. H. (2001).
#'   Greedy Function Approximation: A Gradient Boosting Machine.
#'   The Annals of Statistics, 29(5), 1189–1232.
#'   \doi{https://doi.org/10.1214/aos/1013203451}
#' }
#'
#' @name logicDT.boosting
#' @method logicDT.boosting default
#' @importFrom utils flush.console
#' @export
logicDT.boosting.default <- function(X, y, Z = NULL, boosting.iter = 500, learning.rate = 0.01,
                                     subsample.frac = 1, replace = TRUE, line.search = "min", ...) {
  y_bin <- setequal(unique(y), c(0, 1))

  total_iter_boosting <- 0
  ret <- list()
  ret$models <- list()
  ret$bags <- list()
  ret$X <- X
  ret$y <- y
  ret$Z <- Z
  ret$learning.rate <- learning.rate
  ret$y_bin <- y_bin

  N <- nrow(X)
  Nadj <- floor(subsample.frac * N)
  rhoVector <- vector()
  priorProb <- mean(y)

  if(y_bin) {
    initialModel <- 0.5 * log(priorProb/(1-priorProb))
  } else {
    initialModel <- priorProb
  }

  currentEstimates <- rep(initialModel, N)
  evaluatedWeakLearners <- array(0, c(boosting.iter, N))

  for(i in 1:boosting.iter) {
    bag <- sample(N, Nadj, replace = replace)
    ret$bags[[i]] <- bag
    if(is.null(Z))
      Z_temp <- NULL
    else
      Z_temp <- Z[bag,,drop=FALSE]

    if(y_bin) {
      probs <- 1/(1+exp(-2 * currentEstimates[bag]))
      gradient <- 2*probs - 2*y[bag]
    } else {
      gradient <- currentEstimates[bag] - y[bag]
    }

    ret$models[[i]] <- logicDT(X[bag,,drop=FALSE], -gradient, Z = Z_temp, ...)

    evaluatedWeakLearners[i,] <- predict(ret$models[[i]], X, Z)

    currentRho <- getRho(y[bag], currentEstimates[bag], evaluatedWeakLearners[i,][bag], line.search, y_bin)

    rhoVector <- c(rhoVector, currentRho)
    currentEstimates <- currentEstimates + learning.rate * currentRho * evaluatedWeakLearners[i,]

    total_iter_boosting <- total_iter_boosting + ret$models[[i]]$total_iter
    cat("\r", sprintf("Iteration %d/%d (%.0f%%)", i, boosting.iter, i/boosting.iter * 100))
    flush.console()
  }
  cat("\n")
  ret$total_iter <- total_iter_boosting
  ret$rho <- rhoVector
  ret$initialModel <- initialModel
  class(ret) <- "logic.boosted"
  return(ret)
}

#' @param formula An object of type \code{formula} describing the
#'   model to be fitted.
#' @param data A data frame containing the data for the corresponding
#'   \code{formula} object. Must also contain quantitative covariables
#'   if they should be included as well.
#' @rdname logicDT.boosting
#' @name logicDT.boosting
#' @method logicDT.boosting formula
#' @importFrom stats model.frame model.response model.matrix
#' @export
logicDT.boosting.formula <- function(formula, data, ...) {
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  predictors <- model.matrix(formula, data = mf)[,-1]
  quant.preds <- apply(predictors, 2, function(x) any(!(x %in% c(0, 1))))
  X <- predictors[,!quant.preds,drop=FALSE]
  Z <- predictors[,quant.preds,drop=FALSE]
  if(ncol(Z) == 0) Z <- NULL
  logicDT.boosting(X, y, Z = Z, ...)
}

#' @importFrom stats optim
getRho <- function(y, oldEstimates, evaluatedWeakLearners, line.search, y_bin) {
  # For line.search = "bin":
  # Largest rho which leads to a smaller score via binary search

  if(y_bin) {
    calcScore <- function(x) {
      return(sum(log(1+exp(2*(oldEstimates + x * evaluatedWeakLearners))) - 2 * y * (oldEstimates + x*evaluatedWeakLearners)))
    }
  } else {
    calcScore <- function(x) {
      return(sum((y - oldEstimates - x * evaluatedWeakLearners)^2)/2)
    }
  }

  if(line.search == "min") {
    if(y_bin) {
      return(optim(1, calcScore, method = "Brent", lower = 0, upper = 8)$par[1])
    } else {
      return(sum((y - oldEstimates) * evaluatedWeakLearners)/sum((evaluatedWeakLearners)^2))
    }
  }

  rho <- 1.0

  oldScore <- calcScore(0)
  newScore <- calcScore(rho)

  if (newScore < oldScore) {
    while ((newScore < oldScore) & (rho < 256)) {
      rho <- 2 * rho
      newScore <- calcScore(rho)
    }
    rho <- rho / 2
  } else {
    while (newScore >= oldScore) {
      rho <- rho / 2
      newScore <- calcScore(rho)
    }
  }
  # Ideal rho lies in [rho, 2*rho[
  return(recursiveRho(calcScore, rho, 2*rho, 0))
}

recursiveRho <- function(calcScore, lowerRho, upperRho, recDepth) {
  if (recDepth >= 8) {
    return(lowerRho)
  }

  width <- upperRho - lowerRho
  newRho <- lowerRho + width/2
  if (calcScore(newRho) < calcScore(0)) {
    return(recursiveRho(calcScore, newRho, upperRho, recDepth + 1))
  } else {
    return(recursiveRho(calcScore, lowerRho, newRho, recDepth + 1))
  }
}

#' @rdname predict.logicDT
#' @name predict.logicDT
#' @export
predict.logic.boosted <- function(object, X, Z = NULL, type = "prob", ...) {
  preds <- rep(object$initialModel, nrow(X))
  for(i in 1:length(object$models)) {
    preds <- preds + object$learning.rate * object$rho[i] * predict(object$models[[i]], X, Z, ...)
  }
  if(object$y_bin) {
    preds <- 1/(1+exp(-2 * preds))
  }
  if(type == "class") {
    preds <- as.numeric(preds > 0.5)
  }
  return(preds)
}

#' Partial prediction for boosted models
#'
#' Alternative prediction function for \code{logic.boosted} models
#' using up to \code{n.iter} boosting iterations.
#' An array of predictions for every number of boosting iterations
#' up to \code{n.iter} is returned.
#'
#' The main purpose of this function is to retrieve the optimal number of
#' boosting iterations (early stopping) using a validation data set and to
#' restrict future predictions on this number of iterations.
#'
#' @param model Fitted \code{logic.boosted} model
#' @param X Matrix or data frame of binary input data.
#'   This object should correspond to the binary matrix for fitting the model.
#' @param Z Optional quantitative covariables supplied as a matrix or
#'   data frame. Only used (and required) if the model was fitted using them.
#' @param n.iter Maximum number of boosting iterations for prediction
#' @param ... Parameters supplied to \code{\link{predict.logicDT}}
#' @return An array of dimension \code{(N, n.iter)} containing the partial
#'   predictions
#'
#' @export
partial.predict <- function(model, X, Z = NULL, n.iter = 1, ...) {
  N <- nrow(X)
  preds <- array(0.0, dim = c(N, n.iter))
  boosting.iter <- length(model$models)
  n.iter <- min(n.iter, boosting.iter)
  for(i in 1:n.iter) {
    if(i == 1)
      preds[, i] <- model$initialModel + model$learning.rate * model$rho[i] * predict(model$models[[i]], X, Z, ...)
    else
      preds[, i] <- preds[, i-1] + model$learning.rate * model$rho[i] * predict(model$models[[i]], X, Z, ...)
  }
  if(model$y_bin) {
    preds <- 1/(1+exp(-2 * preds))
  }
  return(preds)
}

#' Get the best number of boosting iterations
#'
#' This function can be used to compute the ideal number of boosting iterations
#' for the fitted \code{logic.boosted} model using independent validation data.
#'
#' If the model performance (on the validation data) cannot be increased for
#' \code{consec.iter} consecutive boosting iterations, the last iteration
#' which increased the validation performance induces the ideal number of
#' boosting iterations.
#'
#' @param model Fitted \code{logic.boosted} model
#' @param X Matrix or data frame of binary validation input data.
#'   This object should correspond to the binary matrix for fitting the model.
#' @param y Validation response vector. 0-1 coding for binary outcomes.
#' @param Z Optional quantitative covariables supplied as a matrix or
#'   data frame. Only used (and required) if the model was fitted using them.
#' @param consec.iter Number of consecutive boosting iterations that do not
#'   increase the validation performance for determining the ideal number of
#'   iterations
#' @param scoring_rule Scoring rule computing the validation performance.
#'   This can either be "auc" for the area under the receiver
#'   operating characteristic curve (default for binary reponses),
#'   "deviance" for the deviance, "nce" for the normalized cross entropy
#'   or "brier" for the Brier score.
#'   For regression purposes, the MSE (mean squared error) is
#'   automatically chosen.
#' @return The ideal number of boosting iterations
#'
#' @export
bestBoostingIter <- function(model, X, y, Z = NULL, consec.iter = 5, scoring_rule = "auc") {
  boosting.iter <- length(model$models)
  partial.preds <- partial.predict(model, X, Z = Z, n.iter = boosting.iter)
  best.perf <- -Inf; k <- 0
  for(j in 1:boosting.iter) {
    current.perf <- eval.performance(partial.preds[,j], y, scoring_rule)
    if(current.perf > best.perf) {
      best.perf <- current.perf; k <- 0
    } else {
      k <- k + 1
      if(k == consec.iter) break
    }
  }
  best.iter <- j - k
  return(best.iter)
}

eval.performance <- function(preds, y, scoring_rule) {
  y_bin <- setequal(unique(y), c(0, 1))
  if(y_bin) {
    y <- as.integer(y)
    if(scoring_rule == "deviance")
      -calcDev(preds, y)
    else if(scoring_rule == "brier")
      -calcBrier(preds, y)
    else if(scoring_rule == "nce")
      1-calcNCE(preds, y)
    else
      calcAUC(preds, y)
  } else {
    -calcMSE(preds, y)
  }
}

#' Linear models based on logic terms
#'
#' This function fits a linear or logistic regression model (based on the
#' type of outcome) using the supplied logic terms, e.g., \code{$disj} from
#' a fitted \code{logicDT} model.
#'
#' For creating sparse final models, the LASSO can be used for shrinking
#' unnecessary term coefficients down to zero (\code{type = "lasso"}).
#' If the complexity penalty \code{s} shall be automatically tuned,
#' cross-validation can be employed (\code{type = "cv.lasso"}).
#' However, since other hyperparameters also have to be tuned when fitting
#' a linear boosting model such as the complexity penalty for restricting
#' the number of variables in the terms, manually tuning the LASSO penalty
#' together with the other hyperparameters is recommended.
#' For every hyperparameter setting of the boosting itself, the best
#' corresponding LASSO penalty \code{s} can be identified by, e.g., choosing
#' the \code{s} that minimizes the validation data error.
#' Thus, this hyperparameter does not have to be explicitly tuned via a grid
#' search but is induced by the setting of the other hyperparameters.
#' For finding the ideal value of \code{s} using independent validation data,
#' the function \code{\link{get.ideal.penalty}} can be used.
#'
#' @param X Matrix or data frame of binary input data.
#'   This object should correspond to the binary matrix for fitting the model.
#' @param y Response vector. 0-1 coding for binary outcomes.
#' @param Z Optional quantitative covariables supplied as a matrix or
#'   data frame. Only used (and required) if the model was fitted using them.
#' @param disj Integer matrix of logic terms. As in \code{\link{logicDT}},
#'   each row corresponds to a term/conjunction. Negative values indicate
#'   negations. The absolute values of an entry correspond to the predictor
#'   index in \code{X}.
#' @param Z.interactions Shall interactions with the continuous covariable
#'   \code{Z} be taken into account by including products of the terms with
#'   \code{Z}?
#' @param type Type of linear model to be fitted. Either \code{"standard"}
#'   (without regularization), \code{"lasso"} (LASSO) or \code{"cv.lasso"}
#'   (LASSO with cross-validation for automatically configuring the complexity
#'   penalty).
#' @param s Regularization parameter. Only used if \code{type = "lasso"} is
#'   set.
#' @param ... Additional parameters passed to \code{glmnet} or \code{cv.glmnet}
#'   if the corresponding model type was chosen.
#' @return A \code{linear.logic} model. This is a list containing
#'   the logic terms used as predictors in the model and the fitted \code{glm}
#'   model.
#' @references
#' \itemize{
#'   \item Tibshirani, R. (1996).
#'   Regression Shrinkage and Selection via the Lasso. Journal of the Royal
#'   Statistical Society. Series B (Methodological), 58(1), 267–288.
#'   \doi{https://doi.org/10.1111/j.2517-6161.1996.tb02080.x}
#'   \item Friedman, J., Hastie, T., & Tibshirani, R. (2010).
#'   Regularization Paths for Generalized Linear Models via Coordinate Descent.
#'   Journal of statistical software, 33(1), 1–22.
#'   \doi{https://doi.org/10.18637/jss.v033.i01}
#' }
#'
#' @importFrom stats glm gaussian binomial
#' @importFrom glmnet cv.glmnet glmnet
#' @export
fitLinearLogicModel <- function(X, y, Z = NULL, disj, Z.interactions = TRUE,
                                type = "standard", s = NULL, ...) {
  y_bin <- setequal(unique(y), c(0, 1))
  X <- as.matrix(X)
  mode(X) <- "integer"
  mode(disj) <- "integer"

  if(nrow(disj) == 0) {
    dm <- matrix(nrow = nrow(X), ncol = 0)
    mode(dm) <- "integer"
  } else {
    dm <- getDesignMatrix(X, disj)
  }

  if(!is.null(Z)) {
    Z <- as.numeric(as.matrix(Z))
    if(Z.interactions)
      dm <- cbind(dm, dm * Z, Z = Z)
    else
      dm <- cbind(dm, Z = Z)
    colnames(dm) <- paste("X", 1:ncol(dm), sep="")
  }
  if(type == "standard") {
    if(y_bin) fam <- binomial(link = "logit") else fam <- gaussian()
    dat <- data.frame(y = y, dm)
    lin.mod <- glm(y ~ ., data = dat, family = fam)
  } else if(type == "cv.lasso") {
    if(y_bin) fam <- "binomial" else fam <- "gaussian"
    if(ncol(dm) == 1) dm <- cbind(dm, 1)
    # alpha = 1 (LASSO) is the default setting
    lin.mod <- cv.glmnet(dm, y, family = fam, ...)
    s <- lin.mod$lambda.1se
  } else {
    if(y_bin) fam <- "binomial" else fam <- "gaussian"
    if(ncol(dm) == 1) dm <- cbind(dm, 1)
    # alpha = 1 (LASSO) is the default setting
    lin.mod <- glmnet(dm, y, family = fam, lambda = s, ...)
  }
  ret <- list(disj = disj, lin.mod = lin.mod, Z.interactions = Z.interactions, type = type, s = s)
  class(ret) <- "linear.logic"
  return(ret)
}

#' Prediction for \code{linear.logic} models
#'
#' Use new input data and a fitted \code{linear.logic} model to
#' predict corresponding outcomes.
#'
#' @param object Fitted \code{linear.logic} model
#' @param X Matrix or data frame of binary input data.
#'   This object should correspond to the binary matrix for fitting the model.
#' @param Z Optional quantitative covariables supplied as a matrix or
#'   data frame. Only used (and required) if the model was fitted using them.
#' @param s Regularization parameter. Only used if \code{type = "lasso"} or
#'   \code{type = "cv.lasso"} was set. Only useful if the penalty saved in
#'   \code{object$s} should be overwritten.
#' @param ... Ignored additional parameters
#' @return A numeric vector of predictions. For binary outcomes,
#'   this is a vector with estimates for
#'   \eqn{P(Y=1 \mid X = x)}.
#'
#' @importFrom stats predict.glm
#' @importFrom glmnet predict.glmnet
#' @export
predict.linear.logic <- function(object, X, Z = NULL, s = NULL, ...) {
  X <- as.matrix(X)
  mode(X) <- "integer"

  if(nrow(object$disj) == 0) {
    dm <- matrix(nrow = nrow(X), ncol = 0)
    mode(dm) <- "integer"
  } else {
    dm <- getDesignMatrix(X, object$disj)
  }

  if(!is.null(Z)) {
    Z <- as.numeric(as.matrix(Z))
    if(object$Z.interactions)
      dm <- cbind(dm, dm * Z, Z = Z)
    else
      dm <- cbind(dm, Z = Z)
    colnames(dm) <- paste("X", 1:ncol(dm), sep="")
  }

  if(is.null(s)) s <- object$s

  if(object$type == "standard")
    as.numeric(predict(object$lin.mod, data.frame(dm), type = "response"))
  else if(object$type == "cv.lasso") {
    if(ncol(dm) == 1) dm <- cbind(dm, 1)
    ret <- predict(object$lin.mod, dm, type = "response", s = s)
    if(ncol(ret) == 1) as.numeric(ret) else ret
  } else {
    if(ncol(dm) == 1) dm <- cbind(dm, 1)
    ret <- predict(object$lin.mod, dm, type = "response", s = s)
    if(ncol(ret) == 1) as.numeric(ret) else ret
  }
}

#' Linear models based on boosted models
#'
#' This function uses a fitted \code{logic.boosted} model for fitting
#' a linear or logistic (depending on the type of outcome) regression model.
#'
#' In this procedure, the logic terms are extracted from the individual
#' \code{logicDT} models and the set of unique terms are used as predictors
#' in a regression model. For incorporating a continuous covariable
#' the covariable itself as well as products of the covariable with the
#' extracted logic terms are included as predictors in the regression model.
#'
#' For more details on the possible types of linear models, see
#' \code{\link{fitLinearLogicModel}}.
#'
#' @param model Fitted \code{logic.boosted} model
#' @param n.iter Number of boosting iterations to be used
#' @param type Type of linear model to be fitted. Either \code{"standard"}
#'   (without regularization), \code{"lasso"} (LASSO) or \code{"cv.lasso"}
#'   (LASSO with cross-validation for automatically configuring the complexity
#'   penalty).
#' @param s Regularization parameter. Only used if \code{type = "lasso"} is
#'   set.
#' @param ... Additional parameters passed to \code{glmnet} or \code{cv.glmnet}
#'   if the corresponding model type was chosen.
#' @return A \code{linear.logic} model. This is a list containing
#'   the logic terms used as predictors in the model and the fitted \code{glm}
#'   model.
#'
#' @export
fitLinearBoostingModel <- function(model, n.iter, type = "standard", s = NULL, ...) {
  new.disj <- model$models[[1]]$disj
  max_vars <- ncol(new.disj)
  if(n.iter > 1) {
    for(j in 2:n.iter) {
      current.disj <- model$models[[j]]$disj
      if(ncol(current.disj) > max_vars) {
        new.disj <- cbind(new.disj, matrix(NA_integer_, nrow=nrow(new.disj), ncol=ncol(current.disj)-max_vars))
        max_vars <- ncol(current.disj)
      } else if(ncol(current.disj) < max_vars) {
        current.disj <- cbind(current.disj, matrix(NA_integer_, nrow=nrow(current.disj), ncol=max_vars-ncol(current.disj)))
      }
      new.disj <- rbind(new.disj, current.disj)
    }
  }

  new.disj <- new.disj[!duplicated(new.disj),,drop=FALSE]
  mode(new.disj) <- "integer"

  NA.ind <- which(rowSums(!is.na(new.disj)) == 0)
  if(length(NA.ind) > 0) new.disj <- new.disj[-NA.ind,,drop=FALSE]

  return(fitLinearLogicModel(model$X, model$y, Z = model$Z, disj = new.disj, Z.interactions = TRUE, type = type, s = s, ...))
}

#' Gene-environment (GxE) interaction test based on boosted linear models
#'
#' This function takes a fitted \code{linear.logic} model and independent test
#' data as input for testing if there is a general GxE interaction.
#' This hypothesis test is based on a likelihood-ratio test.
#'
#' In detail, the null hypothesis
#' \deqn{H_0: \delta_1 = \ldots = \delta_B = 0}
#' using the supplied linear model
#' \deqn{g(E[Y]) = \beta_0 + \sum_{i=1}^B \beta_i \cdot 1[C_i] + \delta_0 \cdot E
#' + \sum_{i=1}^B \delta_i \cdot 1[C_i] \cdot E}
#' is tested.
#'
#' @param model A fitted \code{linear.logic} model (i.e., a model created via
#'   \code{\link{fitLinearLogicModel}} or \code{\link{fitLinearBoostingModel}})
#' @param X Matrix or data frame of binary input data.
#'   This object should correspond to the binary matrix for fitting the model.
#' @param y Response vector. 0-1 coding for binary outcomes.
#' @param Z Quantitative covariable supplied as a matrix or data frame
#' @return A list containing
#'   \item{\code{Deviance}}{The deviance used for performing the
#'     likelihood-ratio test}
#'   \item{\code{p.value}}{The p-value of the test}
#'
#' @importFrom stats anova
#' @export
gxe.test.boosting <- function(model, X, y, Z) {
  if(!inherits(model, "linear.logic")) stop("The supplied model has to be of class 'linear.logic'!")

  if(model$type == "standard") {
    mod.complete <- fitLinearLogicModel(X, y, Z, model$disj, Z.interactions = TRUE, type = "standard")$lin.mod
    mod.reduced <- fitLinearLogicModel(X, y, Z, model$disj, Z.interactions = FALSE, type = "standard")$lin.mod
  } else {
    included.vars <- get.included.vars(model)
    mod.complete <- fitRestrictedLinearLogicModel(X, y, Z, included.vars$main.disj, included.vars$Z.disj)
    mod.reduced <- fitRestrictedLinearLogicModel(X, y, Z, included.vars$main.disj, NULL)
  }

  res <- anova(mod.reduced, mod.complete, test="LRT")
  ret <- list(Deviance = res$Deviance[2], p.value = res$`Pr(>Chi)`[2])
  return(ret)
}

#' Term importance test based on boosted linear models
#'
#' This function takes a fitted \code{linear.logic} model and independent test
#' data as input for testing if the included terms are influential with respect
#' to the outcome.
#' This hypothesis test is based on a likelihood-ratio test.
#'
#' In detail, the null hypotheses
#' \deqn{H_0: \beta_j = \delta_j = 0}
#' using the linear model
#' \deqn{g(E[Y]) = \beta_0 + \sum_{i=1}^B \beta_i \cdot 1[C_i] + \delta_0 \cdot E
#' + \sum_{i=1}^B \delta_i \cdot 1[C_i] \cdot E}
#' are tested for each \eqn{j \in \lbrace 1,\ldots,B \rbrace}
#' if \code{Z.interactions} is set to \code{TRUE}.
#' Otherwise, the null hypotheses
#' \deqn{H_0: \beta_j = 0}
#' using the linear model
#' \deqn{g(E[Y]) = \beta_0 + \sum_{i=1}^B \beta_i \cdot 1[C_i] + \delta_0 \cdot E}
#' are tested.
#'
#' @param model A fitted \code{linear.logic} model (i.e., a model created via
#'   \code{\link{fitLinearLogicModel}} or \code{\link{fitLinearBoostingModel}})
#' @param X Matrix or data frame of binary input data.
#'   This object should correspond to the binary matrix for fitting the model.
#' @param y Response vector. 0-1 coding for binary outcomes.
#' @param Z Optional quantitative covariables supplied as a matrix or
#'   data frame. Only used (and required) if the model was fitted using them.
#' @param Z.interactions A Boolean value determining whether interactions with
#'   quantitative covaraible \code{Z} shall be taken into account
#' @return A data frame consisting of three columns,
#'   \item{\code{var}}{The tested term,}
#'   \item{\code{vim}}{The associated variable importance, and}
#'   \item{\code{p.value}}{The corresponding p-value for testing if the term
#'     is influential.}
#'
#' @importFrom stats anova
#' @export
importance.test.boosting <- function(model, X, y, Z, Z.interactions = TRUE) {
  if(!inherits(model, "linear.logic")) stop("The supplied model has to be of class 'linear.logic'!")

  if(model$type == "standard") {
    disj <- model$disj
    n_terms <- nrow(disj)
    mod.complete <- fitLinearLogicModel(X, y, Z, disj, Z.interactions = Z.interactions, type = "standard")$lin.mod
  } else {
    included.vars <- get.included.vars(model)
    if(!Z.interactions) included.vars$Z.disj <- NULL
    main.disj <- included.vars$main.disj; Z.disj <- included.vars$Z.disj

    disj <- unique(rbind(main.disj, Z.disj))
    n_terms <- nrow(disj)

    mod.complete <- fitRestrictedLinearLogicModel(X, y, Z, main.disj, Z.disj)
  }

  vims <- data.frame(matrix(nrow = n_terms, ncol = 3))
  colnames(vims) <- c("var", "vim", "p.value")
  if(n_terms == 0) return(vims)
  real_disj <- translateLogicPET(disj, X)
  vims$var <- getPredictorNames(real_disj, sort_conj = TRUE)
  vims$vim <- 0 -> vims$p.value

  for(i in 1:n_terms) {
    if(model$type == "standard") {
      mod.reduced <- fitLinearLogicModel(X, y, Z, disj[-i,,drop=FALSE], Z.interactions = Z.interactions, type = "standard")$lin.mod
    } else {
      main.rem <- which(apply(main.disj, 1, function(x) all.equal(x, disj[i,])) == "TRUE")
      if(length(main.rem) > 0)
        main.disj.rem <- main.disj[-main.rem,,drop=FALSE]
      else
        main.disj.rem <- main.disj

      if(!is.null(Z.disj)) {
        Z.rem <- which(apply(Z.disj, 1, function(x) all.equal(x, disj[i,])) == "TRUE")
        if(length(Z.rem) > 0)
          Z.disj.rem <- Z.disj[-Z.rem,,drop=FALSE]
        else
          Z.disj.rem <- Z.disj
      } else {
        Z.disj.rem <- NULL
      }

      mod.reduced <- fitRestrictedLinearLogicModel(X, y, Z, main.disj.rem, Z.disj.rem)
    }

    res <- anova(mod.reduced, mod.complete, test="LRT")
    vims$vim[i] <- res$Deviance[2]
    vims$p.value[i] <- res$`Pr(>Chi)`[2]
  }

  return(vims)
}

#' Tuning the LASSO regularization parameter
#'
#' This function takes a fitted \code{linear.logic} model and independent
#' validation data as input for finding the ideal LASSO complexity penalty
#' \code{s}.
#'
#' @param model A fitted \code{linear.logic} model (i.e., a model created via
#'   \code{\link{fitLinearLogicModel}} or \code{\link{fitLinearBoostingModel}})
#' @param X Matrix or data frame of binary input data.
#'   This object should correspond to the binary matrix for fitting the model.
#' @param y Response vector. 0-1 coding for binary outcomes.
#' @param Z Optional quantitative covariables supplied as a matrix or
#'   data frame. Only used (and required) if the model was fitted using them.
#' @param scoring_rule The scoring rule for evaluating the validation
#'   error and its standard error. For classification tasks, \code{"deviance"}
#'   or \code{"Brier"} should be used.
#' @param choose Model selection scheme. If the model that minimizes the
#'   validation error should be chosen, \code{choose = "min"} should be
#'   set. Otherwise, \code{choose = "1se"} leads to simplest model in the range
#'   of one standard error of the minimizing model.
#' @return A list containing
#'   \item{\code{val.res}}{A data frame containing the penalties, the
#'   validation scores and the corresponding standard errors}
#'   \item{\code{best.s}}{The ideal penalty value}
#'
#' @export
get.ideal.penalty <- function(model, X, y, Z = NULL, scoring_rule = "deviance", choose = "min") {
  if(model$type != "lasso") stop("Only applicable for LASSO fits!")
  lambda <- model$lin.mod$lambda
  beta <- model$lin.mod$beta
  preds <- predict(model, X, Z = Z)
  y_bin <- setequal(unique(y), c(0, 1))
  if(!y_bin) scoring_rule <- "mse"

  val.res <- data.frame()
  for(i in 1:ncol(beta)) {
    scores <- calcScorePerObservation(as.numeric(preds[,i]), y, scoring_rule)
    m <- mean(scores)
    se <- sd(scores)/sqrt(length(scores))
    val.res <- rbind(val.res, data.frame(s = lambda[i], score = m, se = se, score.plus.1se = m + se))
  }

  min.ind <- min(which(val.res$score == min(val.res$score)))
  if(choose == "1se") {
    max.val <- val.res$score.plus.1se[min.ind]
    min.ind <- min(which(val.res$score <= max.val))
  }

  best.s <- lambda[min.ind]
  return(list(val.res = val.res, best.s = best.s))
}

get.included.vars <- function(model) {
  if(is.null(model$s)) stop("The regularization parameter $s has to be properly set!")
  lin.mod <- model$lin.mod
  if(model$type == "cv.lasso") lin.mod <- lin.mod$glmnet.fit
  lambda.ind <- match(model$s, lin.mod$lambda)[1]
  beta <- as.numeric(lin.mod$beta[,lambda.ind])

  n_conj <- sum(rowSums(!is.na(model$disj)) > 0)
  main.disj <- model$disj[beta[1:n_conj] != 0,,drop=FALSE]
  Z.disj <- NULL
  if(length(beta) > n_conj + 1)
    Z.disj <- model$disj[beta[(n_conj+1):(2*n_conj)] != 0,,drop=FALSE]

  return(list(main.disj = main.disj, Z.disj = Z.disj))
}

fitRestrictedLinearLogicModel <- function(X, y, Z = NULL, main.disj, Z.disj) {
  y_bin <- setequal(unique(y), c(0, 1))
  X <- as.matrix(X)
  mode(X) <- "integer"
  mode(main.disj) <- "integer"
  dm <- getDesignMatrix(X, main.disj)

  Z.interactions <- !is.null(Z.disj) && nrow(Z.disj) > 0 && !is.null(Z)
  if(Z.interactions) mode(Z.disj) <- "integer"

  if(!is.null(Z)) {
    Z <- as.numeric(as.matrix(Z))
    if(Z.interactions)
      dm <- cbind(dm, getDesignMatrix(X, Z.disj) * Z, Z = Z)
    else
      dm <- cbind(dm, Z = Z)
  }

  if(nrow(main.disj) > 0)
    colnames(dm)[1:nrow(main.disj)] <- paste("X", 1:nrow(main.disj), sep="")
  if(Z.interactions)
    colnames(dm)[(nrow(main.disj)+1):(nrow(main.disj)+nrow(Z.disj))] <- paste("Z", 1:nrow(Z.disj), sep="")

  if(y_bin) fam <- binomial(link = "logit") else fam <- gaussian()
  dat <- data.frame(y = y, dm)
  lin.mod <- glm(y ~ ., data = dat, family = fam)
  return(lin.mod)
}




