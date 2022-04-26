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
  if(model$y_bin) {
    y <- as.integer(y)
    if(scoring_rule == "deviance")
      eval.performance <- function(preds, y) -calcDev(preds, y)
    else if(scoring_rule == "brier")
      eval.performance <- function(preds, y) -calcBrier(preds, y)
    else if(scoring_rule == "nce")
      eval.performance <- function(preds, y) 1-calcNCE(preds, y)
    else
      eval.performance <- function(preds, y) calcAUC(preds, y)
  } else {
    eval.performance <- function(preds, y) -calcMSE(preds, y)
  }
  for(j in 1:boosting.iter) {
    current.perf <- eval.performance(partial.preds[,j], y)
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

#' Linear models based on logic terms
#'
#' This function fits a linear or logistic regression model (based on the
#' type of outcome) using the supplied logic terms, e.g., \code{$disj} from
#' a fitted \code{logicDT} model.
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
#' @return A \code{linear.logic} model. This is a list containing
#'   the logic terms used as predictors in the model and the fitted \code{glm}
#'   model.
#'
#' @importFrom stats glm gaussian binomial
#' @export
fitLinearLogicModel <- function(X, y, Z = NULL, disj, Z.interactions = TRUE) {
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
  if(y_bin) fam <- binomial(link = "logit") else fam <- gaussian()
  dat <- data.frame(y = y, dm)
  lin.mod <- glm(y ~ ., data = dat, family = fam)
  ret <- list(disj = disj, lin.mod = lin.mod, Z.interactions = Z.interactions)
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
#' @param ... Ignored additional parameters
#' @return A numeric vector of predictions. For binary outcomes,
#'   this is a vector with estimates for
#'   \eqn{P(Y=1 \mid X = x)}.
#'
#' @importFrom stats predict.glm
#' @export
predict.linear.logic <- function(object, X, Z = NULL, ...) {
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
  as.numeric(predict(object$lin.mod, data.frame(dm), type = "response"))
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
#' @param model Fitted \code{logic.boosted} model
#' @param n.iter Number of boosting iterations to be used
#' @return A \code{linear.logic} model. This is a list containing
#'   the logic terms used as predictors in the model and the fitted \code{glm}
#'   model.
#'
#' @export
fitLinearBoostingModel <- function(model, n.iter) {
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

  return(fitLinearLogicModel(model$X, model$y, Z = model$Z, disj = new.disj, Z.interactions = TRUE))
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
  if(class(model) != "linear.logic") stop("The supplied model has to be of class 'linear.logic'!")

  mod.complete <- fitLinearLogicModel(X, y, Z, model$disj, Z.interactions = TRUE)
  mod.reduced <- fitLinearLogicModel(X, y, Z, model$disj, Z.interactions = FALSE)

  res <- anova(mod.reduced$lin.mod, mod.complete$lin.mod, test="LRT")
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
  if(class(model) != "linear.logic") stop("The supplied model has to be of class 'linear.logic'!")

  Z.here <- length(model$lin.mod$coefficients) > (nrow(model$disj) + 1)
  disj <- model$disj
  n_terms <- nrow(disj)
  mod.complete <- fitLinearLogicModel(X, y, Z, disj, Z.interactions = Z.interactions)

  vims <- data.frame(matrix(nrow = n_terms, ncol = 3))
  colnames(vims) <- c("var", "vim", "p.value")
  real_disj <- translateLogicPET(disj, X)
  vims$var <- getPredictorNames(real_disj, sort_conj = TRUE)
  vims$vim <- 0 -> vims$p.value

  for(i in 1:n_terms) {
    mod.reduced <- fitLinearLogicModel(X, y, Z, disj[-i,,drop=FALSE], Z.interactions = Z.interactions)
    res <- anova(mod.reduced$lin.mod, mod.complete$lin.mod, test="LRT")
    vims$vim[i] <- res$Deviance[2]
    vims$p.value[i] <- res$`Pr(>Chi)`[2]
  }

  return(vims)
}




