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




