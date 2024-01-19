#' @export
#' @rawNamespace S3method(logicDT,bagging)
logicDT.bagging <- function(X, ...) UseMethod("logicDT.bagging")

#' Fitting bagged logicDT models
#'
#' Function for fitting bagged logicDT models.
#'
#' Details on single logicDT models can be found in \code{\link{logicDT}}.
#'
#' @param X Matrix or data frame of binary predictors coded as 0 or 1.
#' @param y Response vector. 0-1 coding for binary responses.
#'  Otherwise, a regression task is assumed.
#' @param Z Optional matrix or data frame of quantitative/continuous
#'   covariables. Multiple covariables allowed for splitting the trees.
#'   If leaf regression models (such as four parameter logistic models) shall
#'   be fitted, only the first given covariable is used.
#' @param bagging.iter Number of bagging iterations
#' @param ... Arguments passed to \code{logicDT}
#' @return An object of class \code{logic.bagged}. This is a list
#'   containing
#'   \item{\code{models}}{A list of fitted \code{logicDT} models}
#'   \item{\code{bags}}{A list of observation indices which were
#'   used to train each model}
#'   \item{\code{...}}{Supplied parameters of the functional call
#'     to \code{\link{logicDT.bagging}}.}
#'
#' @name logicDT.bagging
#' @method logicDT.bagging default
#' @importFrom utils flush.console
#' @export
logicDT.bagging.default <- function(X, y, Z = NULL, bagging.iter = 500, ...) {
  ret <- list()
  ret$models <- list()
  ret$bags <- list()
  ret$X <- X
  ret$y <- y
  ret$Z <- Z
  ret$total_iter <- 0
  for(i in 1:bagging.iter) {
    bag <- sample(nrow(X), nrow(X), replace=TRUE)
    ret$bags[[i]] <- bag
    if(is.null(Z))
      Z_temp <- NULL
    else
      Z_temp <- Z[bag,,drop=FALSE]
    ret$models[[i]] <- logicDT(X[bag,,drop=FALSE], y[bag], Z = Z_temp, ...)
    ret$total_iter <- ret$total_iter + ret$models[[i]]$total_iter
    cat("\r", sprintf("Iteration %d/%d (%.0f%%)", i, bagging.iter, i/bagging.iter * 100))
    flush.console()
  }
  cat("\n")
  class(ret) <- "logic.bagged"
  return(ret)
}

#' @param formula An object of type \code{formula} describing the
#'   model to be fitted.
#' @param data A data frame containing the data for the corresponding
#'   \code{formula} object. Must also contain quantitative covariables
#'   if they should be included as well.
#' @rdname logicDT.bagging
#' @name logicDT.bagging
#' @method logicDT.bagging formula
#' @importFrom stats model.frame model.response model.matrix
#' @export
logicDT.bagging.formula <- function(formula, data, ...) {
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  predictors <- model.matrix(formula, data = mf)[,-1]
  quant.preds <- apply(predictors, 2, function(x) any(!(x %in% c(0, 1))))
  X <- predictors[,!quant.preds,drop=FALSE]
  Z <- predictors[,quant.preds,drop=FALSE]
  if(ncol(Z) == 0) Z <- NULL
  logicDT.bagging(X, y, Z = Z, ...)
}

#' @rdname predict.logicDT
#' @name predict.logicDT
#' @export
predict.logic.bagged <- function(object, X, Z = NULL, type = "prob", ...) {
  preds <- rep(0, nrow(X))
  for(i in 1:length(object$models)) {
    preds <- preds + predict(object$models[[i]], X, Z, ...)
  }
  preds <- preds/length(object$models)
  if(type == "class") {
    preds <- as.numeric(preds > 0.5)
  }
  return(preds)
}




