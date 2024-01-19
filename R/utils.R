# Tune classification threshold
tune.threshold <- function(preds, y) {
  thresholds <- seq(0, 1, 0.001)
  best_eval <- 0
  best_thresholds <- numeric()
  for(border in thresholds) {
    eval <- mean(as.numeric(preds > border) == y)
    if(eval == best_eval) {
      best_thresholds <- c(best_thresholds, border)
    } else if(eval > best_eval) {
      best_thresholds <- border
      best_eval <- eval
    }
  }
  best_thresholds <- best_thresholds[order(abs(best_thresholds - 0.5))]
  return(list(best_thresholds = best_thresholds, best_eval = best_eval))
}

#' Calculate the deviance
#'
#' Computation of the deviance, i.e., two times the negative
#' log likelihood for risk estimates in a binary classification
#' problem.
#'
#' @param preds Numeric vector of predictions
#' @param y True outcomes
#' @return The deviance
#'
#' @export
calcDev <- function(preds, y) {
  return(-2*sum(log(y * preds + (1-y) * (1-preds))))
}

#' Calculate the Brier score
#'
#' Computation of the Brier score, i.e., the mean squared
#' error for risk estimates in a binary classification
#' problem.
#'
#' @param preds Numeric vector of predictions
#' @param y True outcomes
#' @return The Brier score
#'
#' @export
calcBrier <- function(preds, y) {
  return(mean((preds - y)^2))
}

#' Calculate the misclassification rate
#'
#' Computation of the misclassification rate for risk
#' estimates in a binary classification problem.
#'
#' @param preds Numeric vector of predictions
#' @param y True outcomes
#' @param cutoff Classification cutoff. By default,
#'   scores above 50% are assigned class 1 and class 0
#'   otherwise.
#' @return The misclassification rate
#'
#' @export
calcMis <- function(preds, y, cutoff = 0.5) {
  return(mean((preds > cutoff) != y))
}

#' Calculate the normalized cross entropy
#'
#' This function computes the normalized cross entropy (NCE)
#' which is given by
#' \deqn{\mathrm{NCE} = \frac{\frac{1}{N} \sum_{i=1}^{N}
#' y_i \cdot \log(p_i) + (1-y_i) \cdot \log(1-p_i)}{
#' p \cdot \log(p) + (1-p) \cdot \log(1-p)}}
#' where (for \eqn{i \in \lbrace 1,\ldots,N \rbrace})
#' \eqn{y_i \in \lbrace 0,1 \rbrace} are the true classes,
#' \eqn{p_i} are the risk/probability predictions and
#' \eqn{p = \frac{1}{N} \sum_{i=1}^{N} y_i} is total unrestricted
#' empirical risk estimate.
#'
#' Smaller values towards zero are generally prefered.
#' A NCE of one or above would indicate that the used model
#' yields comparable or worse predictions than the naive mean
#' model.
#'
#' @references
#' \itemize{
#'   \item He, X., Pan, J., Jin, O., Xu, T., Liu, B., Xu, T.,
#'     Shi, Y., Atallah, A., Herbrich, R., Bowers, S., Candela, J. Q.
#'     (2014). Practical Lessons from Predicting Clicks on Ads at
#'     Facebook. Proceedings of the Eighth International Workshop on
#'     Data Mining for Online Advertising 1-9.
#'     \doi{https://doi.org/10.1145/2648584.2648589}
#' }
#'
#' @param preds Numeric vector of risk estimates
#' @param y Vector of true binary outcomes
#' @return The normalized cross entropy
#'
#' @export
calcNCE <- function(preds, y) {
  ce <- mean(log(y * preds + (1-y) * (1-preds)))
  p <- mean(y)
  denominator <- p * log(p) + (1-p) * log(1-p)
  return(ce/denominator)
}

#' Fast computation of the AUC w.r.t. to the ROC
#'
#' This function computes the area under the receiver operating
#' characteristic curve.
#'
#' @param preds Numeric vector of predicted scores
#' @param y True binary outcomes coded as 0 or 1. Must be an integer
#'   vector.
#' @param fast Shall the computation be as fast as possible?
#' @param sorted Are the predicted scores already sorted
#'   increasingly? If so, this can slightly speed up the computation.
#' @return The AUC between 0 and 1
#'
#' @export
calcAUC <- function(preds, y, fast = TRUE, sorted = FALSE) {
  if(length(preds) > 90000) {
    stop("AUC computation is only supported for up to N = 90000!")
  }
  if(fast)
    return(calcAUCFast(preds, y, sorted))
  else {
    score.pairs <- merge(preds[y == 1], preds[y == 0])
    auc <- mean((score.pairs[,1] > score.pairs[,2]) + 0.5 * (score.pairs[,1] == score.pairs[,2]))
    return(auc)
  }
}

#' Calculate the MSE
#'
#' Computation of the mean squared error.
#'
#' @param preds Numeric vector of predictions
#' @param y True outcomes
#' @return The MSE
#'
#' @export
calcMSE <- function(preds, y) {
  return(mean((preds - y)^2))
}

#' Calculate the NRMSE
#'
#' Computation of the normalized root mean squared error.
#'
#' @param preds Numeric vector of predictions
#' @param y True outcomes
#' @param type \code{"sd"} uses the standard deviation of \code{y} for
#'   normalization. \code{"range"} uses the whole span of \code{y}.
#' @return The NRMSE
#'
#' @importFrom stats sd
#' @export
calcNRMSE <- function(preds, y, type = "sd") {
  if(type == "range") {
    return(sqrt(calcMSE(preds, y))/diff(range(y)))
  } else {
    return(sqrt(calcMSE(preds, y))/sd(y))
  }
}

#' Split biallelic SNPs into binary variables
#'
#' This function takes a matrix or data frame of SNPs coded as
#' 0, 1, 2 or 1, 2, 3 and returns a data frame with twice as many
#' columns. SNPs are splitted into dominant and recessive modes,
#' i.e., for a \eqn{\mathrm{SNP} \in \lbrace 0,1,2 \rbrace}, two variables
#' \eqn{\mathrm{SNP}_D = (\mathrm{SNP} \neq 0)} and
#' \eqn{\mathrm{SNP}_R = (\mathrm{SNP} = 2)} are generated.
#'
#' @param data A matrix or data frame only consisting of SNPs to
#'   be splitted
#' @return A data frame of the splitted SNPs
#'
#' @export
splitSNPs <- function(data) {
  modified_data <- data.frame(matrix(ncol = 0, nrow = nrow(data)))
  coding <- range(data)

  for (i in 1:ncol(data)) {
    xd <- ifelse(data[,i] != coding[1], 1, 0)
    xr <- ifelse(data[,i] == coding[2], 1, 0)
    modified_data[paste(colnames(data)[i], "D", sep="")] <- xd
    modified_data[paste(colnames(data)[i], "R", sep="")] <- xr
  }
  return(modified_data)
}

