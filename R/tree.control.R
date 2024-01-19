#' Control parameters for fitting decision trees
#'
#' Configure the fitting process of individual decision trees.
#'
#' For the Gini or MSE splitting criterion,
#' if any considered split \eqn{s} leads to
#' \deqn{P(t) \cdot \Delta I(s,t) > \texttt{cp}}
#' for a node \eqn{t}, the empirical node probability
#' \eqn{P(t)} and the impurity reduction \eqn{\Delta I(s,t)},
#' then the node is further splitted. If not, the node is
#' declared as a leaf.
#' For continuous outcomes, \code{cp} will be scaled by the
#' empirical variance of \code{y} to ensure the right scaling,
#' i.e., \code{cp <- cp * var(y)}. Since the impurity measure
#' for continuous outcomes is the mean squared error, this can
#' be interpreted as controlling the minimum reduction of the
#' normalized mean squared error (NRMSE to the power of two).
#'
#' If one chooses the 4pL or linear splitting criterion, likelihood
#' ratio tests testing the alternative of better fitting individual
#' models are employed. The corresponding test statistic
#' asymptotically follows a \eqn{\chi^2} distribution where
#' the degrees of freedom are given by the difference in the
#' number of model parameters, i.e., leading to
#' \eqn{2 \cdot 4 - 4 = 4} degrees of freedom in the case of 4pL
#' models and to \eqn{2 \cdot 2 - 2 = 2} degrees of freedom in
#' the case of linear models.
#'
#' For binary outcomes, choosing to fit linear models for evaluating
#' the splits or for modeling the leaves actually leads to fitting
#' LDA (linear discriminant analysis) models.
#'
#' @param nodesize Minimum number of samples contained in a
#'   terminal node. This parameter ensures that enough samples
#'   are available for performing predictions which includes
#'   fitting regression models such as 4pL models.
#' @param split_criterion Splitting criterion for deciding
#'   when and how to split. The default is \code{"gini"}/\code{"mse"} which
#'   utilizes the Gini splitting criterion for binary risk
#'   estimation tasks and the mean squared error as impurity
#'   measure in regression tasks. Alternatively, \code{"4pl"} can be
#'   used if a quantitative covariable is supplied and
#'   the parameter \code{covariable} is chosen such that 4pL
#'   model fitting is enabled, i.e., \code{covariable = "final_4pl"}
#'   or \code{covariable = "full_4pl"}.
#'   A fast modeling alternative is given by \code{"linear"} which also
#'   requires the parameter \code{covariable} to be properly
#'   chosen, i.e., \code{covariable = "final_linear"}
#'   or \code{covariable = "full_linear"}.
#' @param alpha Significance threshold for the likelihood ratio
#'   tests when using \code{split_criterion = "4pl"} or \code{"linear"}.
#'   Only splits that achieve a p-value smaller than \code{alpha} are eligible.
#' @param cp Complexity parameter. This parameter determines
#'   by which amount the impurity has to be reduced to further
#'   split a node. Here, the total tree impurity is considered.
#'   See details for a specific formula. Only used if
#'   \code{split_criterion = "gini"} or \code{"mse"}.
#' @param smoothing Shall the leaf predictions for risk
#'   estimation be smoothed? \code{"laplace"} yields Laplace smoothing.
#'   The default is \code{"none"} which does not employ smoothing.
#' @param mtry Shall the tree fitting process be randomized
#'   as in random forests? Currently, only \code{"sqrt"} for using
#'   \eqn{\sqrt{p}} random predictors at each node for splitting
#'   and \code{"none"} (default) for fitting conventional decision trees
#'   are supported.
#' @param covariable How shall optional quantitative covariables
#'   be handled? \code{"constant"} ignores them. Alternatively,
#'   they can be considered as splitting variables (\code{"_split"}),
#'   used for fitting 4pL models in each leaf (\code{"_4pl"}), or used
#'   for fitting linear models in each leaf (\code{"_linear"}). If either
#'   splitting or model fitting is chosen, one should state if this
#'   should be handled over the whole search (\code{"full_"},
#'   computationally expensive) or just the final trees
#'   (\code{"final_"}). Thus, \code{"final_4pl"} would lead to fitting
#'   4pL models in each leaf but only for the final tree fitting.
#' @return An object of class \code{tree.control} which is a list
#'   of all necessary tree parameters.
#'
#' @export
tree.control <- function(nodesize = 10, split_criterion = "gini", alpha = 0.05, cp = 0.001,
                         smoothing = "none", mtry = "none", covariable = "final_4pl") {
  if(smoothing == "laplace")
    smoothing <- 1
  else
    smoothing <- 0
  if(mtry == "none")
    mtry <- -1
  else if(mtry == "sqrt")
    mtry <- 0

  if(split_criterion == "linear")
    split_criterion <- 2
  else if(split_criterion == "4pl")
    split_criterion <- 1
  else
    split_criterion <- 0

  if(covariable == "full_linear")
    covariable <- 3
  else if(covariable == "final_linear")
    covariable <- -3
  else if(covariable == "full_4pl")
    covariable <- 2
  else if(covariable == "final_4pl")
    covariable <- -2
  else if(covariable == "full_split")
    covariable <- 1
  else if(covariable == "final_split")
    covariable <- -1
  else
    covariable <- 0

  if(covariable == 0 && split_criterion > 0)
    covariable <- split_criterion + 1

  covariable_search <- max(covariable, 0)
  covariable_final <- abs(covariable)
  # Logistic models instead of LDA models in the final fitting:
  if(covariable_final == 3)
    covariable_final <- 4
  tc <- list(nodesize = as.integer(nodesize),
             split_criterion = as.integer(split_criterion), alpha = as.numeric(alpha),
             cp = as.numeric(cp), cp_orig = as.numeric(cp),
             smoothing = as.integer(smoothing),
             mtry = as.integer(mtry), covariable = as.integer(covariable),
             covariable_search = as.integer(covariable_search), covariable_final = as.integer(covariable_final))
  class(tc) <- "tree.control"
  return(tc)
}

