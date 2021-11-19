#' Control parameters for fitting decision trees
#'
#' Configure the fitting process of individual decision trees.
#'
#' If any considered split \eqn{s} leads to
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
#' @param nodesize Minimum number of samples contained in a
#'   terminal node. This parameter ensures that enough samples
#'   are available for performing predictions which includes
#'   fitting 4pL models.
#' @param cp Complexity parameter. This parameter determines
#'   by which amount the impurity has to be reduced to further
#'   split a node. Here, the total tree impurity is considered.
#'   See details for a concrete formula.
#' @param smoothing Shall the leaf predictions for risk
#'   estimation be smoothed? "laplace" yields Laplace smoothing.
#'   The default is "none" which does not employ smoothing.#'
#' @param mtry Shall the tree fitting process be randomized
#'   as in random forests? Currently, only "sqrt" for using
#'   \eqn{\sqrt{p}} random predictors at each node for splitting
#'   and "none" (default) for fitting conventional decision trees
#'   are supported.
#' @param covariable How shall optional quantitative covariables
#'   be handled? "constant" ignores them. Alternatively,
#'   they can be considered as splitting variables ("_split") or
#'   used for fitting 4pL in each leaf ("_4pl"). If either
#'   splitting or 4pL is chosen, one should state if this
#'   should be handled over the whole search ("full_",
#'   computationally expensive) or just the final trees
#'   ("final_"). Thus, "final_4pl" would lead to fitting 4pL in
#'   each leaf but only for the final fitting of trees.
#' @return An object of class \code{tree.control} which is a list
#'   of all necessary tree parameters.
#'
#' @export
tree.control <- function(nodesize = 10, cp = 0.001, smoothing = "none", mtry = "none", covariable = "final_4pl") {
  if(smoothing == "laplace")
    smoothing <- 1
  else
    smoothing <- 0
  if(mtry == "none")
    mtry <- -1
  else if(mtry == "sqrt")
    mtry <- 0

  if(covariable == "full_4pl")
    covariable <- 2
  else if(covariable == "final_4pl")
    covariable <- -2
  else if(covariable == "full_split")
    covariable <- 1
  else if(covariable == "final_split")
    covariable <- -1
  else
    covariable <- 0
  covariable_search <- max(covariable, 0)
  covariable_final <- abs(covariable)
  tc <- list(nodesize = as.integer(nodesize), cp = as.numeric(cp), cp_orig = as.numeric(cp),
             smoothing = as.integer(smoothing),
             mtry = as.integer(mtry), covariable = as.integer(covariable),
             covariable_search = as.integer(covariable_search), covariable_final = as.integer(covariable_final))
  class(tc) <- "tree.control"
  return(tc)
}

