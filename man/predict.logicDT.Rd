% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logic.bagging.R, R/logic.boosting.R,
%   R/logicDT.R
\name{predict.logicDT}
\alias{predict.logicDT}
\alias{predict.logic.bagged}
\alias{predict.logic.boosted}
\alias{predict.geneticLogicPET}
\title{Prediction for logicDT models}
\usage{
\method{predict}{logic.bagged}(object, X, Z = NULL, type = "prob", ...)

\method{predict}{logic.boosted}(object, X, Z = NULL, type = "prob", ...)

\method{predict}{logicDT}(
  object,
  X,
  Z = NULL,
  type = "prob",
  ensemble = FALSE,
  leaves = "4pl",
  ...
)

\method{predict}{geneticLogicPET}(
  object,
  X,
  Z = NULL,
  models = "best",
  n_models = 10,
  ensemble = NULL,
  leaves = "4pl",
  ...
)
}
\arguments{
\item{object}{Fitted logicDT model. Usually a product of a call
to \code{\link{logicDT}}.}

\item{X}{Matrix or data frame of binary input data. This
object should correspond to the binary matrix for fitting
the model.}

\item{Z}{Optional quantitative covariables supplied as a
matrix or data frame. Only used (and required) if the
model was fitted using them.}

\item{type}{Prediction type. This can either be "prob" for
probability estimates or "class" for
classification in binary responses. Ignored for regression.}

\item{...}{Parameters supplied to \code{\link{predict.logicDT}}}

\item{ensemble}{If the model was fitted using the inner
validation approach, shall the prediction be constructed
using the final validated ensemble (TRUE) or using the
single final tree (FALSE)?}

\item{leaves}{If four parameter logistic models were fitted
for each leaf, shall they be used for the prediction
("4pl") or shall the constant leaf means be used
("constant")?}

\item{models}{Which models of logicDT model fitted with
genetic programming shall be used for prediction?
"best" leads to the single best model in the final
generation, "all" uses the average over the final
generation and "n_models" uses the n_models best models.}

\item{n_models}{How many models shall be used if
models = "n_models" and genetic programming was employed?}
}
\value{
A numeric vector of predictions. For binary outcomes,
  this is a vector with estimates for \eqn{P(Y=1 \mid X = x)}.
}
\description{
Supply new input data for predicting the outcome with a fitted
logicDT model.
}
