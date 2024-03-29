% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logic.bagging.R
\name{logicDT.bagging}
\alias{logicDT.bagging}
\alias{logicDT.bagging.default}
\alias{logicDT.bagging.formula}
\title{Fitting bagged logicDT models}
\usage{
\method{logicDT.bagging}{default}(X, y, Z = NULL, bagging.iter = 500, ...)

\method{logicDT.bagging}{formula}(formula, data, ...)
}
\arguments{
\item{X}{Matrix or data frame of binary predictors coded as 0 or 1.}

\item{y}{Response vector. 0-1 coding for binary responses.
Otherwise, a regression task is assumed.}

\item{Z}{Optional matrix or data frame of quantitative/continuous
covariables. Multiple covariables allowed for splitting the trees.
If leaf regression models (such as four parameter logistic models) shall
be fitted, only the first given covariable is used.}

\item{bagging.iter}{Number of bagging iterations}

\item{...}{Arguments passed to \code{logicDT}}

\item{formula}{An object of type \code{formula} describing the
model to be fitted.}

\item{data}{A data frame containing the data for the corresponding
\code{formula} object. Must also contain quantitative covariables
if they should be included as well.}
}
\value{
An object of class \code{logic.bagged}. This is a list
  containing
  \item{\code{models}}{A list of fitted \code{logicDT} models}
  \item{\code{bags}}{A list of observation indices which were
  used to train each model}
  \item{\code{...}}{Supplied parameters of the functional call
    to \code{\link{logicDT.bagging}}.}
}
\description{
Function for fitting bagged logicDT models.
}
\details{
Details on single logicDT models can be found in \code{\link{logicDT}}.
}
