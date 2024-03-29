% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logicDT.R
\name{fitLinearModel}
\alias{fitLinearModel}
\title{Fitting linear models}
\usage{
fitLinearModel(y, Z, logistic = TRUE)
}
\arguments{
\item{y}{Response vector. 0-1 coding for binary outcomes,
otherwise conventional regression is performed.}

\item{Z}{Numeric vector of (univariate) input samples.}

\item{logistic}{Logical indicating whether, in the case of a binary
outcome, a logistic regression model should be fitted
(\code{TRUE}) or a LDA model should be fitted (\code{FALSE})}
}
\value{
An object of class \code{"linear"} which contains a numeric
  vector of the fitted parameters b and c.
}
\description{
Method for fitting linear models.
In the fashion of this package, only binary and quantitative
outcomes are supported.
}
\details{
For binary outcomes, predictions are cut at 0 or 1 for generating
proper probability estimates.
}
