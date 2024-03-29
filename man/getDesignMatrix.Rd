% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logicDT.R
\name{getDesignMatrix}
\alias{getDesignMatrix}
\title{Design matrix for the set of conjunctions}
\usage{
getDesignMatrix(X, disj)
}
\arguments{
\item{X}{The original (binary) predictor matrix. This has to be of type
\code{integer}.}

\item{disj}{The conjunction matrix which can, e.g., be extracted from a
fitted \code{logicDT} model via \code{$disj}.}
}
\value{
The transformed design matrix.
}
\description{
Transform the original predictor matrix X into the conjunction design matrix
which contains for each conjunction a corresponding column.
}
