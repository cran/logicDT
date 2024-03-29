% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fancy.plot.R, R/plot.R
\name{plot.logicDT}
\alias{plot.logicDT}
\alias{fancy.plot}
\title{Plot a logic decision tree}
\usage{
fancy.plot(x, cdot = FALSE, ...)

\method{plot}{logicDT}(
  x,
  fancy = TRUE,
  x_scaler = 0.5,
  margin_scaler = 0.2,
  cex = 1,
  cdot = FALSE,
  ...
)
}
\arguments{
\item{x}{An object of the class \code{logicDT}}

\item{cdot}{Should a centered dot be used instead of a logical and
for depicting interactions?}

\item{...}{Arguments passed to fancy plotting function}

\item{fancy}{Should the fancy mode be used for plotting? Default is
\code{TRUE}.}

\item{x_scaler}{Scaling factor on the horizontal axis for deeper trees,
i.e., \code{x_scaler = 0.5} means that the horizontal distance
between two adjacent nodes is halved for every vertical level.}

\item{margin_scaler}{Margin factor. Smaller values lead to smaller
margins.}

\item{cex}{Scaling factor for the plotted text elements.}
}
\value{
No return value, called for side effects
}
\description{
This function plots a logicDT model on the active graphics device.
}
\details{
There are two plotting modes:
\itemize{
  \item \code{fancy = FALSE} which draws a tree with direct edges between
  the nodes. Leaves are represented by their prediction value which is
  obtained by the (observed) conditional mean.
  \item \code{fancy = TRUE} plots a tree similar to those in the \code{rpart}
  (Therneau and Atkinson, 2019) and \code{splinetree} (Neufeld and Heggeseth,
  2019) \code{R} packages. The trees are drawn in an angular manner and
  if leaf regression models were fitted, appropriate plots of the fitted
  curves are depicted in the leaves. Otherwise, the usual prediction values
  are shown.
}
}
\references{
\itemize{
  \item Therneau, T. & Atkinson, B. (2019). rpart: Recursive Partitioning
  and Regression Trees. \url{https://CRAN.R-project.org/package=rpart}
  \item Neufeld, A. & Heggeseth, B. (2019). splinetree: Longitudinal
  Regression Trees and Forests.
  \url{https://CRAN.R-project.org/package=splinetree}
}
}
