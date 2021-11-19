#' @export
print.logicDT <- function(x, ...) {
  cat(getPredictorNames(x$real_disj))
  cat("\n")
}

#' @export
print.logic.bagged <- function(x, ...) {
  for(m in x$models) {
    print(m); cat("\n")
  }
}

#' @export
print.logic.boosted <- function(x, ...) {
  for(m in x$models) {
    print(m); cat("\n")
  }
}

#' @export
print.vim <- function(x, ...) {
  print(x$vims)
}

