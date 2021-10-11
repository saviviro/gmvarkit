#' @title Deprecated S3 methods for the deprecated class 'gmvar'
#'
#' @description Deprecated S3 methods for the deprecated class 'gmvar'. From
#'   the gmvarkit version 2.0.0 onwards, class 'gsmvar' is used instead.
#'
#' @param gmvar a class 'gmvar' object. THIS CLASS IS DEPRECATED FROM THE VERSION
#'   2.0.0 ONWARDS.
#' @param digits number of digits to be printed.
#' @details These methods exists so that models estimated with earlier versions
#'   of the package can be used.
#' @export

print.gmvar <- function(gmvar, ..., digits=2) {
  class(gmvar) <- "gsmvar"
  print(gmvar, digits=digits)
}

#' @rdname print.gmvar
#' @export

summary.gmvar <- function(gmvar, ...., digits=2) {
  class(gmvar) <- "gsmvar"
  summary(gmvar, digits=digits)
}

#' @rdname print.gmvar
#' @export

summary.gmvar <- function(gmvar) {
  class(gmvar) <- "gsmvar"
  plot(gmvar)
}

# To do: gvmarpred.plot/print, print.gmvarsum,
