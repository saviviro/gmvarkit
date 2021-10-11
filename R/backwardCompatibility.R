#' @title Makes class 'gmvar' objects compatible with the functions using class 'gsmvar' objects
#'
#' @description \code{gmvar_to_gsmvar} class 'gmvar' objects compatible with the functions using
#' s class 'gsmvar' objects
#'
#' @param gsmvar a class 'gmvar' or 'gsmvar' object.
#' @details This exists so that models estimated with earlier versions
#'   of the package can be used conveniently.
#' @return If the provided object has the class 'gsmvar', the provided object
#'   is returned without modifications. If the provided object has the class 'gmvar',
#'   its element \code{$model} is given a new subelement called also model and this is
#'   set to be "GMVAR". Also, the class of this object is changes to 'gsmvar' and then
#'   it is returned.
#' @export

gmvar_to_gsmvar <- function(gsmvar) {
  stopifnot(class(gsmvar) %in% c("gsmvar", "gmvar"))
  if(class(gsmvar) == "gmvar") {
    class(gsmvar) <- "gsvmar"
    gsmvar$mode$model <- "GMVAR"
  }
  gsmvar
}



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
  gmvar$model$model <- "GMVAR"
  print(gmvar, digits=digits)
}

#' @rdname print.gmvar
#' @export

summary.gmvar <- function(gmvar, ...., digits=2) {
  class(gmvar) <- "gsmvar"
  gmvar$model$model <- "GMVAR"
  summary(gmvar, digits=digits)
}

#' @rdname print.gmvar
#' @export

plot.gmvar <- function(gmvar) {
  class(gmvar) <- "gsmvar"
  gmvar$model$model <- "GMVAR"
  plot(gmvar)
}

# To do: gvmarpred.plot/print, print.gmvarsum,
