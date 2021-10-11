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
#' @param x a class 'gmvar' object. THIS CLASS IS DEPRECATED FROM THE VERSION
#'   2.0.0 ONWARDS.
#' @param object  a class 'gmvar' object. THIS CLASS IS DEPRECATED FROM THE VERSION
#'   2.0.0 ONWARDS.
#' @param digits number of digits to be printed.
#' @param ... See the usage from the documentation of the appropriate class 'gsmvar' S3 method.
#' @details These methods exist so that models estimated with earlier versions
#'   of the package can be used normally.
#' @export

print.gmvar <- function(x, ..., digits=2) {
  gsmvar <- gmvar_to_gsmvar(x)
  print(gsmvar, ..., digits=digits)
}


#' @rdname print.gmvar
#' @export

summary.gmvar <- function(object, ..., digits) {
  gsmvar <- gmvar_to_gsmvar(object)
  summary(gsmvar, ..., digits=digits)
}

#' @rdname print.gmvar
#' @export

plot.gmvar <- function(x, ...) {
  gsmvar <- gmvar_to_gsmvar(x)
  plot(gsmvar, ...)
}


#' @rdname print.gmvar
#' @param object object of class \code{'gmvar'}. THIS CLASS IS DEPRECATED FROM THE VERSION
#'   2.0.0 ONWARDS.
#' @export

logLik.gmvar <- function(object, ...) {
  gsmvar <- gmvar_to_gsmvar(object)
  gsmvar$loglik
}


#' @rdname print.gmvar
#' @export

residuals.gmvar <- function(object, ...) {
  gsmvar <- gmvar_to_gsmvar(object)
  res <- gsmvar$quantile_residuals
  colnames(res) <- colnames(object$data)
  res
}


#' @title plot method for class 'gmvarpred' objects
#'
#' @description \code{plot.gmvarpred} is plot method for gsmvarpred objects.
#'  EXISTS FOR BACKWARD COMPATIBILITY. THE CLASS 'gmvarpred' IS DEPRECATED FROM
#'  THE VERSION 2.0.0 ONWARD: WE USE THE CLASS 'gsmvarpred' NOW.
#'
#' @inheritParams plot.gsmvarpred
#' @param digits how many digits to print?
#' @details These methods exist so that objects created with earlier versions
#'   of the package can be used normally.
#' @export

plot.gmvarpred <- function(x, ..., nt, mix_weights=TRUE, add_grid=TRUE) {
  class(x) <- "gsmvarped"
  plot(x, ..., mix_weights=mix_weights, add_grid=add_grid)
}


#' @rdname plot.gmvarpred
#' @export

print.gmvarpred <- function(x, ..., digits=2) {
  class(x) <- "gsmvarped"
  print(x, ..., digits=digits)
}



#' @title Summary print method from objects of class 'gmvarsum'
#'
#' @description \code{print.gmvarsum} is a print method for object \code{'gmvarsum'}.
#'   EXISTS FOR BACKWARD COMPATIBILITY. CLASS 'gmvarsum' IS DEPRECATED FROM THE VERSION
#'   2.0.0. ONWARDS. NOW, WE USE THE CLASS 'gsmvarsum'.
#' @inheritParams print.gsmvarsum
#' @export

print.gmvarsum <- function(x, ..., digits) {
  class(x) <- "gsmvarsum"
  print(x, ..., digits=2)
}

