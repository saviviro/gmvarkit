#' @title Makes the old class 'gmvar' objects compatible with the functions using class 'gsmvar' objects
#'
#' @description \code{gmvar_to_gsmvar} makes class 'gmvar' objects compatible with the functions using
#' class 'gsmvar' objects
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
  stopifnot(methods::is(gsmvar) %in% c("gsmvar", "gmvar"))
  if(inherits(gsmvar, what="gmvar")) {
    class(gsmvar) <- "gsmvar"
    gsmvar$model$model <- "GMVAR"
  }
  gsmvar
}
