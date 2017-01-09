#' @useDynLib nandb
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr "%>%" "%T>%"
NULL

## quiet concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".", "Var1", "Var2", "value", "x", "y", "colour"))
}

.onLoad <- function(libname, pkgname) {
  ## Workaround needed for BiocParallel MulticoreParam to work on mac
  options(bphost = "localhost")
}

#' nandb: Number and brightness in R.
#'
#' The nandb package provides functions for performing number and brightness
#' analysis in R.
#'
#' @docType package
#' @name nandb
#' @references Digman MA, Dalal R, Horwitz AF, Gratton E. Mapping the Number of
#'   Molecules and Brightness in the Laser Scanning Microscope. Biophysical
#'   Journal. 2008;94(6):2320-2332. doi:10.1529/biophysj.107.114645.
NULL
