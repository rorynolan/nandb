#' Detrend an image series
#'
#' Apply detrending to an image time series using the bleaching correction
#' method described in Digman et al. 2008.
#'
#' @inheritParams Brightness
#'
#' @return The detrended image series.
#' @export
CorrectForBleaching <- function(mat3d, tau) {
  d <- dim(mat3d)
  if (length(d) != 3) stop("mat3d must be a three-dimensional array")
  smoothed <- ExpSmoothPillars(mat3d, tau)
  filtered <- mat3d - smoothed
  means <- MeanPillars(mat3d)
  corrected <- filtered + as.vector(means)  # as it so happens, this will add the means to the pillars as we desire
  corrected[corrected < 0] <- 0
  corrected <- round(corrected)  # return the array back to integer counts
}

#' Find the best tau for exponential filtering detrend.
#'
#' Say you have an image series that you wish to detrend before performing a
#' brightness calculation. This function finds the best \code{tau} for an
#' exponential filtering detrend. See \code{vignette("Adaptive Detrending",
#' package = "nandb")} for more details.
#'
#' @param img.arr A 3-dimensional array of images (where the ith frame is
#'   \code{img.arr[, , i]}).
#' @param mst Do you want to apply an intensity threshold prior to calculating
#'   brightness (via \code{\link{MedStackThresh}})? If so, set your thresholding
#'   \emph{method} here.
#' @param tol What size of error in the estimate of the \emph{ideal} \code{tau}
#'   (aside from the error introduced by the random image simulation, see
#'   \code{vignette("Adaptive Detrending", package = "nandb")}) are you willing
#'   to tolerate?
#'
#' @return A number. The estimate of the ideal \code{tau} to use, with an
#'   attribute "\code{brightness.immobile}" giving the brightness of the
#'   simulated (from all immobile particles) image series after detrending with
#'   this \code{tau} (this should be very close to 1).
#' @export
BestTau <- function(img.arr, mst = NULL, tol = 1) {
  if (!is.null(mst)) img.arr <- MedStackThresh(img.arr, mst, fail = NA)
  raw.brightness.mean <- Brightness(img.arr) %>% mean(na.rm = TRUE)
  if (raw.brightness.mean < 1) {
    stop("Your raw brightness mean is below 1,",
         "there's probably something wrong with your acquisition.")
  }
  means <- apply(img.arr, 3, mean, na.rm = TRUE)
  p <- sum(!is.na(img.arr[, , 1]))
  sim.img.arr <- sapply(means, rpois, n = p)
  BrightnessMeanSimMatTau <- function(sim.img.arr, tau) {
    detrended <- sim.img.arr - ExpSmoothRows(sim.img.arr, tau) +
      rowMeans(sim.img.arr)
    brightnesses <- matrixStats::rowVars(detrended) / rowMeans(detrended)
    mean(brightnesses)
  }
  tau2.sim.brightness.mean <- BrightnessMeanSimMatTau(sim.img.arr, 2)
  if (tau2.sim.brightness.mean > 1) {
    stop("Even with a savage detrend of tau = 2, ",
         "the brightnesses still have mean greater than 1. ",
         "There's probably a problem with your data. ",
         "You should check this out, ",
         "and if you want to work with the data as is, ",
         "then you'll have to choose your own detrend.")
  }
  lower <- 2
  upper <- length(means)
  while (BrightnessMeanSimMatTau(sim.img.arr, upper) < 1) {
    lower <- upper
    upper <- 2 * upper
  }
  TauFarFromOne <- function(tau, sim.img.arr) {
    BrightnessMeanSimMatTau(sim.img.arr, tau) - 1
  }
  root <- uniroot(TauFarFromOne, c(lower, upper), sim.img.arr,
                  tol = tol, extendInt = "upX")
  tau <- root$root
  attr(tau, "brightness.immobile") <- 1 + root$f.root
  tau
}
