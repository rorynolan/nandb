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
