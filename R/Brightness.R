#' Calculate brightness from image series.
#'
#' Given a time stack of images, this function performs a calculation of the
#' molecular brightness for each pixel. If \code{tau} is specified, bleaching
#' correction is performed via \code{\link{CorrectForBleaching}}.
#'
#' @param mat3d A 3-dimensional array (the image stack) where the \eqn{n}th
#'   slice is the \eqn{n}th image in the time series.
#' @param tau The time constant for the exponential filtering.
#'
#' @return A matrix, the brightness image.
#' @export
Brightness <- function(mat3d, tau = NA) {
  d <- dim(mat3d)
  if (length(d) != 3) stop("mat3d must be a three-dimensional array")
  if (!is.na(tau)) mat3d <- CorrectForBleaching(mat3d, tau)
  brightness <- VarPillars(mat3d) / MeanPillars(mat3d)
  brightness
}

#' Create a Brightness time-series.
#'
#' Given a stack of images, use the first \code{frames.per.set} of them to
#' create one brightness image, the next \code{frames.per.set} of them to create
#' the next brightness image and so on to get a time-series of brightness
#' images. If \code{tau} is specified, bleaching correction is performed via
#' \code{\link{CorrectForBleaching}}.
#'
#' This may discard some images, for example if 175 frames are in the input and
#' \code{frames.per.set = 50}, then the last 25 are discarded. If bleaching
#' correction is selected, it is performed on the whole image stack before the
#' sectioning is done for calculation of brightnesses.
#'
#' @inheritParams Brightness
#' @param frames.per.set The number of frames with which to calculate the
#'   successive brightnesses.
#'
#' @return An array where the \eqn{i}th slice is the \eqn{i}th brightness image.
#' @seealso \code{\link{Brightness}}.
BrightnessTimeSeries <- function(mat3d, frames.per.set, tau = NA) {
  d <- dim(mat3d)
  if (length(d) != 3) stop("mat3d must be a three-dimensional array")
  if (frames.per.set > d[3]) stop("frames.per.set must not be greater than the depth of mat3d")
  if (!is.na(tau)) mat3d <- CorrectForBleaching(mat3d, tau)
  n.sets <- floor(d[3] / frames.per.set)
  set.indices <- lapply(seq_len(n.sets),
                        function(x) ((x - 1) * frames.per.set + 1):(x * frames.per.set))
  sets <- lapply(set.indices, Slices, detrended)
  brightnesses <- lapply(sets, Brightness) %>%
    Reduce(function(x, y) abind::abind(x, y, along = 3), .)
  if (length(dim(brightnesses)) == 2) {
    brightnesses <- abind::abind(brightnesses, along = 3)
  }
  brightnesses
}
