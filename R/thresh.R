#' Threshold every image in a stack based on their mean.
#'
#' Via the \code{\link[autothresholdr]{auto_thresh}} function from the
#' \code{autothresholdr} package, apply a threshold to the image stack before
#' calculating brightness, such as to only have nonzero brightnesses in areas
#' where the cell is.
#'
#' This function finds a threshold based on all of the frames and then uses
#' takes the mean of all the frames in the stack image to create a mask and
#' then applies this mask to every frame in the stack (so for a given pillar in
#' the image stack, either it is wiped away to zero, or it is untouched).
#'
#' @param mat3d A 3-dimensional array (the image stack) where the \eqn{n}th
#'   slice is the \eqn{n}th image in the time series.
#' @param method The thresholding method to use. See
#'   \code{\link[autothresholdr]{auto_thresh}}.
#' @param fail To which value should pixels not exceeeding the threshold be set?
#' @param skip.consts An array with only one value (a 'constant array') won't
#'   threshold properly. By default the function would give an error, but by
#'   setting this parameter to \code{TRUE}, the array would instead be skipped
#'   (the function will return the original array) and give a warning.
#'
#' @return The thresholded stack, pillars not exceeding the threshold are set to
#'   zero. The attribute 'threshold' gives the value used for thresholding.
#'
#' @examples
#' library(EBImage)
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' display(normalize(img[, , 1]), method = 'raster')
#' img_thresh_mask <- MeanStackThresh(img, 'Otsu')
#' display(img_thresh_mask[, , 1] > 0, method = 'r')
#' display(normalize(img[, , 1]), method = 'raster')
#' img_thresh_mask <- MeanStackThresh(img, 'Triangle')
#' display(img_thresh_mask[, , 1] > 0, method = 'r')
#'
#' @export
MeanStackThresh <- function(mat3d, method, fail = NA, skip.consts = FALSE) {
  stopifnot(length(dim(mat3d)) == 3)
  if (length(unique(as.vector(mat3d))) == 1 && skip.consts) {
    warning("Constant array, skipping thresholding.")
    return(mat3d)
  }
  thresh <- autothresholdr::auto_thresh(mat3d, method)
  med.stack <- MeanPillars(mat3d)
  med.stack.mask <- med.stack > thresh
  set.indices <- rep(!as.vector(med.stack.mask), dim(mat3d)[3])
  mat3d[set.indices] <- fail
  attr(mat3d, "threshold") <- thresh
  mat3d
}

