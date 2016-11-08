#' Threshold every image in a stack based on their median.
#'
#' Via the \code{\link[autothresholdr]{auto_thresh_mask}} function from the
#' \code{autothresholdr} package, apply a threshold to the image stack before
#' calculating brightness, such as to only have nonzero brightnesses in areas
#' where the cell is.
#'
#' This function takes the median of all the slices in the stack and uses this
#' image to create a mask via \code{\link[autothresholdr]{auto_thresh_mask}} and
#' then applies this mask to every slice in the stack (so for a given pillar in
#' the image stack, either it is wiped away to zero, or it is untouched).
#'
#' @param mat3d A 3-dimensional array (the image stack) where the \eqn{n}th
#'   slice is the \eqn{n}th image in the time series.
#' @param method The thresholding method to use. See
#'   \code{\link[autothresholdr]{auto_thresh_mask}}.
#'
#' @return The thresholded stack, pillars not exceeding the threshold are set to
#'   zero.
#'
#' @examples
#' library(EBImage)
#' img <- ReadImageData(system.file("extdata",
#' "FKBP-mClover_before_0.5nM_AP1510.tif",
#' package = "nandb"))
#' display(normalize(img[, , 1]), method = "raster")
#' img_thresh_mask <- MedStackThresh(img, "Otsu")
#' display(img_thresh_mask[, , 1] > 0, method = "r")
#' display(normalize(img[, , 1]), method = "raster")
#' img_thresh_mask <- MedStackThresh(img, "Triangle")
#' display(img_thresh_mask[, , 1] > 0, method = "r")
#'
#' @export
MedStackThresh <- function(mat3d, method = "Huang") {
  medpillared.rounded <- MedianPillars(mat3d) %>% round
  attr(medpillared.rounded, "bits") <- attr(mat3d, "bits")
  thresh.mask <- autothresholdr::auto_thresh_mask(medpillared.rounded, method)
  mat3d * as.vector(thresh.mask)
}
