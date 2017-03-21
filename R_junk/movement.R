#' Which pixels have too many consecutive low values?
#'
#' If a cell moves into a pixel during the course of a multi-frame acquisition,
#' as the frames go on, the intensity at that pixel will go from low (absence of
#' cell) to high (presence of cell). We detect this movement in a given pixel by
#' detecting the series of consecutive low values. The 'low' values are those
#' failing to exceed a certain threshold (set by thresholding the whole image
#' via [autothresholdr::auto_thresh()]). Of course this also works for
#' cells moving out of pixels.
#'
#' This assumes that you've already thresholded out pixels where the cell
#' \emph{never} was (perhaps via [autothresholdr::auto_thresh()]),
#' setting those to `NA`.
#'
#' @param mat3d A 3d matrix representing the image series.
#' @param thresh The threshold below which values are considered \emph{low}.
#' @param method We want to decide how many consecutive \emph{low} values are
#'   too many. This will be done by using
#'   [autothresholdr::auto_thresh()], so here you set the method for
#'   that function.
#'
#' @return A logical matrix (or vector or array) with `TRUE` for pixels
#'   that have too many consecutive low values and `FALSE` for the others
#'   (so `TRUE` for the \emph{bad} pixels).
CellMoved <- function(mat3d, thresh, method = "Huang") {
  stopifnot(length(dim(mat3d)) == 3)
  mclep <- MostConsecutiveLEsPillars(mat3d, thresh)
  autothresholdr::auto_thresh_mask(max(mclep, na.rm = TRUE) - 
    mclep, method)
}
