#' Calculate mean intensity from image series.
#'
#' Given a time stack of images, \code{MeanIntensity} calculates the mean
#' intensity for each pixel, returning a matrix. \code{MeanIntensityTxtFolder}
#' does this for every image in a folder, writing the results as text files via
#' \code{\link{WriteImageTxt}}.
#'
#' @param mat3d A 3-dimensional array (the image stack) where the \eqn{n}th
#'   slice is the \eqn{n}th image in the time series. To perform this on a file
#'   that has not yet been read in, set this argument to the path to that file
#'   (a string).
#' @param mst Do you want to apply an intensity threshold prior to calculating
#'   mean intensities (via \code{\link{MedStackThresh}})? If so, set your
#'   thresholding \emph{method} here.
#' @param filt Do you want to median/smooth filter (with a radius of 1) the
#'   resulting image using \code{\link{MedianFilterB}} or
#'   \code{\link{SmoothFilterB}}?
#' @param verbose If mat3d is specified as a file name, print a message to tell
#'   the user that that file is now being processed (useful for
#'   \code{MeanIntensityTxtFolder}, does not work with multiple cores) and to
#'   tell when \code{MeanIntensityTxtFolder} is done.
#'
#' @return \code{Brightness} returns a matrix, the brightness image. The result
#'   of \code{BrightnessTxtFolder} is the text csv files written to disk (in the
#'   same folder as the input images).
#'
#' @examples
#' library(EBImage)
#' img <- ReadImageData(system.file("extdata",
#' "low_oligomers.tif",
#' package = "nandb"))
#' display(normalize(img[, , 1]), method = "raster")
#' mean.intensity <- MeanIntensity(img, mst = "Huang", filt = "median")
#' display(normalize(mean.intensity), method = "r")
#' @export
MeanIntensity <- function(mat3d, mst = NULL, filt = NULL, verbose = TRUE) {
  if (is.character(mat3d)) {
    if (verbose) print(paste0("Now processing: ", mat3d, "."))
    mat3d <- ReadImageData(mat3d)
  }
  d <- dim(mat3d)
  if (length(d) != 3) stop("mat3d must be a three-dimensional array")
  if (!is.null(mst)) mat3d <- MedStackThresh(mat3d, method = mst)
  mean.intensity <- MeanPillars(mat3d)
  if (!is.null(filt)) {
    allowed <- c("median", "smooth")
    filt <- tolower(filt)
    sw <- startsWith(allowed, filt)
    if (!any(sw)) stop("filt must be either 'median' or 'smooth'")
    filt <- allowed[sw]
    if (filt == "median") {
      mean.intensity <- MedianFilterB(mean.intensity, na_count = TRUE)
    } else {
      mean.intensity <- SmoothFilterB(mean.intensity, na_count = TRUE)
    }
  }
  attr(mean.intensity, "frames") <- d[3]
  attr(mean.intensity, "thresh") <- ifelse(is.null(mst), NA, mst)
  attr(mean.intensity, "filter") <- ifelse(is.null(filt), NA, filt)
  mean.intensity
}

#' @rdname MeanIntensity
#'
#' @param folder.path The path (relative or absolute) to the folder you wish to
#'   process.
#' @param ext the file extension of the images in the folder that you wish to
#'   process (can be rooted in regular expression for extra-safety, as in the
#'   default). You must wish to process all files with this extension; if there
#'   are files that you don't want to process, take them out of the folder.
#' @param mcc The number of parallel cores to use for the processing.
#'
#' @export
MeanIntensityTxtFolder <- function(folder.path = ".", tau = NA,
                                mst = NULL, ext = "\\.tif$",
                                mcc = parallel::detectCores(),
                                verbose = TRUE) {
  init.dir <- getwd()
  on.exit(setwd(init.dir))
  setwd(folder.path)
  file.names <- list.files(pattern = ext)
  mean.intensities <- MCLapply(file.names, MeanIntensity, mst = mst,
                           mcc = mcc, verbose = verbose)
  frames <- sapply(mean.intensities, function(x) attr(x, "frames"))
  msts <- sapply(mean.intensities, function(x) attr(x, "thresh"))
  filters <- sapply(mean.intensities, function(x) attr(x, "filter"))
  names.noext.brightness <- sapply(file.names, filesstrings::StrBeforeNth,
                                   stringr::coll("."), -1) %>%
    paste0("_MeanIntensity_frames=", frames, "_mst=", msts, "_filter=", filters)
  mapply(WriteImageTxt, mean.intensities, names.noext.brightness) %>% invisible
  if (verbose) print("Done. Please check folder.")
}
