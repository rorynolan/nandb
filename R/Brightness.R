#' Calculate brightness from image series.
#'
#' Given a time stack of images, this \code{Brightness} performs a calculation
#' of the molecular brightness for each pixel. If \code{tau} is specified,
#' bleaching correction is performed via \code{\link{CorrectForBleaching}}.
#' \code{BrightnessTxtFolder} does this calculation for an entire folder,
#' writing the results as text files via \code{\link{WriteImageTxt}}.
#'
#' @param mat3d A 3-dimensional array (the image stack) where the \eqn{n}th
#'   slice is the \eqn{n}th image in the time series. To perform this on a file
#'   that has not yet been read in, set this argument to the path to that file
#'   (a string).
#' @param tau The time constant for the exponential filtering.
#' @param mst Do you want to apply an intensity threshold prior to calculating
#'   brightness (via \code{\link{MedStackThresh}})? If so, set your thresholding
#'   \emph{method} here.
#' @param filt Do you want to median/smooth filter (with a radius of 1) the
#'   brightness image using \code{\link{MedianFilterB}} or
#'   \code{\link{SmoothFilterB}}? If \code{filt = "median"}, \code{\link{MedianFilterB}} is invoked with the option \code{na_count = TRUE}.
#' @param verbose If mat3d is specified as a file name, print a message to tell
#'   the user that that file is now being processed (useful for
#'   \code{BrightnessFolder}, does not work with multiple cores) and to tell
#'   when \code{MeanIntensityTxtFolder} is done.
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
#' brightness <- Brightness(img, tau = 10, mst = "Huang", filt = "median")
#' KmerPlot(brightness, 1.16)
#' @export
Brightness <- function(mat3d, tau = NA, mst = NULL, filt = NULL, verbose = TRUE) {
  if (is.character(mat3d)) {
    if (verbose) print(paste0("Now processing: ", mat3d, "."))
    mat3d <- ReadImageData(mat3d)
  }
  d <- dim(mat3d)
  if (length(d) != 3) stop("mat3d must be a three-dimensional array")
  if (!is.null(mst)) mat3d <- MedStackThresh(mat3d, method = mst)
  if (!is.na(tau)) mat3d <- CorrectForBleaching(mat3d, tau)
  brightness <- VarPillars(mat3d) / MeanPillars(mat3d)
  if (!is.null(filt)) {
    allowed <- c("median", "smooth")
    filt <- tolower(filt)
    sw <- startsWith(allowed, filt)
    if (!any(sw)) stop("filt must be either 'median' or 'smooth'")
    filt <- allowed[sw]
    if (filt == "median") {
      brightness <- MedianFilterB(brightness, na_count = TRUE)
    } else {
      brightness <- SmoothFilterB(brightness, na_count = TRUE)
    }
  }
  attributes(brightness) <- c(attributes(brightness),
                              list(frames = d[3], tau = tau,
                                 filter = ifelse(is.null(filt), NA, filt),
                                 mst = ifelse(is.null(mst), NA, mst)
                                 )
  )
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
BrightnessTimeSeries <- function(mat3d, frames.per.set, pbt = NULL, tau = NA) {
  d <- dim(mat3d)
  if (length(d) != 3) stop("mat3d must be a three-dimensional array")
  if (frames.per.set > d[3]) stop("frames.per.set must not be greater than the depth of mat3d")
  if (!is.na(tau)) mat3d <- CorrectForBleaching(mat3d, tau)
  n.sets <- floor(d[3] / frames.per.set)
  set.indices <- lapply(seq_len(n.sets),
                        function(x) ((x - 1) * frames.per.set + 1):(x * frames.per.set))
  sets <- lapply(set.indices, Slices, detrended)
  brightnesses <- lapply(sets, Brightness, tau = tau, pbt = pbt) %>%
    Reduce(function(x, y) abind::abind(x, y, along = 3), .)
  if (length(dim(brightnesses)) == 2) {
    brightnesses <- abind::abind(brightnesses, along = 3)
  }
  brightnesses
}

#' @rdname Brightness
#'
#' @param folder.path The path (relative or absolute) to the folder you wish to
#'   process.
#' @param ext the file extension of the images in the folder that you wish to
#'   process (can be rooted in regular expression for extra-safety, as in the
#'   default). You must wish to process all files with this extension; if there
#'   are files that you don't want to process, take them out of the folder.
#'
#' @export
BrightnessTxtFolder <- function(folder.path = ".", tau = NA, mst = NULL,
                                filt = NULL, ext = "\\.tif$",
                                mcc = parallel::detectCores(), verbose = TRUE) {
  init.dir <- getwd()
  on.exit(setwd(init.dir))
  setwd(folder.path)
  file.names <- list.files(pattern = ext)
  brightnesses <- MCLapply(file.names, Brightness, mst = mst, tau = tau,
                           filt = filt, mcc = mcc, verbose = verbose)
  frames <- sapply(brightnesses, function(x) attr(x, "frames"))
  names.noext.brightness <- sapply(file.names, filesstrings::StrBeforeNth,
                                   stringr::coll("."), -1) %>%
    paste0("_brightness_frames=", frames, "_tau=", tau,
           "_mst=", mst, "_filter=", filt)
  mapply(WriteImageTxt, brightnesses, names.noext.brightness) %>% invisible
  if (verbose) print("Done. Please check folder.")
}
