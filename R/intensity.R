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
#'   mean intensities (via \code{\link{MeanStackThresh}})? If so, set your
#'   thresholding \emph{method} here. Pixels failing to exceed the threshold are
#'   set to \code{NA}.
#' @param skip.consts An image array with only one value (a 'constant array')
#'   won't threshold properly. By default the function would give an error, but
#'   by setting this parameter to \code{TRUE}, the array would instead be
#'   skipped (the function will return the original array) and give a warning.
#' @param fail If thresholding is done, to which value should pixels not
#'   exceeeding the threshold be set?
#' @param filt Do you want to smooth (\code{filt = 'smooth'}) or median
#'   (\code{filt = 'median'}) filter the mean intensity image using
#'   \code{\link{SmoothFilterB}} or \code{\link{MedianFilterB}} respectively? If
#'   selected, these are invoked here with a filter radius of 1 and with the
#'   option \code{na_count = TRUE}. If you want to smooth/median filter the mean
#'   intensity image in a different way, first calculate the mean intensities
#'   without filtering (\code{filt = NULL}) using this function and then perform
#'   your desired filtering routine on the result.
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
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' display(normalize(img[, , 1]), method = 'raster')
#' mean.intensity <- MeanIntensity(img, mst = 'Huang', filt = 'median')
#' display(normalize(mean.intensity), method = 'r')
#' @export
MeanIntensity <- function(mat3d, mst = NULL, filt = NULL, skip.consts = FALSE,
  fail = NA, verbose = FALSE) {
  if (is.character(mat3d)) {
    if (verbose)
      print(paste0("Now processing: ", mat3d, "."))
    mat3d <- ReadImageData(mat3d)
  }
  d <- dim(mat3d)
  if (length(d) != 3)
    stop("mat3d must be a three-dimensional array")
  if (!is.null(mst)) {
    mat3d <- MeanStackThresh(mat3d, method = mst, fail = fail,
      skip.consts = skip.consts)
  }
  mean.intensity <- MeanPillars(mat3d)
  if (!is.null(filt)) {
    allowed <- c("median", "smooth")
    filt <- tolower(filt)
    sw <- startsWith(allowed, filt)
    if (!any(sw))
      stop("filt must be either 'median' or 'smooth'")
    filt <- allowed[sw]
    if (filt == "median") {
      mean.intensity <- MedianFilterB(mean.intensity, na_count = TRUE)
    } else {
      mean.intensity <- SmoothFilterB(mean.intensity, na_count = TRUE)
    }
  }
  attr(mean.intensity, "frames") <- d[3]
  attr(mean.intensity, "thresh") <- ifelse(is.null(mst), NA_character_, mst)
  attr(mean.intensity, "filter") <- ifelse(is.null(filt), NA_character_, filt)
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
#' @examples
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' dir.create('tempdir')
#' WriteIntImage(img, 'tempdir/50.tif')
#' WriteIntImage(img, 'tempdir/50again.tif')
#' MeanIntensityTxtFolder('tempdir', mcc = 2)
#' filesstrings::RemoveDirs('tempdir')
#'
#' @export
MeanIntensityTxtFolder <- function(folder.path = ".", mst = NULL,
  ext = "\\.tif$", mcc = parallel::detectCores(), verbose = TRUE) {
  init.dir <- getwd()
  on.exit(setwd(init.dir))
  setwd(folder.path)
  file.names <- list.files(pattern = ext)
  mean.intensities <- MeanIntensities(file.names, MeanIntensity,
    mst = mst, verbose = verbose, mcc = mcc)
  frames <- vapply(mean.intensities, function(x) attr(x, "frames"), integer(1))
  msts <- vapply(mean.intensities, function(x) attr(x, "thresh"), character(1))
  filters <- vapply(mean.intensities, function(x) attr(x, "filter"),
                    character(1))
  names.noext.mean.intensity <- vapply(file.names,
                                  filesstrings::BeforeLastDot, character(1)) %>%
    paste0("_MeanIntensity_frames=", frames, "_mst=", msts, "_filter=", filters)
  invisible(mapply(WriteImageTxt, mean.intensities, names.noext.mean.intensity))
  if (verbose)
    print("Done. Please check folder.")
}

#' @rdname MeanIntensity
#'
#' @param mat3d.list A list of 3-dimensional arrays. To perform this on files
#'   that have not yet been read in, set this argument to the path to these
#'   files (a character vector).
#'
#' @examples
#' img.paths <- rep(system.file('extdata', '50.tif', package = 'nandb'), 2)
#' mean.intensities <- MeanIntensities(img.paths, mst = 'Huang', mcc = 2)
#'
#' @export
MeanIntensities <- function(mat3d.list, mst = NULL, skip.consts = FALSE,
  fail = NA, filt = NULL, verbose = FALSE, mcc = parallel::detectCores()) {
  if (is.null(mst)) {
    mean.intensities <- BiocParallel::bplapply(mat3d.list,
      MeanIntensity, filt = filt, verbose = verbose,
      BPPARAM = suppressWarnings(BiocParallel::MulticoreParam(workers = mcc)))
  } else if (is.list(mat3d.list)) {
    mat3d.list <- lapply(mat3d.list, MeanStackThresh, method = mst,
      fail = fail, skip.consts = skip.consts)
    mean.intensities <- BiocParallel::bplapply(mat3d.list,
      MeanIntensity, filt = filt, verbose = verbose,
      BPPARAM = suppressWarnings(BiocParallel::MulticoreParam(workers = mcc)))
  } else {
    if (!is.character(mat3d.list)) {
      stop("mat3d.list must either be a list of 3d arrays, ",
        "or a character vector of paths to the locations ",
        "of 3d arrays on disk.")
    }
    mean.intensities <- list()
    sets <- seq_along(mat3d.list) %>% {
      split(., ((. - 1)%/%mcc) + 1)
    }
    for (i in sets) {
      arrays <- lapply(mat3d.list[i], ReadImageData)
      threshed <- lapply(arrays, MeanStackThresh, method = mst,
        fail = fail, skip.consts = skip.consts)
      mean.intensities.i <- BiocParallel::bplapply(threshed,
        MeanIntensity, filt = filt, verbose = verbose,
        BPPARAM = suppressWarnings(BiocParallel::MulticoreParam(workers = mcc)))
      mean.intensities[i] <- mean.intensities.i
    }
  }
  mean.intensities <- lapply(mean.intensities, function(x) {
    attr(x, "mst") <- mst
    x
  })
  mean.intensities
}
