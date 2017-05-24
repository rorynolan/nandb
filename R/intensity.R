#' Calculate mean intensity from image series.
#'
#' Given a time stack of images, `MeanIntensity` calculates the mean intensity
#' for each pixel, returning a matrix. `MeanIntensityTxtFolder` does this for
#' every image in a folder, writing the results as text files via
#' [WriteImageTxt()].
#'
#' @param arr An array, can be 3- or 4-dimensional. The first two slots give the
#'   x- and y-cordinates of pixels respectively. If the array is 3-dimensional,
#'   the third slot gives the index of the frame. If it is 4-dimensional, the
#'   third slot indexes the channel and the fourth indexes the frame in the
#'   stack.. To perform this on a file that has not yet been read in, set this
#'   argument to the path to that file (a string).
#' @param mst Do you want to apply an intensity threshold prior to calculating
#'   mean intensities (via [autothresholdr::mean_stack_thresh()])? If so, set your
#'   thresholding \emph{method} here. Pixels failing to exceed the threshold are
#'   set to `NA`.
#' @param skip.consts An image array with only one value (a 'constant array')
#'   won't threshold properly. By default the function would give an error, but
#'   by setting this parameter to `TRUE`, the array would instead be skipped
#'   (the function will return the original array) and give a warning.
#' @param fail If thresholding is done, to which value should pixels not
#'   exceeeding the threshold be set?
#' @param filt Do you want to smooth (`filt = 'smooth'`) or median (`filt =
#'   'median'`) filter the mean intensity image using [SmoothFilterB()] or
#'   [MedianFilterB()] respectively? If selected, these are invoked here with a
#'   filter radius of 1 and with the option `na_count = TRUE`. If you want to
#'   smooth/median filter the mean intensity image in a different way, first
#'   calculate the mean intensities without filtering (`filt = NULL`) using this
#'   function and then perform your desired filtering routine on the result.
#' @param verbose If arr is specified as a file name, print a message to tell
#'   the user that that file is now being processed (useful for
#'   `MeanIntensityTxtFolder`, does not work with multiple cores) and to tell
#'   when `MeanIntensityTxtFolder` is done.
#'
#' @return `MeanIntensity` returns a matrix, the mean-intensity image.
#'   `MeanIntensities` returns a list of these. The result of
#'   `MeanIntensityTxtFolder` is csv files written to disk (in the same folder
#'   as the input images).
#'
#' @examples
#' library(EBImage)
#' library(magrittr)
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' display(normalize(img[, , 1]), method = 'raster')
#' mean.intensity <- MeanIntensity(img, mst = 'Huang', filt = 'median')
#' display(normalize(mean.intensity), method = 'r')
#' two.channel.img <- abind::abind(img, img, along = 4) %>% aperm(c(1, 2, 4, 3))
#' mint.2ch <- MeanIntensity(two.channel.img, mst = "h", filt = "med")
#' @export
MeanIntensity <- function(arr, mst = NULL, filt = NULL, skip.consts = FALSE,
                          verbose = FALSE) {
  if (is.character(arr)) {
    if (verbose)
      message(paste0("Now processing: ", arr, "."))
    arr <- ReadImageData(arr)
  }
  d <- dim(arr)
  if (length(d) == 3) {
    return(MeanIntensity_(arr, mst = mst, filt = filt,
                          skip.consts = skip.consts))
  }
  mint.args <- list(arr3d = ListChannels(arr), mst = mst, filt = filt,
                    skip.consts = skip.consts)
  for (i in seq_along(mint.args)) {
    if (is.null(mint.args[[i]])) {
      mint.args[[i]] <- list(NULL)[rep(1, length(mint.args$arr3d))]
    }
  }
  mint.args %>%
    purrr::pmap(MeanIntensity_) %>%
    ChannelList2Arr
}

MeanIntensity_ <- function(arr3d, mst = NULL,
                           filt = NULL, skip.consts = FALSE) {
  d <- dim(arr3d)
  if (length(d) != 3)
    stop("arr3d must be a three-dimensional array ",
         "or the path to an image on disk which can be read in as a 3d array.")
  if (!is.null(mst)) {
    arr3d <- autothresholdr::mean_stack_thresh(arr3d, method = mst)
  }
  mean.intensity <- MeanPillars(arr3d)
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
#' @param ext The file extension of the images in the folder that you wish to
#'   process. You must wish to process all files with this extension; if there
#'   are files that you don't want to process, take them out of the folder. The
#'   default is for tiff files. Do not use regular expression in this argument.
#' @param mcc The number of parallel cores to use for the processing.
#'
#' @examples
#' setwd(tempdir())
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' WriteIntImage(img, '50.tif')
#' WriteIntImage(img, '50again.tif')
#' MeanIntensityTxtFolder(mcc = 2)
#'
#' @export
MeanIntensityTxtFolder <- function(folder.path = ".", mst = NULL,
  ext = "tif", mcc = parallel::detectCores(), verbose = TRUE) {
  init.dir <- getwd()
  on.exit(setwd(init.dir))
  setwd(folder.path)
  if (filesstrings::StrElem(ext, 1) != ".") ext <- paste0(".", ext)
  ext <- ore::ore_escape(ext) %>% paste0("$")
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
  written <- mapply(WriteImageTxt, mean.intensities, names.noext.mean.intensity)
  if (verbose) ("Done. Please check folder.")
  written
}

#' @rdname MeanIntensity
#'
#' @param arr3d.list A list of 3-dimensional arrays. To perform this on files
#'   that have not yet been read in, set this argument to the path to these
#'   files (a character vector).
#'
#' @examples
#' img.paths <- rep(system.file('extdata', '50.tif', package = 'nandb'), 2)
#' mean.intensities <- MeanIntensities(img.paths, mst = 'Huang', mcc = 2)
#'
#' @export
MeanIntensities <- function(arr3d.list, mst = NULL, skip.consts = FALSE,
  fail = NA, filt = NULL, verbose = FALSE, mcc = parallel::detectCores()) {
  if (is.null(mst)) {
    mean.intensities <- BiocParallel::bplapply(arr3d.list,
      MeanIntensity, filt = filt, verbose = verbose,
      BPPARAM = bpp(mcc))
  } else if (is.list(arr3d.list)) {
    arr3d.list <- lapply(arr3d.list, autothresholdr::mean_stack_thresh,
                         method = mst, fail = fail)
    mean.intensities <- BiocParallel::bplapply(arr3d.list,
      MeanIntensity, filt = filt, verbose = verbose, BPPARAM = bpp(mcc))
  } else {
    if (!is.character(arr3d.list)) {
      stop("arr3d.list must either be a list of 3d arrays, ",
        "or a character vector of paths to the locations ",
        "of 3d arrays on disk.")
    }
    mean.intensities <- list()
    sets <- seq_along(arr3d.list) %>% {
      split(., ((. - 1)%/%mcc) + 1)
    }
    for (i in sets) {
      arrays <- lapply(arr3d.list[i], ReadImageData)
      threshed <- lapply(arrays, autothresholdr::mean_stack_thresh, method = mst,
                         fail = fail)
      mean.intensities.i <- BiocParallel::bplapply(threshed,
        MeanIntensity, filt = filt, verbose = verbose,
        BPPARAM = bpp(mcc))
      mean.intensities[i] <- mean.intensities.i
    }
  }
  mean.intensities <- lapply(mean.intensities, function(x) {
    attr(x, "mst") <- mst
    x
  })
  mean.intensities
}
