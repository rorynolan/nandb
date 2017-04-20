#' Calculate brightness from image series.
#'
#' Given a time stack of images, `Brightness` performs a calculation of the
#' brightness for each pixel. `BrightnessTxtFolder` does this calculation
#' for an entire folder, writing the results as text files via
#' [WriteImageTxt()]. `Brightnesses` calculates the brightnesses
#' for several image series in parallel.
#'
#' Do not try to parallelize the use of `Brightness` and friends yourself
#' (e.g. with [mclapply()]) because it will throw an error (this is
#' because the `autothresholdr` package does not work in parallel (indeed
#' anything run using the `rJava` package won't)). Always use
#' `Brightnesses` for this purpose (this has a workaround whereby it does
#' the thresholding in series and does the rest in parallel).
#'
#' @param arr An array, can be 3- or 4-dimensional. The first two slots give the
#'   x- and y-cordinates of pixels respectively. If the array is 3-dimensional,
#'   the third slot gives the index of the frame. If it is 4-dimensional, the
#'   third slot indexes the channel and the fourth indexes the frame in the
#'   stack. To perform this on a file that has not yet been read in, set this
#'   argument to the path to that file (a string).
#' @param tau If this is specified, bleaching correction is performed with
#'   [CorrectForBleaching()] which uses exponential filtering with
#'   time constant `tau` (where the unit of time is the time between
#'   frames). If this is set to `'auto'`, then the value of `tau` is
#'   calculated automatically via [BestTau()].
#' @param mst Do you want to apply an intensity threshold prior to calculating
#'   brightness (via [MeanStackThresh()])? If so, set your
#'   thresholding \emph{method} here.
#' @param skip.consts An image array with only one value (a 'constant array')
#'   won't threshold properly. By default the function would give an error, but
#'   by setting this parameter to `TRUE`, the array would instead be
#'   skipped (the thresholding will return the original array) and give a
#'   warning.
#' @param filt Do you want to smooth (`filt = 'smooth'`) or median
#'   (`filt = 'median'`) filter the brightness image using
#'   [SmoothFilterB()] or [MedianFilterB()] respectively? If
#'   selected, these are invoked here with a filter radius of 1 and with the
#'   option `na_count = TRUE`. If you want to smooth/median filter the
#'   brightness image in a different way, first calculate the brightnesses
#'   without filtering (`filt = NULL`) using this function and then perform
#'   your desired filtering routine on the result.
#' @param verbose If arr3d is specified as a file name, print a message to tell
#'   the user that that file is now being processed (useful for
#'   `BrightnessTxtFolder`, does not work with multiple cores).
#'
#' @return `Brightness` returns a matrix, the brightness image;
#'   `Brightnesses` returns a list of these. The result of
#'   `BrightnessTxtFolder` is the text csv files written to disk (in the
#'   same folder as the input images).
#'
#' @examples
#' library(magrittr)
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' EBImage::display(EBImage::normalize(img[, , 1]), method = 'raster')
#' brightness <- Brightness(img, tau = 'auto', mst = 'Huang', filt = 'median')
#' MatrixRasterPlot(brightness, log.trans = TRUE)
#' two.channel.img <- abind::abind(img, img, along = 4) %>% aperm(c(1, 2, 4, 3))
#' brightness.2ch <- Brightness(two.channel.img)
#' @export
Brightness <- function(arr, tau = NA, mst = NULL, skip.consts = FALSE,
                       filt = NULL, verbose = FALSE) {
  if (is.character(arr)) {
    if (verbose)
      message(paste0("Now processing: ", arr, "."))
    arr <- ReadImageData(arr)
  }
  d <- dim(arr)
  if (length(d) == 3) {
    return(Brightness_(arr, tau = tau, mst = mst,
                       skip.consts = skip.consts, filt = filt))
  }
  brightness.args <- list(arr3d = ListChannels(arr), tau = tau, mst = mst,
                          skip.consts = skip.consts, filt = filt)
  for (i in seq_along(brightness.args)) {
    if (is.null(brightness.args[[i]])) {
      brightness.args[[i]] <- list(NULL)[rep(1, length(brightness.args$arr3d))]
    }
  }
  brightness.args %>%
    purrr::pmap(Brightness_) %>%
    ChannelList2Arr
}

Brightness_ <- function(arr3d, tau = NA, mst = NULL, skip.consts = FALSE,
                        filt = NULL) {
  d <- dim(arr3d)
  if (length(d) != 3)
    stop("arr3d must be a three-dimensional array")
  if (!is.null(mst)) {
    arr3d <- MeanStackThresh(arr3d, method = mst, fail = NA,
      skip.consts = skip.consts)
  }
  tau.auto <- FALSE
  if (!is.na(tau)) {
    if (is.character(tau)) {
      tau <- tolower(tau)
      if (startsWith("auto", tau)) {
        tau <- BestTau(arr3d)
        tau.auto <- TRUE
      } else {
        stop("If tau is a string, it must be 'auto'.")
      }
    } else if ((!is.numeric(tau)) && (!is.na(tau))) {
      stop("If tau is not numeric, then it must be NA or 'auto'.")
    }
    arr3d <- CorrectForBleaching(arr3d, tau)
  }
  brightness <- VarPillars(arr3d)/MeanPillars(arr3d)
  if (!is.null(filt)) {
    allowed <- c("median", "smooth")
    filt <- tolower(filt)
    sw <- startsWith(allowed, filt)
    if (!any(sw))
      stop("filt must be either 'median' or 'smooth'")
    filt <- allowed[sw]
    if (filt == "median") {
      brightness <- MedianFilterB(brightness, na_count = TRUE)
    } else {
      brightness <- SmoothFilterB(brightness, na_count = TRUE)
    }
  }
  tau <- ifelse(tau.auto, ifelse(is.na(tau), "auto=NA", stringr::str_c("auto=",
    round(tau))), tau)
  filter <- ifelse(is.null(filt), NA, filt)
  mst <- ifelse(is.null(mst), NA, mst)
  new.brightness.attrs <- list(frames = d[3], tau = tau, filter = filter,
    mst = mst)
  attributes(brightness) <- c(attributes(brightness), new.brightness.attrs)
  brightness
}

#' Create a Brightness time-series.
#'
#' Given a stack of images `arr`, use the first `frames.per.set` of
#' them to create one brightness image, the next `frames.per.set` of them
#' to create the next brightness image and so on to get a time-series of
#' brightness images. If `tau` is specified, bleaching correction is
#' performed via [CorrectForBleaching()].
#'
#' This may discard some images, for example if 175 frames are in the input and
#' `frames.per.set = 50`, then the last 25 are discarded. If bleaching
#' correction is selected, it is performed on the whole image stack before the
#' sectioning is done for calculation of brightnesses.
#'
#' @inheritParams Brightness
#' @param frames.per.set The number of frames with which to calculate the
#'   successive brightnesses.
#' @param mcc The number of cores to use for the parallel processing.
#' @param seed A seed for the random number generation for [BestTau]. Don't use
#'   [set.seed], it won't work.
#'
#' @return An array where the \eqn{i}th slice is the \eqn{i}th brightness image.
#' @seealso [Brightness()].
#'
#' @examples
#' library(magrittr)
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' EBImage::display(EBImage::normalize(img[, , 1]), method = 'raster')
#' bts <- BrightnessTimeSeries(img, 10, tau = 'auto', mst = 'tri',
#' filt = 'median', mcc = 2)
#' two.channel.img <- abind::abind(img, img, along = 4) %>% aperm(c(1, 2, 4, 3))
#' bts.2ch <- BrightnessTimeSeries(two.channel.img, 10)
#' @export
BrightnessTimeSeries <- function(arr, frames.per.set, tau = NA,
                                 mst = NULL, skip.consts = FALSE, filt = NULL,
                                 verbose = FALSE,
                                 mcc = parallel::detectCores(), seed = NULL) {
  if (is.character(arr)) {
    if (verbose)
      message(paste0("Now processing: ", arr, "."))
    arr <- ReadImageData(arr)
  }
  d <- dim(arr)
  if (length(d) == 3) {
    return(BrightnessTimeSeries_(arr, frames.per.set = frames.per.set,
                                 tau = tau, mst = mst,
                                 skip.consts = skip.consts, filt = filt,
                                 verbose = verbose, mcc = mcc, seed = seed))
  }
  bts.args <- list(arr3d = ListChannels(arr), frames.per.set = frames.per.set,
                   tau = tau, mst = mst, skip.consts = skip.consts, filt = filt,
                   verbose = verbose, mcc = mcc, seed = seed)
  for (i in seq_along(bts.args)) {
    if (is.null(bts.args[[i]])) {
      bts.args[[i]] <- list(NULL)[rep(1, length(bts.args$arr3d))]
    }
  }
  bts.args %>%
    purrr::pmap(BrightnessTimeSeries_) %>%
    ChannelList2Arr
}
BrightnessTimeSeries_ <- function(arr3d, frames.per.set, tau = NA,
  mst = NULL, skip.consts = FALSE, filt = NULL, verbose = FALSE,
  mcc = parallel::detectCores(), seed = NULL) {
  d <- dim(arr3d)
  if (length(d) != 3)
    stop("arr3d must be a three-dimensional array")
  if (frames.per.set > d[3]) {
    stop("frames.per.set must not be greater than the depth of arr3d")
  }
  if (!is.null(mst)) arr3d <- MeanStackThresh(arr3d, mst)
  if (!is.na(tau)) arr3d <- CorrectForBleaching(arr3d, tau)
  n.sets <- d[3] %/% frames.per.set
  set.indices <- lapply(seq_len(n.sets), function(x) {
    ((x - 1) * frames.per.set + 1):(x * frames.per.set)
  })
  sets <- lapply(set.indices, Slices, arr3d)
  brightnesses <- Brightnesses(sets, tau = NA, mst = NULL,
    skip.consts = skip.consts, filt = filt, mcc = mcc, seed = seed) %>%
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
#' @param ext The file extension of the images in the folder that you wish to
#'   process. You must wish to process all files with this extension; if there
#'   are files that you don't want to process, take them out of the folder. The
#'   default is for tiff files. Do not use regular expression in this argument.
#' @param mcc The number of cores to use for the parallel processing.
#' @param seed A seed for the random number generation for [BestTau]. Don't use
#'   [set.seed], it won't work.
#'
#' @examples
#' setwd(tempdir())
#' WriteIntImage(img, '50.tif')
#' WriteIntImage(img, '50again.tif')
#' BrightnessTxtFolder(tau = 'auto', mst = 'tri', mcc = 2)
#' list.files()
#' file.remove(list.files())  # cleanup
#' @export
BrightnessTxtFolder <- function(folder.path = ".", tau = NA,
  mst = NULL, skip.consts = FALSE, filt = NULL, ext = "tif",
  mcc = parallel::detectCores(), verbose = FALSE, seed = NULL) {
  init.dir <- getwd()
  on.exit(setwd(init.dir))
  setwd(folder.path)
  if (filesstrings::StrElem(ext, 1) != ".") ext <- paste0(".", ext)
  ext <- ore::ore_escape(ext) %>% paste0("$")
  file.names <- list.files(pattern = ext)
  brightnesses <- Brightnesses(file.names, tau = tau, mst = mst,
    skip.consts = skip.consts, filt = filt, seed = seed)
  frames <- vapply(brightnesses, function(x) attr(x, "frames"), integer(1))
  if (skip.consts) {
    # this may seem inefficient but it's cost is negligible
    # relative to that of the brightness calculations
    CheckConst <- function(file.path) {
      arr <- ReadImageData(file.path)
      length(unique(as.vector(arr))) == 1
    }
    const <- vapply(file.names, CheckConst, TRUE)
    mst <- ifelse(const, "fail", mst)
  }
  tau <- simplify2array(lapply(brightnesses, function(x) attr(x, "tau")))
  names.noext.brightness <- vapply(file.names, filesstrings::BeforeLastDot,
                                   character(1)) %>%
    paste0("_brightness_frames=", frames, "_tau=", tau, "_mst=",
      mst, "_filter=", ifelse(is.null(filt), NA, filt))
  mapply(WriteImageTxt, brightnesses, names.noext.brightness) %>%
    invisible
  if (verbose)
    message("Done. Please check folder.")
}

#' @rdname Brightness
#'
#' @param arr3d.list A list of 3-dimensional arrays. To perform this on files
#'   that have not yet been read in, set this argument to the path to these
#'   files (a character vector).
#'
#' @examples
#' img.paths <- rep(system.file('extdata', '50.tif', package = 'nandb'), 2)
#' brightnesses <- Brightnesses(img.paths, mst = 'Huang', tau = 'auto', mcc = 2)
#'
#' @export
Brightnesses <- function(arr3d.list, tau = NA, mst = NULL, skip.consts = FALSE,
  filt = NULL, verbose = FALSE, mcc = parallel::detectCores(), seed = NULL) {
  if (is.null(mst)) {
    brightnesses <- BiocParallel::bplapply(arr3d.list, Brightness, tau = tau,
      filt = filt, verbose = verbose, BPPARAM = bpp(mcc, seed = seed))
  } else if (is.list(arr3d.list)) {
    arr3d.list <- lapply(arr3d.list, MeanStackThresh, method = mst,
      fail = NA, skip.consts = skip.consts)
    brightnesses <- BiocParallel::bplapply(arr3d.list, Brightness, tau = tau,
      filt = filt, verbose = verbose, BPPARAM = bpp(mcc, seed = seed))
  } else {
    if (!is.character(arr3d.list)) {
      stop("arr3d.list must either be a list of 3d arrays, ",
        "or a character vector of paths to the locations ",
        "of 3d arrays on disk.")
    }
    brightnesses <- list()
    sets <- seq_along(arr3d.list) %>% {
      split(., ((. - 1) %/% mcc) + 1)
    }
    for (i in sets) {
      arrays <- lapply(arr3d.list[i], ReadImageData)
      threshed <- lapply(arrays, MeanStackThresh, method = mst,
        fail = NA, skip.consts = skip.consts)
      brightnesses.i <- BiocParallel::bplapply(threshed,
        Brightness, tau = tau, filt = filt, verbose = verbose,
        BPPARAM = bpp(mcc, seed = seed))
      brightnesses[i] <- brightnesses.i
    }
  }
  brightnesses <- lapply(brightnesses, function(x) {
    attr(x, "mst") <- mst
    x
  })
  brightnesses
}
