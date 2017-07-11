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
#' @param n.ch The number of channels in the image (default 1).
#' @param tau If this is specified, bleaching correction is performed with
#'   [CorrectForBleaching()] which uses exponential filtering with
#'   time constant `tau` (where the unit of time is the time between
#'   frames). If this is set to `'auto'`, then the value of `tau` is
#'   calculated automatically via [BestTau()].
#' @param mst Do you want to apply an intensity threshold prior to calculating
#'   brightness (via [autothresholdr::mean_stack_thresh()])? If so, set your
#'   thresholding \emph{method} here.
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
#' brightness <- Brightness(img, tau = "auto", mst = 'Huang', filt = 'median')
#' MatrixRasterPlot(brightness, log.trans = TRUE)
#' two.channel.img <- abind::abind(img, img, along = 4) %>% aperm(c(1, 2, 4, 3))
#' brightness.2ch <- Brightness(two.channel.img, n.ch = 2)
#' @export
Brightness <- function(arr, tau = NA, mst = NULL,
                       filt = NULL, n.ch = 1, verbose = FALSE) {
  if (is.character(arr)) {
    if (verbose)
      message(paste0("Now processing: ", arr, "."))
    arr <- ReadImageData(arr)
  }
  d <- dim(arr)
  if (n.ch == 1) {
    if (length(d) != 3) {
      stop("Expected a three-dimensional image array image but got a ",
           length(d), " dimensional one.")
    }
  }
  if (n.ch > 1) {
    ld <- length(d)
    if (! ld %in% c(3, 4)) {
      stop("There is a problem with your image. It was read in as a ", ld,
           "-dimensional image. If read in as 4-dimensional, it is left alone, ",
           "or if read in as 3-dimensional, it is coerced to 4-dimansional ",
           "via nandb::ForceChannels", "but for ", ld, "-dimensional image, ",
           "there's nothing to be done.")
    }
    if (ld == 3) arr <- ForceChannels(arr, n.ch)
    d <- dim(arr)
  }
  if (length(d) == 3) {
    return(Brightness_(arr, tau = tau, mst = mst, filt = filt))
  }
  brightness.args <- list(arr3d = ListChannels(arr, n.ch), tau = tau, mst = mst,
                          filt = filt)
  for (i in seq_along(brightness.args)) {
    if (is.null(brightness.args[[i]])) {
      brightness.args[[i]] <- list(NULL)[rep(1, length(brightness.args$arr3d))]
    }
  }
  brightness.args %>%
    purrr::pmap(Brightness_) %>%
    ChannelList2Arr
}

Brightness_ <- function(arr3d, tau = NA, mst = NULL, filt = NULL) {
  d <- dim(arr3d)
  if (length(d) != 3)
    stop("arr3d must be a three-dimensional array")
  if (!is.null(mst)) {
    arr3d <- autothresholdr::mean_stack_thresh(arr3d, method = mst, fail = NA)
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
  brightness <- VarPillars(arr3d) / MeanPillars(arr3d)
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
#' `frames.per.set = 50`, then the last 25 are discarded. If bleaching or/and
#' thresholding are selected, they are performed on the whole image stack before
#' the sectioning is done for calculation of brightnesses.
#'
#' @inheritParams Brightness
#' @param frames.per.set The number of frames with which to calculate the
#'   successive brightnesses.
#' @param mcc The number of cores to use for the parallel processing.
#' @param seed If using parallel processing (`mcc` > 1), a seed for the random
#'   number generation for [BestTau]. Don't use [set.seed], it won't work.
#'
#' @return An array where the \eqn{i}th slice is the \eqn{i}th brightness image.
#' @seealso [Brightness()].
#'
#' @examples
#' library(magrittr)
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' EBImage::display(EBImage::normalize(img[, , 1]), method = 'raster')
#' bts <- BrightnessTimeSeries(img, 20, tau = NA, mst = "Huang", mcc = 2)
#' @export
BrightnessTimeSeries <- function(arr, frames.per.set, tau = NA,
                                 mst = NULL, filt = NULL,
                                 n.ch = 1, verbose = FALSE,
                                 mcc = parallel::detectCores(), seed = NULL) {
  if (is.character(arr)) {
    if (verbose)
      message(paste0("Now processing: ", arr, "."))
    arr <- ReadImageData(arr)
  }
  d <- dim(arr)
  if (n.ch == 1) {
    if (length(d) != 3) {
      stop("Expected a three-dimensional image array image but got a ",
           length(d), " dimensional one.")
    }
  }
  if (n.ch > 1) {
    ld <- length(d)
    if (! ld %in% c(3, 4)) {
      stop("There is a problem with your image. It was read in as a ", ld,
           "-dimensional image. If read in as 4-dimensional, it is left alone, ",
           "or if read in as 3-dimensional, it is coerced to 4-dimansional ",
           "via nandb::ForceChannels", "but for ", ld, "-dimensional image, ",
           "there's nothing to be done.")
    }
    if (ld == 3) arr <- ForceChannels(arr, n.ch)
    d <- dim(arr)
  }
  if (length(d) == 3) {
    return(BrightnessTimeSeries_(arr, frames.per.set = frames.per.set,
                                 tau = tau, mst = mst, filt = filt,
                                 verbose = verbose, mcc = mcc, seed = seed))
  }
  bts.args <- list(arr3d = ListChannels(arr, n.ch), frames.per.set = frames.per.set,
                   tau = tau, mst = mst, filt = filt,
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
                                  mst = NULL, filt = NULL, verbose = FALSE,
                                  mcc = parallel::detectCores(), seed = NULL) {
  d <- dim(arr3d)
  if (length(d) != 3) stop("arr3d must be a three-dimensional array")
  if (frames.per.set > d[3]) {
    stop("frames.per.set must not be greater than the depth of arr3d")
  }
  if (!is.null(mst)) arr3d <- autothresholdr::mean_stack_thresh(arr3d, mst)
  tau.auto <- FALSE
  if (is.character(tau)) if (startsWith("auto", tolower(tau))) tau.auto <- TRUE
  if (tau.auto) tau <- BestTau(arr3d)
  if (!is.na(tau)) arr3d <- CorrectForBleaching(arr3d, tau)
  n.sets <- d[3] %/% frames.per.set
  set.indices <- lapply(seq_len(n.sets), function(x) {
    ((x - 1) * frames.per.set + 1):(x * frames.per.set)
  })
  sets <- lapply(set.indices, Slices, arr3d)
  brightnesses <- Brightnesses(sets, tau = NA, mst = NULL,
                               filt = filt, mcc = mcc, seed = seed)
  atts <- attributes(brightnesses[[1]])
  brightnesses <- Reduce(function(x, y) abind::abind(x, y, along = 3),
                         brightnesses)
  if (length(dim(brightnesses)) == 2) {
    brightnesses <- abind::abind(brightnesses, along = 3)
  }
  atts <- atts[! names(atts) %in% names(attributes(brightnesses))]
  attributes(brightnesses) <- c(attributes(brightnesses), atts)
  if (tau.auto) tau <- paste0("tau=", tau)
  attr(brightnesses, "tau") <- tau
  brightnesses
}

#' @rdname BrightnessTimeSeries
#'
#' @param folder.path The path to the folder to process (defaults to current location).
#' @param ext The file extension of the images in the folder that you wish to
#'   process. You must wish to process all files with this extension; if there
#'   are files that you don't want to process, take them out of the folder. The
#'   default is for tiff files. Do not use regular expression in this argument.
#'
#' @examples
#' setwd(tempdir())
#' WriteIntImage(img, '50.tif')
#' WriteIntImage(img, '50again.tif')
#' BrightnessTimeSeriesTxtFolder(tau = NA, mcc = 2, frames.per.set = 20)
#' list.files()
#' file.remove(list.files())  # cleanup
#' @export
BrightnessTimeSeriesTxtFolder <- function(folder.path = ".", frames.per.set,
  ext = "tif", tau = NA, mst = NULL, filt = NULL, n.ch = 1, verbose = FALSE,
  mcc = parallel::detectCores(), seed = NULL) {
  init.dir <- getwd()
  on.exit(setwd(init.dir))
  setwd(folder.path)
  if (filesstrings::str_elem(ext, 1) != ".") ext <- paste0(".", ext)
  ext <- ore::ore_escape(ext) %>% paste0("$")
  file.names <- list.files(pattern = ext)
  btss <- BrightnessTimeSeriess(file.names, tau = tau, mst = mst,
                                filt = filt, seed = seed, n.ch = n.ch,
                                frames.per.set = frames.per.set)
  frames <- vapply(btss, function(x) attr(x, "frames"), integer(1))
  tau <- simplify2array(lapply(btss, function(x) attr(x, "tau")))
  names.noext.btss <- vapply(file.names, filesstrings::before_last_dot,
                                   character(1)) %>%
    paste0("_brightness_frames=", frames, "_tau=", tau, "_mst=",
           mst, "_filter=", ifelse(is.null(filt), NA, filt))
  mapply(WriteImageTxt, btss, names.noext.btss) %>%
    invisible
  if (verbose) message("Done. Please check folder.")
}
BrightnessTimeSeriess <- function(arr3d.list, frames.per.set, tau = NA,
                                  mst = NULL, filt = NULL,
                                  n.ch = 1, verbose = FALSE,
                                  mcc = parallel::detectCores(), seed = NULL) {
  if (is.null(mst)) {
    btss <- BiocParallel::bplapply(arr3d.list, BrightnessTimeSeries, tau = tau,
                                   filt = filt, n.ch = n.ch, verbose = verbose,
                                   frames.per.set = frames.per.set, BPPARAM = bpp(mcc, seed = seed))
  } else if (is.list(arr3d.list)) {
    arr3d.list <- lapply(arr3d.list, autothresholdr::mean_stack_thresh,
                         method = mst, fail = NA)
    btss <- BiocParallel::bplapply(arr3d.list, BrightnessTimeSeries, tau = tau,
                                   filt = filt, n.ch = n.ch, verbose = verbose,
                                   frames.per.set = frames.per.set,
                                   BPPARAM = bpp(mcc, seed = seed))
  } else {
    if (!is.character(arr3d.list)) {
      stop("arr3d.list must either be a list of 3d arrays, ",
           "or a character vector of paths to the locations ",
           "of 3d arrays on disk.")
    }
    btss <- list()
    sets <- seq_along(arr3d.list) %>% {
      split(., ((. - 1) %/% mcc) + 1)
    }
    for (i in sets) {
      arrays <- lapply(arr3d.list[i], ReadImageData)
      threshed <- lapply(arrays, autothresholdr::mean_stack_thresh, method = mst,
                         fail = NA)
      bts.i <- BiocParallel::bplapply(threshed, BrightnessTimeSeries, tau = tau,
        filt = filt, n.ch = n.ch, verbose = verbose,
        frames.per.set = frames.per.set, BPPARAM = bpp(mcc, seed = seed))
      btss[i] <- bts.i
    }
  }
  btss <- lapply(btss, function(x) {
    attr(x, "mst") <- mst
    x
  })
  btss
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
#' @param seed If using parallel processing (`mcc` > 1), a seed for the random
#'   number generation for [BestTau]. Don't use [set.seed], it won't work.
#'
#' @examples
#' setwd(tempdir())
#' WriteIntImage(img, '50.tif')
#' WriteIntImage(img, '50again.tif')
#' BrightnessTxtFolder(tau = NA, mcc = 2)
#' list.files()
#' file.remove(list.files())  # cleanup
#' @export
BrightnessTxtFolder <- function(folder.path = ".", tau = NA,
  mst = NULL, filt = NULL, ext = "tif", n.ch = 1,
  mcc = parallel::detectCores(), verbose = FALSE, seed = NULL) {
  init.dir <- getwd()
  on.exit(setwd(init.dir))
  setwd(folder.path)
  if (filesstrings::str_elem(ext, 1) != ".") ext <- paste0(".", ext)
  ext <- ore::ore_escape(ext) %>% paste0("$")
  file.names <- list.files(pattern = ext)
  brightnesses <- Brightnesses(file.names, tau = tau, mst = mst,
    filt = filt, seed = seed, n.ch = n.ch)
  frames <- vapply(brightnesses, function(x) attr(x, "frames"), integer(1))
  tau <- simplify2array(lapply(brightnesses, function(x) attr(x, "tau")))
  mst <- purrr::map(brightnesses, ~ attr(., "mst")) %>% unlist
  names.noext.brightness <- vapply(file.names, filesstrings::before_last_dot,
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
#' brightnesses <- Brightnesses(img.paths, tau = NA, mcc = 2)
#'
#' @export
Brightnesses <- function(arr3d.list, tau = NA, mst = NULL,
                         filt = NULL, n.ch = 1, verbose = FALSE,
                         mcc = parallel::detectCores(), seed = NULL) {
  if (is.null(mst)) {
    brightnesses <- BiocParallel::bplapply(arr3d.list, Brightness, tau = tau,
                                           filt = filt, verbose = verbose,
                                           n.ch = n.ch,
                                           BPPARAM = bpp(mcc, seed = seed))
  } else if (is.list(arr3d.list)) {
    arr3d.list <- lapply(arr3d.list, autothresholdr::mean_stack_thresh,
                         method = mst, fail = NA)
    brightnesses <- BiocParallel::bplapply(arr3d.list, Brightness, tau = tau,
                                           filt = filt, verbose = verbose,
                                           n.ch = n.ch,
                                           BPPARAM = bpp(mcc, seed = seed))
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
      threshed <- lapply(arrays, autothresholdr::mean_stack_thresh, method = mst,
                         fail = NA)
      brightnesses.i <- BiocParallel::bplapply(threshed, Brightness, tau = tau,
                                               filt = filt, verbose = verbose,
                                               n.ch = n.ch,
                                               BPPARAM = bpp(mcc, seed = seed))
      brightnesses[i] <- brightnesses.i
    }
  }
  if (is.null(mst)) mst <- NA
  brightnesses <- lapply(brightnesses, function(x) {
    attr(x, "mst") <- mst
    x
  })
  brightnesses
}

