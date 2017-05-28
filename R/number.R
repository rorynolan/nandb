#' Calculate number from image series.
#'
#' Given a time stack of images, `Number` performs a calculation of the
#' number for each pixel. `NumberTxtFolder` does this
#' calculation for an entire folder, writing the results as text files via
#' [WriteImageTxt()]. `Numbers` calculates the numbers for
#' several image series in parallel.
#'
#' Do not try to parallelize the use of `Number` and friends yourself (e.g.
#' with [mclapply()]) because it will throw an error (this is because
#' the `autothresholdr` package does not work in parallel (indeed anything
#' run using the `rJava` package won't)). Always use `Numbers` for
#' this purpose (this has a workaround whereby it does the thresholding in
#' series and does the rest in parallel).
#'
#' @param arr A 3-dimensional array (the image stack) where the \eqn{n}th
#'   slice is the \eqn{n}th image in the time series. To perform this on a file
#'   that has not yet been read in, set this argument to the path to that file
#'   (a string).
#' @param n.ch The number of channels in the image (default 1).
#' @param tau If this is specified, bleaching correction is performed with
#'   [CorrectForBleaching()] which uses exponential filtering with
#'   time constant `tau` (where the unit of time is the time between
#'   frames). If this is set to `'auto'`, then the value of `tau` is
#'   calculated automatically via [BestTau()].
#' @param mst Do you want to apply an intensity threshold prior to calculating
#'   number (via [autothresholdr::mean_stack_thresh()])? If so, set your thresholding
#'   \emph{method} here.
#' @param fail If thresholding is done, to which value should pixels not
#'   exceeeding the threshold be set?
#' @param filt Do you want to smooth (`filt = 'smooth'`) or median
#'   (`filt = 'median'`) filter the number image using
#'   [SmoothFilterB()] or [MedianFilterB()] respectively? If
#'   selected, these are invoked here with a filter radius of 1 and with the
#'   option `na_count = TRUE`. If you want to smooth/median filter the
#'   number image in a different way, first calculate the numbers without
#'   filtering (`filt = NULL`) using this function and then perform your
#'   desired filtering routine on the result.
#' @param verbose If arr3d is specified as a file name, print a message to tell
#'   the user that that file is now being processed (useful for
#'   `NumberTxtFolder`, does not work with multiple cores).
#'
#' @return `Number` returns a matrix, the number image; `Numbers`
#'   returns a list of these. The result of `NumberTxtFolder` is the text
#'   csv files written to disk (in the same folder as the input images).
#'
#' @examples
#' library(EBImage)
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' display(normalize(img[, , 1]), method = 'raster')
#' number <- Number(img, tau = NA, mst = "Huang")
#' @export
Number <- function(arr, tau = NA, mst = NULL,
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
    return(Number_(arr, tau = tau, mst = mst, filt = filt))
  }
  number.args <- list(arr3d = ListChannels(arr), tau = tau, mst = mst,
                      filt = filt)
  for (i in seq_along(number.args)) {
    if (is.null(number.args[[i]])) {
      number.args[[i]] <- list(NULL)[rep(1, length(number.args$arr3d))]
    }
  }
  number.args %>%
    purrr::pmap(Number_) %>%
    ChannelList2Arr
}

Number_ <- function(arr3d, tau = NA, mst = NULL,
  fail = NA, filt = NULL) {
  d <- dim(arr3d)
  if (length(d) != 3)
    stop("arr3d must be a three-dimensional array")
  if (!is.null(mst)) {
    arr3d <- autothresholdr::mean_stack_thresh(arr3d, method = mst, fail = fail)
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
  number <- MeanPillars(arr3d) ^ 2 / VarPillars(arr3d)
  if (!is.null(filt)) {
    allowed <- c("median", "smooth")
    filt <- tolower(filt)
    sw <- startsWith(allowed, filt)
    if (!any(sw))
      stop("filt must be either 'median' or 'smooth'")
    filt <- allowed[sw]
    if (filt == "median") {
      number <- MedianFilterB(number, na_count = TRUE)
    } else {
      number <- SmoothFilterB(number, na_count = TRUE)
    }
  }
  attributes(number) <- c(attributes(number), list(frames = d[3],
    tau = as.character(ifelse(tau.auto, ifelse(is.na(tau), "auto=NA",
      stringr::str_c("auto=", round(tau))), tau)),
    filter = ifelse(is.null(filt), NA, filt),
    mst = ifelse(is.null(mst), NA, mst)))
  number
}

#' Create a number time-series.
#'
#' Given a stack of images `arr3d`, use the first `frames.per.set` of them to
#' create one number image, the next `frames.per.set` of them to create
#' the next number image and so on to get a time-series of number
#' images. If `tau` is specified, bleaching correction is performed via
#' [CorrectForBleaching()].
#'
#' This may discard some images, for example if 175 frames are in the input and
#' `frames.per.set = 50`, then the last 25 are discarded. If bleaching
#' correction is selected, it is performed on the whole image stack before the
#' sectioning is done for calculation of numbers.
#'
#' @inheritParams Number
#' @param frames.per.set The number of frames with which to calculate the
#'   successive numbers.
#' @param n.ch The number of channels in the image (default 1).
#' @param mcc The number of cores to use for the parallel processing.
#' @param seed A seed for the random number generation for [BestTau]. Don't use
#'   [set.seed], it won't work.
#'
#' @return An array where the \eqn{i}th slice is the \eqn{i}th number image.
#' @seealso [Number()].
#'
#' @examples
#' library(EBImage)
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' bts <- NumberTimeSeries(img, 20, tau = NA, mst = "Huang", mcc = 2)
#' @export
NumberTimeSeries <- function(arr, frames.per.set, tau = NA,
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
    return(NumberTimeSeries_(arr, frames.per.set = frames.per.set,
                                 tau = tau, mst = mst,
                                 filt = filt, mcc = mcc, seed = seed))
  }
  nts.args <- list(arr3d = ListChannels(arr), frames.per.set = frames.per.set,
                   tau = tau, mst = mst, filt = filt, mcc = mcc, seed = seed)
  for (i in seq_along(nts.args)) {
    if (is.null(nts.args[[i]])) {
      nts.args[[i]] <- list(NULL)[rep(1, length(nts.args$arr3d))]
    }
  }
  nts.args %>%
    purrr::pmap(NumberTimeSeries_) %>%
    ChannelList2Arr
}

NumberTimeSeries_ <- function(arr3d, frames.per.set, tau = NA,
                              mst = NULL, filt = NULL,
                              mcc = parallel::detectCores(), seed = NULL) {
  d <- dim(arr3d)
  if (length(d) != 3)
    stop("arr3d must be a three-dimensional array")
  if (frames.per.set > d[3]) {
    stop("frames.per.set must not be greater than the depth of arr3d")
  }
  if (!is.na(tau))
    arr3d <- CorrectForBleaching(arr3d, tau)
  n.sets <- floor(d[3]/frames.per.set)
  set.indices <- lapply(seq_len(n.sets), function(x) {
    ((x - 1) * frames.per.set + 1):(x * frames.per.set)
  })
  sets <- lapply(set.indices, Slices, arr3d)
  numbers <- Numbers(sets, tau = NA, mst = NULL, filt = filt,
                     mcc = mcc, seed = seed) %>%
    Reduce(function(x, y) abind::abind(x, y, along = 3), .)
  if (length(dim(numbers)) == 2) {
    numbers <- abind::abind(numbers, along = 3)
  }
  numbers
}

#' @rdname Number
#'
#' @param folder.path The path (relative or absolute) to the folder you wish to
#'   process.
#' @param ext The file extension of the images in the folder that you wish to
#'   process. You must wish to process all files with this extension; if there
#'   are files that you don't want to process, take them out of the folder. The
#'   default is for tiff files. Do not use regular expression in this argument.
#'
#' @examples
#' setwd(tempdir())
#' WriteIntImage(img, '50.tif')
#' WriteIntImage(img, '50again.tif')
#' NumberTxtFolder(tau = NA, mst = "Huang", mcc = 2)
#' file.remove(list.files())  # cleanup
#' @export
NumberTxtFolder <- function(folder.path = ".", tau = NA, mst = NULL,
                            filt = NULL, ext = "tif",
                            mcc = parallel::detectCores(), seed = NULL,
                            verbose = FALSE) {
  init.dir <- getwd()
  on.exit(setwd(init.dir))
  setwd(folder.path)
  if (filesstrings::StrElem(ext, 1) != ".") ext <- paste0(".", ext)
  ext <- ore::ore_escape(ext) %>% paste0("$")
  file.names <- list.files(pattern = ext)
  numbers <- Numbers(file.names, tau = tau, mst = mst,
                     filt = filt, mcc = mcc, seed = seed)
  frames <- vapply(numbers, function(x) attr(x, "frames"), integer(1))
  tau <- vapply(numbers, function(x) attr(x, "tau"), character(1))
  mst <- purrr::map(numbers, ~ attr(., "mst")) %>% unlist
  names.noext.number <- vapply(file.names, filesstrings::BeforeLastDot,
                               character(1)) %>%
    paste0("_number_frames=", frames, "_tau=", tau, "_mst=",
      mst, "_filter=", ifelse(is.null(filt), NA, filt))
  mapply(WriteImageTxt, numbers, names.noext.number) %>% invisible
  if (verbose)
    message("Done. Please check folder.")
}

#' @rdname Number
#'
#' @param arr3d.list A list of 3-dimensional arrays. To perform this on files
#'   that have not yet been read in, set this argument to the path to these
#'   files (a character vector).
#' @param mcc The number of cores to use for the parallel processing.
#' @param seed A seed for the random number generation for [BestTau]. Don't use
#'   [set.seed], it won't work.
#'
#' @examples
#' img.paths <- rep(system.file('extdata', '50.tif', package = 'nandb'), 2)
#' numbers <- Numbers(img.paths, mst = 'otsu', mcc = 2)
#'
#' @export
Numbers <- function(arr3d.list, tau = NA, mst = NULL, fail = NA, filt = NULL,                          verbose = FALSE, mcc = parallel::detectCores(),
                    seed = NULL) {
  if (is.null(mst)) {
    numbers <- BiocParallel::bplapply(arr3d.list, Number, tau = tau,
                filt = filt, verbose = verbose, BPPARAM = bpp(mcc, seed = seed))
  } else if (is.list(arr3d.list)) {
    arr3d.list <- lapply(arr3d.list, autothresholdr::mean_stack_thresh,
                         method = mst, fail = fail)
    numbers <- BiocParallel::bplapply(arr3d.list, Number, tau = tau,
                filt = filt, verbose = verbose, BPPARAM = bpp(mcc, seed = seed))
  } else {
    if (!is.character(arr3d.list)) {
      stop("arr3d.list must either be a list of 3d arrays, ",
        "or a character vector of paths to the locations ",
        "of 3d arrays on disk.")
    }
    numbers <- list()
    sets <- seq_along(arr3d.list) %>% {
      split(., ((. - 1) %/% mcc) + 1)
    }
    for (i in sets) {
      arrays <- lapply(arr3d.list[i], ReadImageData)
      threshed <- lapply(arrays, autothresholdr::mean_stack_thresh, method = mst,
                         fail = fail)
      numbers.i <- BiocParallel::bplapply(threshed, Number, tau = tau,
        filt = filt, verbose = verbose, BPPARAM = bpp(mcc, seed = seed))
      numbers[i] <- numbers.i
    }
  }
  if (is.null(mst)) mst <- NA
  numbers <- lapply(numbers, function(x) {
    attr(x, "mst") <- mst
    x
  })
  numbers
}
