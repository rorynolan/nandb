#' Detrend an image series
#'
#' `CorrectForBleaching()` applies detrending to an image time series using the
#' method described in Nolan et al. 2017. `CorrectForBleachingFolder()` performs
#' this correction on all images in a folder, writing the corrected images to
#' disk.
#'
#' If you wish to apply thresholding and bleaching correction, aplpluy the
#' thresholding first. `CorrectforBleachingFolder` takes care of this for you.
#'
#' @param arr An array, can be 3- or 4-dimensional. The first two slots give the
#'   x- and y-cordinates of pixels respectively. If the array is 3-dimensional,
#'   the third slot gives the index of the frame. If it is 4-dimensional, the
#'   third slot indexes the channel and the fourth indexes the frame in the
#'   stack.. To perform this on a file that has not yet been read in, set this
#'   argument to the path to that file (a string).
#' @param tau The time constant for the exponential filtering. If this is set to
#'   `'auto'`, then the value of `tau` is calculated automatically via
#'   [BestTau()].
#' @param n.ch The number of channels in the image (default 1).
#' @param mst Do you want to apply an intensity threshold prior to calculating
#'   brightness (via [autothresholdr::mean_stack_thresh()])? If so, set your
#'   thresholding \emph{method} here.
#' @param fail When using `auto_thresh_apply_mask`, to what value do you wish to
#'   set the pixels which fail to exceed the threshold. `fail = 'saturate'` sets
#'   them to saturated value (see "Details"). `fail = 'zero'` sets them to zero.
#'   You can also specify directly here a natural number (must be between 0 and
#'   2 ^ 16 - 1) to use in place of `NA`s.
#' @param ignore_black Ignore black pixels/elements (zeros) when performing the
#'   thresholding?
#' @param ignore_white Ignore white pixels when performing the thresholding? If
#'   set to `TRUE`, the function makes a good guess as to what the white
#'   (saturated) value would be (see "Details"). If this is set to a number, all
#'   pixels with value greater than or equal to that number are ignored.
#' @param mcc The number of cores to use for the parallel processing (this
#'   function can be parallelized over channels.
#' @param seed A seed for the random number generation for [BestTau]. Don't use
#'   [set.seed], it won't work.
#'
#' @references Stroud, P. D.: A recursive exponential filter for time-sensitive
#'   data, Los Alamos national Laboratory, LAUR-99-5573,
#'   \url{public.lanl.gov/stroud/ExpFilter/ExpFilter995573.pdf}, 1999.
#'
#' @return `CorrectForBleaching()` returns the detrended image series.
#'   `CorrectForBleachingFolder()` returns a character vector of the input tifs
#'   invisibly.
#'
#' @examples
#' library(magrittr)
#' img <- ReadImageData(system.file("extdata", "50.tif", package = "nandb"))
#' autotau <- CorrectForBleaching(img, "auto", mst = "huang")
#' img <- ReadImageData(system.file("extdata", "two_ch.tif", package = "nandb"))
#' tau10_tri <- CorrectForBleaching(img, 10, n.ch = 2, mst = "triangle")
#' @export
CorrectForBleaching <- function(arr, tau, n.ch = 1,
                                mst = NULL, fail = NA,
                                ignore_black = FALSE,
                                ignore_white = FALSE,
                                mcc = 1, seed = NULL) {
  if (is.character(arr)) arr <- ReadImageData(arr)
  checkmate::assert_array(arr, min.d = 3, max.d = 4)
  checkmate::assert_scalar(tau, na.ok = TRUE)
  if (is.character(tau)) {
    if (!startsWith("auto", tolower(tau))) {
      stop("If tau is a string, it must be 'auto'.")
    }
  } else {
    checkmate::assert_number(tau, na.ok = TRUE, lower = 0)
  }
  checkmate::assert_int(n.ch, lower = 1)
  if (n.ch > 1) {
    if (length(dim(arr)) == 3) arr <- ForceChannels(arr, n.ch)
  }
  if (!is.null(mst)) {
    if (anyNA(arr)) {
      stop("Cannot perform thresholding if arr has NA values to begin with.")
    }
  }
  fail <- TranslateFail(arr, fail)
  if (length(dim(arr)) == 3) {
    if (!is.null(mst)) {
      arr <- autothresholdr::mean_stack_thresh(arr, method = mst, fail = fail,
                                               ignore_black = ignore_black,
                                               ignore_white = ignore_white)
    }
    CorrectForBleaching_(arr, tau)
  } else {
    channel.list <- ListChannels(arr, n.ch = n.ch)
    if (length(tau) != 1) {
      if (length(tau) != length(channel.list)) {
        stop("Tau must either be length 1 or ",
             "have the same length as the number of channels in arr.")
      }
    }
    if (!is.null(mst)) {
      channel.list <- mapply(autothresholdr::mean_stack_thresh, channel.list,
                      method = mst, fail = NA,
                      ignore_black = ignore_black,
                      ignore_white = ignore_white,
                      SIMPLIFY = FALSE) %>%
        BiocParallel::bpmapply(CorrectForBleaching_, ., tau, SIMPLIFY = FALSE,
                               BPPARAM = bpp(mcc, seed))
      for (i in seq_along(channel.list)) {
        channel.list[[i]][is.na(channel.list[[i]])] <- fail
      }
    } else {
      channel.list <- BiocParallel::bpmapply(CorrectForBleaching_,
                                             channel.list, tau,
                                             BPPARAM = bpp(mcc, seed),
                                             SIMPLIFY = FALSE)
    }
    ChannelList2Arr(channel.list) %>%
      radiant.data::set_attr("tau",
                             unlist(purrr::map(channel.list, attr, "tau")))
  }
}
CorrectForBleaching_ <- function(arr3d, tau) {
  d <- dim(arr3d)
  if (length(d) != 3) stop("arr3d must be a three-dimensional array")
  auto.tau <- FALSE
  if (is.character(tau)) {
    tau <- tolower(tau)
    if (startsWith("auto", tau)) {
      tau <- BestTau(arr3d)
      auto.tau <- TRUE
    } else {
      stop("If tau is a string, it must be 'auto'.")
    }
  } else if ((!is.numeric(tau)) && (!is.na(tau))) {
    stop("If tau is not numeric, then it must be NA or 'auto'.")
  }
  if (is.na(tau)) {
    corrected <- arr3d
  } else {
    smoothed <- ExpSmoothPillars(arr3d, tau)
    filtered <- arr3d - smoothed
    means <- MeanPillars(arr3d)
    corrected <- filtered + as.vector(means)
    # as it so happens, this will add the means to the pillars as we desire
    corrected[corrected < 0] <- 0
    corrected <- round(corrected)  # return the array back to integer counts
  }
  attr(corrected, "tau") <- ifelse(auto.tau, paste0("auto=", round(tau)),
                                   as.character(tau))
  corrected
}


CorrectForBleachings <- function(arr.list, tau, n.ch = 1,
                                 mst = NULL, fail = NA,
                                 ignore_black = FALSE,
                                 ignore_white = FALSE,
                                 mcc = 1, seed = NULL) {
  if (is.character(arr.list)) arr.list <- purrr:::map(arr.list, ReadImageData)
  arr.list <- purrr::map(arr.list, ForceChannels, n.ch)
  if (!is.null(mst)) {
    if (! length(mst) %in% c(1, n.ch)) {
      stop("If mst is specified, it must be length 1 or have length ",
           "equal to the number of channels n.ch.")
    }
    if (length(mst) == 1) mst <- as.list(mst)[rep(1, n.ch)]
    for (i in seq_along(arr.list)) {
      if (n.ch == 1) {
        arr.list[[i]] <- autothresholdr::mean_stack_thresh(arr.list[[i]],
                                                           mst[[1]])
      } else {
        for (j in seq_len(dim(arr.list[[i]])[3])) {
          arr.list[[i]][, , j, ] <- autothresholdr::mean_stack_thresh(
                                      arr.list[[i]][, , j, ], mst[[j]])
        }
      }
    }
    for (i in seq_along(mst)) {
      if (is.null(mst[[i]])) mst[[i]] <- NA
    }
  }
  BiocParallel::bplapply(arr.list, CorrectForBleaching, tau = tau,
                         n.ch = n.ch, mst = NULL, fail = fail,
                         ignore_black = ignore_black,
                         ignore_white = ignore_white,
                         BPPARAM = bpp(mcc, seed)) %>%
    lapply(radiant.data::set_attr, "mst", unlist(mst))
}

#' @rdname CorrectForBleaching
#' @param folder.path The path (relative or absolute) to the folder you wish to
#'   process.
#' @param ext The file extension of the images in the folder that you wish to
#'   process. You must wish to process all files with this extension; if there
#'   are files that you don't want to process, take them out of the folder. The
#'   default is for tiff files. Do not use regular expression in this argument.
#' @param na How do you want to treat `NA` values? R can only write integer
#'   values (and hence not `NA`s) to tiff pixels. `na = 'saturate'` sets them to
#'   saturated value. `na = 'zero'` sets them to zero, while `na = 'error'` will
#'   give an error if the image contains `NA`s. Note that if you threshold, you
#'   are almost certain to get `NA`s.
#'
#' @examples
#' setwd(tempdir())
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' WriteIntImage(img, '50.tif')
#' WriteIntImage(img, '50again.tif')
#' set.seed(33)
#' CorrectForBleachingFolder(tau = 1000, mst = "Huang", mcc = 2, na = "s")
#' list.files()
#' file.remove(list.files())
#'
#' @export
CorrectForBleachingFolder <- function(folder.path = ".", n.ch = 1,
                                      tau = NA, mst = NULL,
                                      ignore_black = FALSE,
                                      ignore_white = FALSE,
                                      ext = "tif", na = "error",
                                      mcc = 1,
                                      seed = NULL) {
  checkmate::assert_scalar(tau, na.ok = TRUE)
  if (is.character(tau)) {
    if (!startsWith("auto", tolower(tau))) {
      stop("If tau is a string, it must be 'auto'.")
    }
  } else {
    checkmate::assert_number(tau, na.ok = TRUE, lower = 0)
  }
  checkmate::assert_int(n.ch, lower = 1)
  checkmate::assert_int(mcc, lower = 1)
  if (!is.null(mst)) {
    if (! length(mst) %in% c(1, n.ch)) {
      stop("If mst is specified, it must be length 1 or have length ",
           "equal to the number of channels n.ch.")
    }
  }
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(folder.path)
  if (filesstrings::str_elem(ext, 1) != ".") ext <- paste0(".", ext)
  ext <- ore::ore_escape(ext) %>% paste0("$")
  file.paths <- list.files(pattern = ext)
  if (is.null(mst)) mst <- list(NULL)[rep(1, length(file.paths))]
  tif.sets <- split(file.paths, (1 + seq_along(file.paths)) %/% mcc)
  for (tifs in tif.sets) {
    corrected <- CorrectForBleachings(tifs, tau = tau, n.ch = n.ch, mst = mst,
                                      ignore_black = ignore_black,
                                      ignore_white = ignore_white,
                                      mcc = mcc, seed = seed)
    corrected.names <- paste0(filesstrings::before_last_dot(tifs),
      "_detrended_mst=", paste(unlist(mst), collapse = ","),
      "_tau=", paste(purrr::map_chr(corrected,
                       ~ paste(unlist(attr(., "tau")), collapse = ","))))
    na.write <- purrr::map_int(corrected, TranslateFail, na)
    writing <- mapply(WriteIntImage, corrected, corrected.names, na = na.write)
  }
  invisible(file.paths)
}

#' Find the best tau for exponential filtering detrend.
#'
#' Say you have an image series that you wish to detrend before performing a
#' brightness calculation. This function finds the best `tau` for an exponential
#' filtering detrend. See \code{vignette('Adaptive Detrending', package =
#' 'nandb')} for more details.
#'
#' @param arr3d A 3-dimensional array (the image stack) where the \eqn{n}th
#'   slice is the \eqn{n}th image in the time series. To perform this on a file
#'   that has not yet been read in, set this argument to the path to that file
#'   (a string).
#' @param mst Do you want to apply an intensity threshold prior to calculating
#'   brightness (via [autothresholdr::mean_stack_thresh()])? If so, set your
#'   thresholding \emph{method} here.
#' @param tol What size of error in the estimate of the \emph{ideal} `tau`
#'   (aside from the error introduced by the random image simulation, see
#'   `vignette('AdaptiveDetrend', package = 'nandb')`) are you willing to
#'   tolerate? The default is 1.
#' @param mcc The number of cores to use for the parallel processing (this
#'   function can be parallelized over channels.
#' @param seed A seed for the random number generation for [BestTau]. Don't use
#'   [set.seed], it won't work.
#'
#' @return A number. The estimate of the ideal `tau` to use, with an attribute
#'   '`brightness.immobile`' giving the brightness of the simulated (from all
#'   immobile particles) image series after detrending with this `tau` (this
#'   should be very close to 1).
#'
#' @examples
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' BestTau(img, tol = 3)
#'
#' @export
BestTau <- function(arr3d, mst = NULL, tol = 1, mcc = 1, seed = NULL) {
  if (is.character(arr3d)) {
    arr3d <- ReadImageData(arr3d[1])
  }
  if (!is.null(mst)) {
    arr3d <- autothresholdr::mean_stack_thresh(arr3d, mst, fail = NA)
  }
  raw.brightness.mean <- Brightness(arr3d) %>% mean(na.rm = TRUE)
  if (raw.brightness.mean < 1) {
    stop("Your raw brightness mean is below 1,",
         "there's probably something wrong with your acquisition.")
  }
  means <- apply(arr3d, 3, mean, na.rm = TRUE)
  p <- sum(!is.na(arr3d[, , 1]))
  sim.img.arr <- vapply(means, stats::rpois, numeric(p), n = p)
  sim.img.arr.extended <- MedReflectExtendRows(sim.img.arr, TRUE, TRUE)
  BrightnessMeanSimMatTau <- function(sim.img.arr, sim.img.arr.extended, tau) {
    if (is.na(tau)) {
      detrended <- sim.img.arr
    } else {
      detrended <- (sim.img.arr - ExpSmoothRows(sim.img.arr.extended,
                                                tau, extended = TRUE) +
                      rowMeans(sim.img.arr)) %>% round
    }
    brightnesses <- matrixStats::rowVars(detrended) / rowMeans(detrended)
    mean(brightnesses)
  }
  notau.sim.brightness.mean <- BrightnessMeanSimMatTau(sim.img.arr,
                                                       sim.img.arr.extended, NA)
  if (notau.sim.brightness.mean <= 1) {
    # This is when the original image had mean brightness greater than 1 but the
    # simulated image had mean brightness <= 1.
    # In this case, no bleaching is detected, therefore no detrend is needed
    return(NA)
  }
  tau2.sim.brightness.mean <- BrightnessMeanSimMatTau(sim.img.arr,
                                                      sim.img.arr.extended, 2)
  if (tau2.sim.brightness.mean > 1) {
    stop("Even with a savage detrend of tau = 2, ",
         "the brightnesses still have mean greater than 1. ",
         "There's probably a problem with your data, ",
         "or else your region of interest wasn't properly selected ",
         "(perhaps by thresholding). ",
         "You should check this out. ",
         "If you want to work with the data as is, ",
         "then you'll have to choose your own detrend.")
  }
  lower <- 2
  upper <- length(means)
  while (BrightnessMeanSimMatTau(sim.img.arr,
                                 sim.img.arr.extended, upper) < 1) {
    lower <- upper
    upper <- 2 * upper
  }
  TauFarFromOne <- function(tau, sim.img.arr, sim.img.arr.extended) {
    BrightnessMeanSimMatTau(sim.img.arr, sim.img.arr.extended, tau) - 1
  }
  root <- stats::uniroot(TauFarFromOne, c(lower, upper), sim.img.arr,
                         sim.img.arr.extended, tol = tol, extendInt = "upX")
  tau <- root$root
  attr(tau, "brightness.immobile") <- 1 + root$f.root
  tau
}
