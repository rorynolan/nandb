#' Detrend an image series
#'
#' `CorrectForBleaching` applies detrending to an image time series using the
#' method described in Nolan et al. 2017. `CorrectForBleachingFolder` performs
#' this correction on all images in a folder, writing the corrected images to
#' disk.
#'
#' If you wish to apply thresholding and bleaching correction, aplpluy the
#' thresholding first. `CorrectforBleachingFolder` takes care of this for you.
#'
#' @param mat3d A 3-dimensional array (the image stack) where the \eqn{n}th
#'   slice is the \eqn{n}th image in the time series. To perform this on a file
#'   that has not yet been read in, set this argument to the path to that file
#'   (a string).
#' @param tau The time constant for the exponential filtering. If this is set to
#'   `'auto'`, then the value of `tau` is calculated automatically via
#'   [BestTau()].
#'
#' @references Stroud, P. D.: A recursive exponential filter for time-sensitive
#'   data, Los Alamos national Laboratory, LAUR-99-5573,
#'   \url{public.lanl.gov/stroud/ExpFilter/ExpFilter995573.pdf}, 1999.
#'
#'   img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#'   tau10 <- CorrectForBleaching(img, 10) autotau <- CorrectForBleaching(img,
#'   'auto')
#'
#' @return `CorrectForBleaching` returns the detrended image series.
#'
#' @examples
#' img <- ReadImageData(system.file("extdata", "50.tif", package = "nandb"))
#' tau10 <- CorrectForBleaching(img, 10)
#' autotau <- CorrectForBleaching(img, "auto")
#'
#' @export
CorrectForBleaching <- function(mat3d, tau) {
  d <- dim(mat3d)
  if (length(d) != 3) stop("mat3d must be a three-dimensional array")
  auto.tau <- FALSE
  if (is.na(tau)) {
    attr(mat3d, "tau") <- NA
    return(mat3d)  # option to do nothing
  }
  if (is.character(tau)) {
    tau <- tolower(tau)
    if (startsWith("auto", tau)) {
      tau <- BestTau(mat3d)
      auto.tau <- TRUE
    } else {
      stop("If tau is a string, it must be 'auto'.")
    }
  } else if (!is.numeric(tau)) {
    stop("If tau is not numeric, then it must be NA or 'auto'.")
  }
  if (is.na(tau)) {
    corrected <- mat3d
  } else {
    smoothed <- ExpSmoothPillars(mat3d, tau)
    filtered <- mat3d - smoothed
    means <- MeanPillars(mat3d)
    corrected <- filtered + as.vector(means)
    # as it so happens, this will add the means to the pillars as we desire
    corrected[corrected < 0] <- 0
    corrected <- round(corrected)  # return the array back to integer counts
  }
  attr(corrected, "tau") <- ifelse(auto.tau, paste0("auto=", round(tau)),
                                   as.character(tau))
  corrected
}

#' @rdname CorrectForBleaching
#' @param folder.path The path (relative or absolute) to the folder you wish to
#'   process.
#' @param mst Do you want to apply an intensity threshold prior to correcting
#'   for bleaching (via [MeanStackThresh()])? If so, set your thresholding
#'   \emph{method} here.
#' @param ext the file extension of the images in the folder that you wish to
#'   process (can be rooted in regular expression for extra-safety, as in the
#'   default). You must wish to process all files with this extension; if there
#'   are files that you don't want to process, take them out of the folder. The
#'   default is for tiff files.
#' @param na How do you want to treat `NA` values? R can only write integer
#'   values (and hence not `NA`s) to tiff pixels. `na = 'saturate'` sets them to
#'   saturated value. `na = 'zero'` sets them to zero, while `na = 'error'` will
#'   give an error if the image contains `NA`s. Note that if you threshold, you
#'   are almost certain to get `NA`s.
#' @param mcc The number of cores to use for the parallel processing.
#' @param seed A seed for the random number generation for [BestTau]. Don't use
#'   [set.seed], it won't work.
#'
#'
#' @examples
#' setwd(tempdir())
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' WriteIntImage(img, '50.tif')
#' WriteIntImage(img, '50again.tif')
#' set.seed(33)
#' CorrectForBleachingFolder(tau = 'auto', mst = 'tri', mcc = 2, na = "s")
#' list.files()
#' file.remove(list.files())
#'
#' @export
CorrectForBleachingFolder <- function(folder.path = ".", tau = NA, mst = NULL,
                                      ext = "\\.tif$", na = "error",
                                      mcc = parallel::detectCores(),
                                      seed = NULL) {
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(folder.path)
  file.paths <- list.files(pattern = ext)
  BiocParallel::bplapply(file.paths, CorrectForBleachingFile,
                         tau = tau, mst = mst, na = na,
                         BPPARAM = bpp(mcc, seed = seed))
}

CorrectForBleachingFile <- function(file.path, tau = NA, mst = NULL,
                                    na = "error") {
  arr3d <- ReadImageData(file.path)
  stopifnot(length(dim(arr3d)) == 3)
  if (!is.null(mst)) arr3d <- MeanStackThresh(arr3d, mst)
  corrected <- CorrectForBleaching(arr3d, tau)
  out.name <- filesstrings::BeforeLastDot(file.path) %>%
    paste0("_tau=", attr(corrected, "tau"), ".",
           filesstrings::StrAfterNth(file.path, stringr::coll("."), -1))
  WriteIntImage(corrected, out.name, na = na)
}

#' Find the best tau for exponential filtering detrend.
#'
#' Say you have an image series that you wish to detrend before performing a
#' brightness calculation. This function finds the best `tau` for an
#' exponential filtering detrend. See \code{vignette('Adaptive Detrending',
#' package = 'nandb')} for more details.
#'
#' @param mat3d A 3-dimensional array (the image stack) where the \eqn{n}th
#'   slice is the \eqn{n}th image in the time series. To perform this on a file
#'   that has not yet been read in, set this argument to the path to that file
#'   (a string).
#' @param mst Do you want to apply an intensity threshold prior to calculating
#'   brightness (via [MeanStackThresh()])? If so, set your
#'   thresholding \emph{method} here.
#' @param tol What size of error in the estimate of the \emph{ideal} `tau`
#'   (aside from the error introduced by the random image simulation, see
#'   `vignette('AdaptiveDetrend', package = 'nandb')`) are you willing
#'   to tolerate? The default is 1.
#'
#' @return A number. The estimate of the ideal `tau` to use, with an
#'   attribute '`brightness.immobile`' giving the brightness of the
#'   simulated (from all immobile particles) image series after detrending with
#'   this `tau` (this should be very close to 1).
#'
#' @examples
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' BestTau(img, mst = 'tri', tol = 3)
#'
#' @export
BestTau <- function(mat3d, mst = NULL, tol = 1) {
  if (is.character(mat3d)) {
    mat3d <- ReadImageData(mat3d[1])
  }
  if (!is.null(mst))
    mat3d <- MeanStackThresh(mat3d, mst, fail = NA)
  raw.brightness.mean <- Brightness(mat3d) %>% mean(na.rm = TRUE)
  if (raw.brightness.mean < 1) {
    stop("Your raw brightness mean is below 1,",
         "there's probably something wrong with your acquisition.")
  }
  means <- apply(mat3d, 3, mean, na.rm = TRUE)
  p <- sum(!is.na(mat3d[, , 1]))
  sim.img.arr <- vapply(means, stats::rpois, numeric(p), n = p)
  BrightnessMeanSimMatTau <- function(sim.img.arr, tau) {
    if (is.na(tau)) {
      detrended <- sim.img.arr
    } else {
      detrended <- sim.img.arr - ExpSmoothRows(sim.img.arr,
        tau) + rowMeans(sim.img.arr)
    }
    brightnesses <- matrixStats::rowVars(detrended) / rowMeans(detrended)
    mean(brightnesses)
  }
  notau.sim.brightness.mean <- BrightnessMeanSimMatTau(sim.img.arr, NA)
  if (notau.sim.brightness.mean <= 1) {
    # This is when the original image had mean brightness greater than 1 but the
    # simulated image had mean brightness <= 1.
    # In this case, no bleaching is detected, therefore no detrend is needed
    return(NA)
  }
  tau2.sim.brightness.mean <- BrightnessMeanSimMatTau(sim.img.arr,
    2)
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
  while (BrightnessMeanSimMatTau(sim.img.arr, upper) < 1) {
    lower <- upper
    upper <- 2 * upper
  }
  TauFarFromOne <- function(tau, sim.img.arr) {
    BrightnessMeanSimMatTau(sim.img.arr, tau) - 1
  }
  root <- stats::uniroot(TauFarFromOne, c(lower, upper), sim.img.arr,
    tol = tol, extendInt = "upX")
  tau <- root$root
  attr(tau, "brightness.immobile") <- 1 + root$f.root
  tau
}
