#' Calculate numbers of kmers based on brightnesses
#'
#' Calculate the number pixels occupied of molecules in each oligomeric state
#' (monomer, dimer, trimer, ..., kmer, ...) based on the brightnesses of those
#' pixels.
#'
#' @param brightnesses Vector or array of pixel brightnesses
#' @param monomer.med The median brightness of a monomer. You must know what
#'   this is to use this function correctly. This must be greater than 1 (this
#'   is the 'apparent brightness' of a monomer). If you're reading the Digman et
#'   al. (2008) paper, this is \eqn{1 + \epsilon}.
#'
#'   This function takes the brightnesses in the range (1, 1 + 2 *
#'   \code{monomer.med} to be monomers, those in the range (1 + 2 *
#'   \code{monomer.med}, 1 + 4 * \code{monomer.med} to be dimers, (1 + 4 *
#'   \code{monomer.med}, 1 + 6 * \code{monomer.med} to be trimers and so on.
#'
#' @return A named vector (named '1mers' '2mers' '3mers' and so on) with each
#'   element detailing the number of that kmer found.
#'
#' @examples
#' brightnesses <- runif(100, 1, 3)
#' monomer.med <- 1.2
#' KmersFromBrightnesses(brightnesses, monomer.med)
#' @export
KmersFromBrightnesses <- function(brightnesses, monomer.med) {
  if (monomer.med <= 1)
    stop("monomer.med must be greater than 1")
  window.size <- 2 * (monomer.med - 1)
  mb <- max(brightnesses, na.rm = TRUE)
  if (mb <= 1) {
    kmers <- 0
    names(kmers) <- "1mers"
  } else {
    window.breaks <- c(seq(1, mb, by = window.size), mb)
    # these values delimit the windows in which we'll count the kmers
    cuts <- cut(brightnesses, window.breaks)
    kmers <- as.vector(table(cuts))
    names(kmers) <- paste0(seq_along(kmers), "mers")
  }
  kmers
}

#' Calculate numbers of kmers based on an image time-series
#'
#' Given an image (as a file path or an array), \code{KmersFromImage} does the
#' brightness calculation via \code{\link{Brightness}} and then counts the
#' numbers of each kmer via \code{\link{KmersFromBrightnesses}}.
#' \code{KmersFromImagesFolder} does this for an entire folder (directory) of
#' images and outputs a csv file of the results.
#'
#' @param mat3d A 3-dimensional array that one would might input to
#'   \code{\link{Brightness}} \emph{or} the path to an image file on disk.
#' @param monomer.med The median brightness of a monomer. You must know what
#'   this is to use this function correctly. This must be greater than 1 (this
#'   is the 'apparent brightness' of a monomer). If you're reading the Digman et
#'   al. (2008) paper, this is \eqn{1 + \epsilon}.
#' @param tau The time constant for the exponential filtering (see
#'   \code{\link{Brightness}}).
#' @param mst Do you want to apply an intensity threshold prior to calculating
#'   brightness (via \code{\link{MeanStackThresh}})? If so, set your thresholding
#'   \emph{method} here.
#' @param skip.consts An image array with only one value (a 'constant array')
#'   won't threshold properly. By default the function would give an error, but
#'   by setting this parameter to \code{TRUE}, the array would instead be
#'   skipped (the function will return the original array) and give a warning.
#' @param filt Do you want to smooth (\code{filt = 'smooth'}) or median
#'   (\code{filt = 'median'}) filter the brightness image using
#'   \code{\link{SmoothFilterB}} or \code{\link{MedianFilterB}} respectively? If
#'   selected, these are invoked here with a filter radius of 1 and with the
#'   option \code{na_count = TRUE}. If you want to smooth/median filter
#'   the brightness image in a different way, first calculate the brightnesses
#'   without filtering (\code{filt = NULL}) using this function and then perform
#'   your desired filtering routine on the result.
#' @param verbose If mat3d is specified as a file name, print a message to tell
#'   the user that that file is now being processed (useful for
#'   \code{BrightnessFolder}, does not work with multiple cores) and to tell
#'   when \code{KmersFromImagesFolder} is done.
#'
#' @return A named vector (named '1mers' '2mers' '3mers' and so on) with each
#'   element detailing the number of that kmer found, or for
#'   \code{KmersFromImagesFolder}, a csv file is written to disk detailing one
#'   of these vectors for each image. This vector also has an attribute
#'   'mean.intensity' giving the mean intensity of the input image.
#' @export
KmersFromImage <- function(mat3d, monomer.med, tau = NA, mst = NULL,
  skip.consts = FALSE, filt = NULL, verbose = TRUE) {
  if (is.character(mat3d)) {
    if (verbose)
      print(paste0("Now processing: ", mat3d, "."))
    mat3d <- ReadImageData(mat3d)
  }
  bright <- Brightness(mat3d, tau = tau, mst = mst, skip.consts = skip.consts,
    fail = NA, filt = filt)
  kmers <- KmersFromBrightnesses(bright, monomer.med)
  attr(kmers, "mean.intensity") <- mean(mat3d)
  kmers
}

#' @rdname KmersFromImage
#'
#' @param folder.path The path (relative or absolute) to the folder you wish to
#'   process.
#' @param ext the file extension of the images in the folder that you wish to
#'   process (can be rooted in regular expression for extra-safety, as in the
#'   default). You must wish to process all files with this extension; if there
#'   are files that you don't want to process, take them out of the folder.
#' @param out.name The name of the results csv file.
#' @param mcc The number of parallel cores to use for the processing.
#' @export
KmersFromImagesFolder <- function(folder.path = ".", monomer.med,
  tau = NA, mst = NULL, filt = NULL, out.name = "results",
  ext = "\\.tif$", mcc = parallel::detectCores(), verbose = TRUE) {
  init.dir <- getwd()
  on.exit(setwd(init.dir))
  setwd(folder.path)
  tif.names <- list.files(pattern = ext)
  kmerss <- tif.names %>% BiocParallel::bplapply(KmersFromImage,
    tau, monomer.med, verbose = verbose,
    BPPARAM = BiocParallel::MulticoreParam(workers = mcc))
  means <- vapply(kmerss, function(x) attr(x, "mean.intensity"), numeric(1))
  max.ks <- vapply(kmerss, length, integer(1))
  # for each image, the maximum k for which that image contains
  # at least one kmer
  kmers.table <- matrix(0, nrow = length(tif.names), ncol = max(max.ks))
  colnames(kmers.table) <- names(kmerss[[which.max(max.ks)]])
  for (i in seq_along(kmerss)) {
    for (j in seq_along(kmerss[[i]])) {
      kmers.table[i, j] <- kmerss[[i]][j]
    }
  }
  results <- data.frame(ImageName = tif.names, MeanIntensity = means) %>%
    cbind(kmers.table)
  out.name <- filesstrings::MakeExtName(out.name, "csv")
  readr::write_csv(results, out.name)
  invisible(results)
}

#' Get an image where each pixel represents a kmer.
#'
#' Based on a brightness image (array, can be more than two-dimensional), create
#' a kmer image (array) where each pixel represents a kmer (0 for immobile, 1
#' for monomer, 2 for dimer and so on).
#'
#' @param brightness.arr The brightness array.
#' @param monomer.brightness The (median) brightness of a monomer.
#'
#' @return A matrix.
#' @export
KmerArray <- function(brightness.arr, monomer.brightness) {
  stopifnot(monomer.brightness > 1)
  max.b <- max(brightness.arr, na.rm = TRUE)
  if (max.b > monomer.brightness) {
    ranges <- c(seq(1 + 0.5 * (monomer.brightness - 1), max.b,
      monomer.brightness - 1), max.b) %>% unique
    # the unique avoids the unlikely possibility of repeating the max at the end
    ranges.mat <- ranges %>% {
      cbind(dplyr::lag(.), .)[-1, ]
    }
    kmer.arr <- WhichInterval(brightness.arr, ranges.mat)
  } else {
    kmer.arr <- brightness.arr %T>% {
      .[!is.na(.)] <- 0
    }
  }
  kmer.arr
}


#' Create kmer tiff files from brightness csvs
#'
#' For each brightness csv image in a folder, given a monomeric brightness,
#' create a tiff file of the kmer positions using \code{\link{KmerArray}}.
#'
#' @param csv.paths The paths to the brightness csv files, defaults to
#'   \code{list.files(pattern = '[Bb]rightness.*\\.csv')}.
#' @param monomer.brightness The (median) brightness of a monomer.
#' @param out.names The names you want the output files to have (will be forced
#'   to .tif).
#' @param na See \code{\link{WriteIntImage}}.
#' @param mcc The number of parallel cores to use for the processing.
#' @param verbose Do you want to print a message when the function is done?
#'
#' @export
KmerTIFFsFromBrightnessCSVs <- function(monomer.brightness, csv.paths = NULL,
  out.names = NULL, na = "saturate", mcc = parallel::detectCores(),
  verbose = TRUE) {
  if (is.null(csv.paths)) {
    csv.paths <- list.files(pattern = "[Bb]rightness.*\\.csv")
  }
  if (is.null(out.names)) {
    out.names <- stringr::str_replace(csv.paths, "[Bb]rightness",
      "kmers")
  } else {
    if (length(out.names) != length(csv.paths)) {
      stop("The number of input files and output names is not equal.")
    }
  }
  out.names <- vapply(out.names, filesstrings::MakeExtName, character(1),
    "tif", replace = TRUE)
  brightnesses <- BiocParallel::bplapply(csv.paths, ReadImageTxt,
    BPPARAM = BiocParallel::MulticoreParam(workers = mcc))
  kmer.mats <- BiocParallel::bplapply(brightnesses, KmerArray,
    monomer.brightness, BPPARAM = BiocParallel::MulticoreParam(workers = mcc))
  BiocParallel::bpmapply(WriteIntImage, kmer.mats, out.names,
    na = na, BPPARAM = BiocParallel::MulticoreParam(workers = mcc)) %>%
    invisible
  if (verbose)
    print("Done. Please check folder.")
}
