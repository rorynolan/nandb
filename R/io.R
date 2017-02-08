#' Read image as array object.
#'
#' Read in an image file from the disk as an array of pixel intensities.
#'
#' This function wraps \code{\link[EBImage]{readImage}} and
#' \code{\link[EBImage]{imageData}} from the \code{\link[EBImage]{EBImage}}
#' package. By default, \code{readImage} reads in pixel intensities in the range
#' [0, 1]. \code{ReadImageData} reads in pixel intensities as integers as they
#' would be represented in a tiff file and displayed therefrom in ImageJ. This
#' is necessary when calculating number and brightness, where we need pixel
#' values to be in units of 'counts'.
#'
#' Thinking of the read image as a matrix \code{mat}, the pixel at \eqn{x =
#' }\code{i}, \eqn{y = }\code{j} has colour based on the value of \code{mat[i,
#' j]} where the \eqn{x} axis points right and the \eqn{y} axis points down.
#' This is in accordance with how \code{\link[EBImage]{EBImage}}'s
#' \code{\link[EBImage]{readImage}} (which this function wraps). However, when
#' one prints a matrix in a console (or views it in a program such as excel),
#' the value in position \eqn{x = }\code{i}, \eqn{y = }\code{j} is from
#' \code{mat[j, i]}, so if you're confused about a phantom transposition, this
#' is why.
#'
#' @param image.name The path to the image file on disk. The file extension must
#'   be one of '.jpeg', '.png', '.tiff' or '.tif'.
#' @param fix.lut When reading in images (via \code{\link[EBImage]{readImage}}),
#'   R can give an array of different dimensionality than you expect. If you
#'   suspect this happening, set the value of this parameter to the \emph{number
#'   of dimensions} that you expect your read image to have and this function
#'   will try to automatically give you the image array in the form you want.
#'   Read \code{\link{FixLUTError}} to find out more.
#'
#' @return An array of integers representing the image.
#'
#' @examples
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#'
#' @export
ReadImageData <- function(image.name, fix.lut = NULL) {
  image.data <- suppressWarnings(EBImage::imageData(
    EBImage::readImage(image.name, as.is = TRUE)))
  if (!is.null(fix.lut)) {
    if (isTRUE(fix.lut)) {
      stop("If fix.lut is not set to false, it must be specified as an integer",
           " (not as TRUE). Read the documentation for ReadImageData.")
    }
    image.data <- FixLUTError(image.data, fix.lut)
  }
  if (!all(CanBeInteger(image.data))) {
    stop("Failed to read in the image as an array of integers.")
  }
  image.data
}

#' Read/write an image array to/from disk as text file(s).
#'
#' If (as with brightness) you wish for the pixel values in an image to be
#' represented by real numbers that aren't necessarily integers, the tiff format
#' won't work. As a workaround we represent images (arrays) as comma-separated
#' value (csv) files on disk, where for image stacks (3-dimensional arrays), we
#' write one file for each slice, numbering it with the slice number.
#'
#' The image slices are transposed prior to being written to disk to ensure that
#' displaying an image with \code{\link[EBImage]{display}} in R will yield the
#' same result (as opposed to a transposed image) as displaying the written text
#' file in ImageJ (i.e. I've made a modification to ensure the files display
#' correctly in ImageJ). These functions do not work for 4-dimensional arrays
#' (e.g. a z stack with several channels).
#'
#' @param img.arr An image, represented by a 2- or 3-dimensional array.
#' @param file.name The name of the input/output output file(s), \emph{without}
#'   a file extension.
#'
#' @examples
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' WriteImageTxt(img, 'temp')
#'
#' @export
WriteImageTxt <- function(img.arr, file.name) {
  d <- dim(img.arr)
  nd <- length(d)
  if (!nd %in% c(2, 3)) {
    stop("img.arr must be 2- or 3-dimensional")
  }
  if (nd == 2) {
    img.arr <- t(img.arr)
    as.data.frame(img.arr) %>%
      readr::write_csv(filesstrings::MakeExtName(file.name, "csv"),
                       col_names = FALSE)
  } else {
    img.arr <- aperm(img.arr, c(2, 1, 3))
    slices.dfs <- lapply(seq_len(d[3]), Slices, img.arr) %>%
      lapply(as.data.frame)
    file.names <- paste0(file.name, "_", seq_len(d[3])) %>%
      vapply(filesstrings::MakeExtName, character(1), "csv") %>%
      (filesstrings::NiceNums)
    mapply(readr::write_csv, slices.dfs, file.names, col_names = FALSE) %>%
      invisible
  }
}

#' @rdname WriteImageTxt
#'
#' @examples
#' img <- ReadImageTxt('temp_01.csv')
#' file.remove(list.files(pattern = '^temp.*\\.csv$'))
#' @export
ReadImageTxt <- function(file.name) {
  suppressMessages(readr::read_csv(file.name, col_names = FALSE,
    progress = FALSE)) %>%
    data.matrix %>%
    magrittr::set_colnames(value = NULL) %>%
    t
}

#' Write an integer array to disk as a tiff image.
#'
#' [EBImage][EBImage::EBImage]'s [writeImage][EBImage::writeImage] truncates all
#' values above 1 to 1 and all below 0 to 0 when writing images. This function
#' allows you to write integer-vlued arrays to disk as tiff files as you would
#' want.
#'
#' @param img.arr An integer array.
#' @param file.name The name of the tif file (with or without the .tif) that you
#'   wish to write.
#' @param na How do you want to treat \code{NA} values? R can only write integer
#'   values (and hence not \code{NA}s) to tiff pixels. \code{na = 'saturate'}
#'   sets them to saturated value. \code{na = 'zero'} sets them to zero, while
#'   \code{na = 'error'} will give an error if th image contains \code{NA}s.
#'
#' @examples
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' dir.create('tempdir')
#' WriteIntImage(img, 'tempdir/50.tif')
#'
#' @export
WriteIntImage <- function(img.arr, file.name, na = "error") {
  if (!all(CanBeInteger(img.arr), na.rm = TRUE)) {
    stop("img.arr must contain only integers")
  }
  if (!all(img.arr >= 0, na.rm = TRUE)) {
    stop("img.arr must not contain negative values")
  }
  na <- RSAGA::match.arg.ext(na, c("saturate", "zero", "error"),
    ignore.case = TRUE)
  any.nas <- anyNA(img.arr)
  if (na == "error" && any.nas)
    stop("img.arr contains NA values.")
  mx <- max(img.arr, na.rm = TRUE)
  if (mx >= 2^32) {
    stop("The maximum value in img.arr must be less than 2^32")
  } else if (mx >= 2^16) {
    bits.per.sample <- 32
  } else if (mx >= 2^8) {
    bits.per.sample <- 16
  } else {
    bits.per.sample <- 8
  }
  img.arr <- img.arr/(2^bits.per.sample - 1)
  if (any.nas) {
    if (na == "saturate")
      img.arr[is.na(img.arr)] <- 1
    if (na == "zero")
      img.arr[is.na(img.arr)] <- 0
  }
  file.name <- filesstrings::MakeExtName(file.name, "tif")
  EBImage::writeImage(img.arr, file.name, bits.per.sample = bits.per.sample)
}

#' Fix an image that didn't recognise channels while reading
#'
#' Sometimes, when you open an image in ImageJ, it displays the channels as you
#' would like, but when you read it into R, it has just mashed all the channels
#' (which you would like to be separated somehow) into a stack. In my
#' expreience, it always does so in a way that, say you have a stack of 3
#' channels and 5 z positions, then the red images would occupy \code{[, , 1]},
#' \code{[, , 6]} and \code{[, , 11]}. This is the type of fixing that this
#' function performs. So in that example it would have a 3d array as input and a
#' 4d as output with dimensions (assuming our images are 256x256 pixels)
#' \code{256, 256, 3, 5}.
#'
#' @param img.arr An array, the read image.
#' @param n.ch The number of channels that you want the read image to have.
#'
#' @return An array, channel indices in the third slot, slice indices in the
#'   fourth.
#'
#' @examples
#' library(magrittr)
#' x <- lapply(1:300, function(x) matrix(runif(4), nrow = 2)) %>%
#'   Reduce(function(x, y) abind::abind(x, y, along = 3), .)
#' str(x)
#' ForceChannels(x, 6) %>% str
#'
#' @export
ForceChannels <- function(img.arr, n.ch) {
  d <- dim(img.arr)
  if (length(d) != 3)
    stop("img.arr must be a 3-dimensional array")
  if (d[3]%%n.ch != 0) {
    message <- paste0("The number of slices in the input array must be a ",
      "multiple of n.ch, however you have ",
      d[3], " slices in your input array and n.ch = ", n.ch, ".")
    stop(message)
  }
  index.groups <- lapply(seq_len(n.ch), function(x) seq(x,
    d[3], n.ch))
  ch.groups <- lapply(index.groups, function(x) img.arr[, ,
    x])
  to.be.apermed <- Reduce(function(x, y) abind::abind(x, y,
    along = 4), ch.groups)
  aperm(to.be.apermed, c(1, 2, 4, 3))
}

#' Put individual tif files into one tif stack.
#'
#' Say you have saved what you would like to be one 3D tif stack as a series of
#' 2D tif files. This helps you stack them into one file.
#'
#' @param file.names The names of the files to stack (in order).
#' @param out.name The name of the output .tif file.
#' @param mcc The number of parallel cores to use for the processing.
#'
#' @examples
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' WriteIntImage(img[, , 1], '50_1.tif')
#' WriteIntImage(img[, , 2], '50_2.tif')
#' Stack2DTifs(c('50_1.tif', '50_2.tif'), '50_1_2')
#' file.remove(c('50_1.tif', '50_2.tif', '50_1_2.tif'))
#' @export
Stack2DTifs <- function(file.names, out.name, mcc = parallel::detectCores()) {
  out.name <- filesstrings::MakeExtName(out.name, "tif")
  images <- BiocParallel::bplapply(file.names, EBImage::readImage,
    BPPARAM = BiocParallel::MulticoreParam(workers = mcc))
  dims <- lapply(images, dim)
  u <- unique(dims)
  if (length(u) != 1) {
    stop("All the images to be stacked must have the same dimensions")
  }
  if (length(u[[1]]) != 2) {
    stop("The images must all be 2-dimensional")
  }
  stacked <- Reduce(function(x, y) abind::abind(x, y, along = 3),
    images)
  EBImage::writeImage(stacked, out.name)
}

