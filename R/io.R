#' Read image as array object.
#'
#' Read in an image file from the disk as an array of pixel intensities.
#'
#' This function wraps [EBImage::readImage()] and [EBImage::imageData()]. By
#' default, `readImage` reads in pixel intensities in the range \code{[0, 1]}.
#' `ReadImageData` reads in pixel intensities as integers as they would be
#' represented in a tiff file and displayed therefrom in ImageJ. This is
#' necessary when calculating number and brightness, where we need pixel values
#' to be in units of 'counts'.
#'
#' Thinking of the read image as a matrix `mat`, the pixel at \eqn{x = }`i`,
#' \eqn{y = }`j` has colour based on the value of \code{mat[i, j]} where the
#' \eqn{x} axis points right and the \eqn{y} axis points down. This is in
#' accordance with how [EBImage::EBImage]'s [EBImage::readImage()] (which this
#' function wraps). However, when one prints a matrix in a console (or views it
#' in a program such as excel), the value in position \eqn{x = }`i`, \eqn{y =
#' }`j` is from `mat[j, i]`, so if you're confused about a phantom
#' transposition, this is why.
#'
#' Sometimes (for whatever reason) the reading of image values as integers with
#' `readTIFF(..., as.is = TRUE)` fails such that the values which should be 1,
#' 2, 3 etc. end up as 256, 512, 778 etc. This function corrects for this as
#' follows: if the greatest common divisor of the elements in the array are one
#' of 2^8, 2^12, 2^16 or 2^32, then divide each element of the array by this
#' greatest common divisor. This will return 256, 512, 778 etc. to 1, 2, 3 etc.
#'
#' @param image.name The path to the image file on disk. The file extension must
#'   be one of '.jpeg', '.png', '.tiff' or '.tif'.
#' @param fix.lut When reading in images (via [EBImage::readImage()]), R can
#'   give an array of different dimensionality than you expect. If you suspect
#'   this happening, set the value of this parameter to the \emph{number of
#'   dimensions} that you expect your read image to have and this function will
#'   try to automatically give you the image array in the form you want. Read
#'   [FixLUTError()] to find out more.
#'
#' @return An array of integers representing the image.
#'
#' @examples
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#'
#' @export
ReadImageData <- function(image.name, fix.lut = NULL) {
  EBIrid <- function(path) {
    suppressWarnings(EBImage::imageData(
      EBImage::readImage(path, as.is = TRUE)))
  }
  image.data <- EBIrid(image.name)
  if (!is.null(fix.lut)) {
    if (isTRUE(fix.lut)) {
      stop("If fix.lut is not set to false, it must be specified as an integer",
           " (not as TRUE). Read the documentation for ReadImageData.")
    }
    image.data <- FixLUTError(image.data, fix.lut)
  }
  vec <- as.vector(image.data)
  if (any(vec > 0)) {
    g0unq <- as.vector(image.data) %>%
      {.[. > 0]} %>%
      unique()
    if (length(g0unq) == 1) {
      gcdg0 <- g0unq
    } else {
      gcdg0 <- DescTools::GCD(g0unq)
    }
    if (gcdg0 %in% (2 ^ c(8, 12, 16, 32))) {
      image.data <- image.data / gcdg0
    }
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
#' displaying an image with [EBImage::display()] in R will yield the
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
#' setwd(tempdir())
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' WriteImageTxt(img, 'temp')
#'
#' @export
WriteImageTxt <- function(img.arr, file.name) {
  d <- dim(img.arr)
  nd <- length(d)
  if (!nd %in% c(2, 3, 4)) {
    stop("img.arr must be 2-, 3- or 4-dimensional")
  }
  if (nd == 2) {
    img.arr <- t(img.arr)
    as.data.frame(img.arr) %>%
      readr::write_csv(filesstrings::give_ext(file.name, "csv"),
                       col_names = FALSE)
  } else if (nd == 3) {
    img.arr <- aperm(img.arr, c(2, 1, 3))
    slices.dfs <- lapply(seq_len(d[3]), Slices, img.arr) %>%
      lapply(as.data.frame)
    file.names <- paste0(file.name, "_", seq_len(d[3])) %>%
      vapply(filesstrings::give_ext, character(1), "csv") %>%
      (filesstrings::nice_nums)
    mapply(readr::write_csv, slices.dfs, file.names, col_names = FALSE,
           SIMPLIFY = FALSE) %>%
      invisible
  } else {
    img.arr <- aperm(img.arr, c(2, 1, 3, 4))
    grid <- expand.grid(seq_len(d[3]), seq_len(d[4]))
    file.names <- paste0(file.name, "_", grid[, 1], "_", grid[, 2]) %>%
      vapply(filesstrings::give_ext, character(1), "csv") %>%
      (filesstrings::nice_nums)
    dfs <- purrr::map(Mat2RowList(as.matrix(grid)),
                       ~ img.arr[, , .[1], .[2]]) %>%
      lapply(as.data.frame)
    mapply(readr::write_csv, dfs, file.names, col_names = FALSE,
           SIMPLIFY = FALSE) %>%
      invisible
  }
}

#' @rdname WriteImageTxt
#'
#' @examples
#' img <- ReadImageTxt('temp_01.csv')
#' file.remove(list.files())  # cleanup
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
#' [EBImage::writeImage] truncates all values above 1 to 1 and all below 0 to 0
#' when writing images. This function allows you to write integer-vlued arrays
#' to disk as tiff files as you would want.
#'
#' @param img.arr An integer array.
#' @param file.name The name of the tif file (with or without the .tif) that you
#'   wish to write.
#' @param na How do you want to treat `NA` values? R can only write integer
#'   values (and hence not `NA`s) to tiff pixels. `na = 'saturate'` sets them to
#'   saturated value. `na = 'zero'` sets them to zero, while `na = 'error'` will
#'   give an error if the image contains `NA`s. You can also specify directly
#'   here a natural number (must be between 0 and 2 ^ 16 - 1) to use in place of
#'   `NA`s.
#'
#' @examples
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' setwd(tempdir())
#' WriteIntImage(img, '50.tif')
#' file.remove(list.files())  # cleanup
#'
#' @export
WriteIntImage <- function(img.arr, file.name, na = "error") {
  to.invisibly.return <- img.arr
  if (!all(CanBeInteger(img.arr), na.rm = TRUE)) {
    stop("img.arr must contain only integers")
  }
  if (!all(img.arr >= 0, na.rm = TRUE)) {
    stop("img.arr must not contain negative values")
  }
  stopifnot(length(na) == 1)
  if (is.numeric(na)) {
    stopifnot(na >= 0, na < 2 ^ 16)
    na <- floor(na)
  } else {
    na <- RSAGA::match.arg.ext(na, c("saturate", "zero", "error"),
                               ignore.case = TRUE)
  }
  any.nas <- anyNA(img.arr)
  if (na == "error" && any.nas) stop("img.arr contains NA values.")
  if (is.numeric(na)) img.arr[is.na(img.arr)] <- na
  mx <- max(img.arr, na.rm = TRUE)
  if (mx >= 2 ^ 16) {
    stop("The maximum value in img.arr must be less than 2^16")
  } else if (mx >= 2 ^ 8) {
    bits.per.sample <- 16
  } else {
    bits.per.sample <- 8
  }
  img.arr <- img.arr / (2 ^ bits.per.sample - 1)
  if (any.nas) {
    if (na == "saturate") img.arr[is.na(img.arr)] <- 1
    if (na == "zero") img.arr[is.na(img.arr)] <- 0
  }
  file.name <- filesstrings::give_ext(file.name, "tif")
  message("Now writing ", file.name, ".")
  EBImage::writeImage(img.arr, file.name, bits.per.sample = bits.per.sample)
  invisible(to.invisibly.return)
}

#' Fix an image that didn't recognise channels while reading
#'
#' Sometimes, when you open an image in ImageJ, it displays the channels as you
#' would like, but when you read it into R, it has just mashed all the channels
#' (which you would like to be separated somehow) into a stack. In my
#' expreience, it always does so in a way that, say you have a stack of 3
#' channels and 5 z positions, then the red images would occupy `[, , 1]`, `[, ,
#' 6]` and `[, , 11]`. This function fixes this kind of confusion. So, in that
#' example it would have a 3d array as input and a 4d as output with dimensions
#' (assuming our images are 256x256 pixels) `256, 256, 3, 5`.
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
#' setwd(tempdir())
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' WriteIntImage(img[, , 1], '50_1.tif')
#' WriteIntImage(img[, , 2], '50_2.tif')
#' Stack2DTifs(c('50_1.tif', '50_2.tif'), '50_1_2')
#' file.remove(list.files())
#' @export
Stack2DTifs <- function(file.names, out.name, mcc = parallel::detectCores()) {
  out.name <- filesstrings::give_ext(out.name, "tif")
  images <- BiocParallel::bplapply(file.names, EBImage::readImage,
    BPPARAM = bpp(mcc))
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

#' Convert a binary image to tiff.
#'
#' Read a binary image (pixels must be integers), specifying its dimension and
#' write it back to disk as a tifff file. The `Bin2TiffFolder` function does
#' this for an entire folder.
#'
#' When using `Bin2TiffFolder`, the images must all have the same dimension.
#'
#' @param bin.path The path to the binary file.
#' @param bits How many bits are there per pixel? This should be a multiple of
#'   8.
#' @param height The height of the image (in pixels).
#' @param width The width of the image (in pixels).
#' @param frames How many frames are there in the stack (default 1).
#' @param endian Eigher `"big"` or `"little"` (see [readBin]).
#' @param signed Logical. Only used for integers of sizes 1 and 2, when it
#'   determines if the quantity on file should be regarded as a signed or
#'   unsigned integer.
#'
#' @export
Bin2Tiff <- function(bin.path, bits, height, width, frames = 1,
                     endian = .Platform$endian, signed = TRUE) {
  file.name.short <- bin.path
  if (stringr::str_detect(file.name.short, "/")) {
    file.dir <- filesstrings::str_before_nth(file.name.short, "/", -1)
    current.wd <- getwd()
    on.exit(setwd(current.wd))
    setwd(file.dir)
    file.name.short <- stringr::str_split(file.name.short, "/") %>%
      unlist %>%
      dplyr::last()
  }
  bin.vector <- readBin(file.name.short, "int", n = height * width * frames,
                        signed = signed, size = bits %/% 8)
  arr <- array(bin.vector, dim = c(height, width, frames))
  if (dim(arr)[3] == 1) arr <- arr[, , 1]
  WriteIntImage(arr, filesstrings::before_last_dot(file.name.short))
}

#' @rdname Bin2Tiff
#'
#' @param folder.path The path to the folder (defaults to current location).
#'
#' @examples
#' setwd(tempdir())
#' dir.create("temp_dir")
#' file.copy(system.file("extdata", "b.bin", package = "nandb"), "temp_dir")
#' Bin2Tiff("temp_dir/b.bin", height = 2, width = 2, bits = 8)
#' Bin2TiffFolder("temp_dir", height = 2, width = 2, bits = 8)
#' list.files("temp_dir")
#' @export
Bin2TiffFolder <- function(folder.path = ".", bits, height, width, frames = 1,
                           endian = .Platform$endian, signed = TRUE) {
  current.wd <- getwd()
  on.exit(setwd(current.wd))
  setwd(folder.path)
  bin.names <- list.files(pattern = "\\.bin$")
  lapply(bin.names, Bin2Tiff, bits, height, width, frames, endian, signed)
}
