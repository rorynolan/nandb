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
#' values to be in units of "counts".
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
#'   be one of ".jpeg", ".png", ".tiff" or ".tif".
#' @param fix.lut When reading in images (via \code{\link[EBImage]{readImage}}),
#'   R can give an array of different dimensionality than you expect. If you
#'   suspect this happening, set the value of this parameter to the \emph{number
#'   of dimensions} that you expect your read image to have and this function
#'   will try to automatically give you the image array in the form you want.
#'   Read \code{\link{FixLUTError}} to find out more.
#'
#' @return An array of integers representing the image.
#' @export
ReadImageData <- function(image.name,
                          fix.lut = NULL) {
  image.data <- suppressWarnings(EBImage::imageData(
    EBImage::readImage(image.name, as.is = TRUE)))
  if (!is.null(fix.lut)) {
    if (isTRUE(fix.lut)) {
      stop("If fix.lut is not set to false, it must be specified as an integer (not as TRUE). Read the documentation for ReadImageData.")
    }
    image.data <- FixLUTError(image.data, fix.lut)
  }
  if (!all(CanBeInteger(image.data))) {
    stop("Failed to read in the image as an array of integers.")
  }
  image.data
}

#' Write an image array to disk as text file(s)
#'
#' If (as with brightness) you wish for the pixel values in an image to be
#' represented by real numbers that aren't necessarily integers, the tiff format
#' won't work. As a workaround we write images (arrays) as comma-separated value
#' (csv) files, where for image stacks (3-dimensional arrays), we write one file
#' for each slice, numbering it with the slice number.
#'
#' This function does not yet work for 4-dimensional arrays (e.g. a z stack with
#' several channels). The image slices are transposed prior to being written to
#' disk to ensure that displaying an image with \code{\link[EBImage]{display}}
#' in R will yield the same result (as opposed to a transposed image) as
#' displaying the written text file in ImageJ (i.e. I've made a modification to
#' ensure the files display correctly in ImageJ).
#'
#' @param img.arr An image, represented by a 2- or 3-dimensional array.
#' @param file.name The name you wish to associate with the output files,
#'   \emph{without} a file extension.
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
    file.names <- paste0(file.name, seq_len(d[3])) %>%
      sapply(filesstrings::MakeExtName, "csv") %>%
      (filesstrings::NiceNums)
    mapply(readr::write_csv, slices.dfs, file.names, col_names = FALSE)
  }
}
