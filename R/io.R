#' Read image as array object.
#'
#' Read in an image file from the disk as an array of pixel intensities. Give it
#' an attribute "bits" detailing how many bits per sample there were in the file
#' on disk. optionally allows the user to read in pixel intensities as integers
#' in the range [0, 2 ^ bits - 1]. See "Details".
#'
#' This function wraps \code{\link[EBImage]{readImage}} and
#' \code{\link[EBImage]{imageData}} from the \code{\link[EBImage]{EBImage}}
#' package. By default, \code{readImage} reads in pixel intensities in the range
#' [0, 1]. \code{ReadImageData} optionally allows the user to read in pixel
#' intensities as integers in the range [0, 2 ^ bits - 1], as would be seen when
#' viewing the file in an application such as imagej. This is sometimes
#' necessary, for example, when calculating number and brightness, where we need
#' pixel values to be in units of "counts". This functionality is enabled by
#' default, but can be disabled with \code{restore.counts = FALSE}.
#'
#' This function deals with reading in image stacks in the same way as
#' \code{\link[EBImage]{EBImage}}.
#'
#' @param image.name The path to the image file on disk. The file extension must
#'   be one of ".jpeg", ".png", ".tiff" or ".tif".
#' @param restore.counts If this is \code{TRUE} (the default), pixel intensities
#'   are read in as integers in the range [0, 2 ^ bits - 1]. If set to
#'   \code{FALSE}, the default behaviour of \code{\link[EBImage]{readImage}} is
#'   restored, however the resulting array still has a "bits" attribute, which
#'   it wouldn't with \code{\link[EBImage]{readImage}}.
#'
#' @return An array with a "bits" attribute.
ReadImageData <- function(image.name, restore.counts = TRUE) {
  image.data <- suppressWarnings(EBImage::imageData(
    EBImage::readImage(image.name)))
  file.info <- paste("file", image.name) %>% system(intern = TRUE)
  bits <- stringr::str_split_fixed(file.info, "bps", 2)[2] %>% filesstrings::NthNumber(1)
  if (restore.counts) {
    if (!bits %in% c(8, 16)) {
      stop("restore.counts only works with 8- or 16-bit images.")
    }
    image.data <- image.data * (2 ^ bits - 1)
  }
  attr(image.data, "bits") <- bits
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
#' several channels).
#'
#' @param img.arr An image, represented by a 2- or 3-dimensional array.
#' @param file.name The name you wish to associate with the output files,
#'   \emph{without} a file extension.
#'
#' @export
WriteImageTxt <- function(img.arr, file.name) {
  d <- dim(img.arr)
  nd <- length(d)
  if (! nd %in% c(2, 3)) {
    stop("img.arr must be 2- or 3-dimensional")
  }
  if (nd == 2) {
    as.data.frame(img.arr) %>%
      readr::write_csv(filesstrings::MakeExtName(file.name, "csv"),
                       col_names = FALSE)
  } else {
    slices.dfs <- lapply(seq_len(d[3]), Slices, img.arr) %>%
      lapply(as.data.frame)
    file.names <- paste0(file.name, seq_len(d[3])) %>%
      sapply(filesstrings::MakeExtName, "csv") %>%
      (filesstrings::NiceNums)
    mapply(readr::write_csv, slices.dfs, file.names, col_names = FALSE)
  }
}
