#' TIFF I/O.
#'
#' These functions are imported from the `detrendr` package. See their full
#' documentation at [detrendr::read_tif()] and [detrendr::write_tif()].
#'
#' @inheritParams detrendr::read_tif
#' @inheritParams detrendr::write_tif
#'
#' @name tiff-io
NULL

#' @rdname tiff-io
#' @export
read_tif <- function(path, n_ch = 1) {
  detrendr::read_tif(path, n_ch = n_ch)
}

#' @rdname tiff-io
#' @export
write_tif <- function(img, file_name, na = "error", rds = FALSE) {
  detrendr::write_tif(img, file_name, na = na, rds = rds)
}

#' Image display.
#'
#' This function is imported from the `detrendr` package. See the full
#' documentation at [detrendr::display()].
#'
#' @inheritParams detrendr::display
#'
#' @export
display <- function(...) detrendr::display(...)


#' Read/write an image array to/from disk as text file(s).
#'
#' If (as with brightness) you wish for the pixel values in an image to be
#' represented by real numbers that aren't necessarily integers, the tiff format
#' won't work. As a workaround we represent images (arrays) as tab-separated
#' `.txt` files on disk, where for image stacks (3-dimensional arrays), we
#' write one file for each slice, numbering it with the slice number.
#'
#' @param img An image, represented by a 2- or 3-dimensional array.
#' @param path The name of the input/output output file(s), \emph{without}
#'   a file extension.
#' @param rds In addition to writing a text file, save the image as an RDS (a
#'   single R object) file?
#'
#' @name text-image-io
NULL

#' @rdname text-image-io
#'
#' @examples
#' setwd(tempdir())
#' img <- read_tif(system.file('extdata', '50.tif', package = 'nandb'))
#' write_txt_img(img, 'temp')
#'
#' @export
write_txt_img <- function(img, path, rds = FALSE) {
  d <- dim(img)
  nd <- length(d)
  if (!nd %in% c(2, 3, 4)) {
    stop("img must be 2-, 3- or 4-dimensional")
  }
  if (rds) saveRDS(img, file = filesstrings::give_ext(path, "rds"))
  if (nd == 2) {
    as.data.frame(img) %>%
      readr::write_tsv(path = filesstrings::give_ext(path, "txt"),
                       col_names = FALSE)
  } else if (nd == 3) {
    slices_dfs <- purrr::map(seq_len(d[3]), ~ img[, , .]) %>%
      purrr::map(as.data.frame)
    paths <- paste0(path, "_", seq_len(d[3])) %>%
      purrr::map_chr(filesstrings::give_ext, "txt") %>%
      filesstrings::nice_nums()
    purrr::map2(slices_dfs, paths,
                ~ readr::write_tsv(.x, .y, col_names = FALSE)) %>%
      invisible()
  } else {
    grid <- expand.grid(seq_len(d[3]), seq_len(d[4])) %>% as.matrix()
    paths <- paste0(path, "_ch", grid[, 1], "_", grid[, 2]) %>%
      purrr::map_chr(filesstrings::give_ext, "txt") %>%
      filesstrings::nice_nums()
    dfs <- purrr::map(enlist_rows(grid), ~ img[, , .[1], .[2]]) %>%
      purrr::map(as.data.frame)
    purrr::map2(dfs, paths, ~ readr::write_tsv(.x, .y, col_names = FALSE))
  }
  invisible(img)
}

#' @rdname text-image-io
#'
#' @examples
#' img <- read_txt_img('temp_01.txt')
#' file.remove(list.files())  # cleanup
#' @export
read_txt_img <- function(path) {
  suppressMessages(readr::read_tsv(path, col_names = FALSE,
                                   progress = FALSE)) %>%
    data.matrix() %>%
    magrittr::set_colnames(value = NULL)
}

nb_get_img <- function(img, n_ch, min_d, max_d) {
  if (is.character(img)) {
    if (length(img) != 1) {
      stop("If 'img' is specified as a character vector ",
           "(file name), it must be length 1")
    }
    img %<>% read_tif(n_ch = n_ch)
  }
  checkmate::assert_array(img, min.d = min_d, max.d = max_d)
  d <- dim(img)
  if (length(d) == 3 && n_ch != 1) {
    stop("You have specified 'img' as having ", n_ch, " channels, ",
         "but it looks like it has 1 channel.")
  }
  if ((length(d) == 4) && (d[3] != n_ch)) {
    stop("You have specified 'img' as having ", n_ch, " channels, ",
         "but it looks like it has ", d[3], ".")
  }
  img
}
