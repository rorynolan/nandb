nb_get_img <- function(img) {
  if (is.character(img)) {
    if (length(img) != 1) {
      stop("If 'img' is specified as a character vector ",
           "(i.e. a file name), it must be length 1")
    }
    img %<>% ijtiff::read_tif()
  }
  checkmate::assert_array(img, min.d = 3, max.d = 4)
  checkmate::assert_numeric(img, lower = 0)
  if (!isTRUE(all.equal(img, floor(img), check.attributes = FALSE))) {
    stop("`img` must be positive integers (and NAs) only")
  }
  img
}
