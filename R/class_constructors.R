#' Number and brightness image classes.
#'
#' The `number_img` and `brightness_img` classes are designed to hold objects
#' which are images calculated from the *number and brightness* technique.
#'
#' An object of class `number_img` or `brightness_img` is a matrix of real
#' numbers with 4 attributes: \describe{\item{`def`}{Are we using the `"N"` or
#' `"n"` definition of number, or the `"B"` or `"epsilon"` definition of
#' brightness?} \item{`thresh`}{A positive integer, possibly an object of class
#' [autothresholdr::th] detailing which threshold and thresholding method was
#' used in preprocessing (in the multi-channel case, one threshold per channel
#' is given).} \item{`tau`}{A positive number indicating the tau parameter used
#' for detrending with an attribute `auto` which is a logical indicating whether
#' or not the tau parameter was chosen automatically.}\item{`filt`}{Was mean or
#' median filtering used in postprocessing?}}
#'
#' @param img The calculated number or brightness image.
#' @param def The number or brightness definition used.
#' @param thresh A positive integer, possibly an object of class
#'   [autothresholdr::th]. If the different channels of the image had different
#'   thresholds, this argument may be specified as a vector or list (of positive
#'   integers, possibly objects of class [autothresholdr::th]), one element for
#'   each channel.
#' @param tau A positive number with an attribute `auto`. If the different
#'   channels of the image had different `tau`s, this argument may be specified
#'   as a list(of positive numbers with attributes `auto`), one element for each
#'   channel.
#' @param filt A string, the filtering method used. Must be either `"mean"` or
#'   `"median"`, or `NA` for no filtering. If the different channels of the
#'   image had different filters, this may be specified as a character vector,
#'   one element for each channel.
#'
#' @return An object of class `number`.
#'
#' @name nb-img-classes
NULL

#' @rdname nb-img-classes
#' @export
number_img <- function(img, def, thresh = NA, tau = NA, filt = NA) {
  checkmate::assert_array(img, min.d = 2, max.d = 3)
  n_ch <- dim(img) %>% {dplyr::if_else(length(.) == 3, .[3], 1L)}
  img %<>% number_img_common(n_ch = n_ch, def = def,
                         thresh = thresh, tau = tau, filt = filt)
  class(img) %<>% c("number_img", .)
  img
}

number_img_common <- function(img, n_ch, def, thresh = NA, tau = NA, filt = NA) {
  checkmate::assert_string(def)
  if (!isTRUE(def %in% c("n", "N"))) stop("def must be one of 'n' or 'N'.")
  if (length(thresh) == 1 && n_ch > 1) thresh %<>% {rep(list(.), n_ch)}
  if (length(tau) == 1 && n_ch > 1) tau %<>% {rep(list(.), n_ch)}
  if (all(is.na(filt))) {
    filt <- rep(NA, n_ch)
  } else {
    filt %<>% fix_filt()
  }
  if (length(filt) == 1 && n_ch > 1) filt %<>% rep(n_ch)
  thresh %<>% c_list_attr_na()
  tau %<>% c_list_attr_na()
  checkmate::assert_numeric(thresh)
  for (i in seq_along(tau)) {
    if (is.na(tau[[i]]) && (! "auto" %in% names(attributes(tau[[i]])))) {
      attr(tau[[i]], "auto") <- FALSE
    }
  }
  checkmate::assert_numeric(tau, lower = 0)
  if (! "auto" %in% names(attributes(tau))) {
    stop("If tau is specified, it must have an attribute 'auto'.")
  } else if (length(tau) != length(attr(tau, "auto"))) {
    if (length(attr(tau, "auto")) == 1) {
      attr(tau, "auto") %<>% rep(length(tau))
    } else {
      stop("The 'auto' attribute of 'tau' ",
           "must be the same length as 'tau' itself.")
    }
  }
  checkmate::assert_logical(attr(tau, "auto"))
  if ((! all(is.na(tau))) && anyNA(attr(tau, "auto"))) {
    stop("Each element of 'tau' must have an associated attribute 'auto' ",
         "which must be TRUE or FALSE and not NA.")
  }
  if (!filesstrings::all_equal(c(length(thresh), length(tau),
                                 length(filt), n_ch))) {
    stop("The lengths of 'thresh', 'tau' and 'filt' must all be the same ",
         "as the number of channel in 'img'.")
  }
  for (att in c("def", "thresh", "tau", "filt")) attr(img, att) <- get(att)
  img
}

#' @rdname nb-img-classes
#' @export
brightness_img <- function(img, def, thresh = NA, tau = NA, filt = NA) {
  checkmate::assert_string(def)
  def %<>% tolower()
  if (def == "b") {
    def <- "B"
  } else if (startsWith("epsilon", def)) {
    def <- "epsilon"
  } else {
    stop("'def' must be one of 'B' or 'epsilon'.")
  }
  if (!isTRUE(def %in% c("B", "epsilon")))
    stop("def must be one of 'B' or 'epsilon'.")
  out <- number_img(img, "n", thresh, tau, filt)
  attr(out, "def") <- def
  class(out)[class(out) == "number_img"] <- "brightness_img"
  out
}

c_list_attr_na <- function(x) {
  l <- length(x)
  if (is.list(x)) {
    x_attr_names <- x %>%
      purrr::map(~ names(attributes(.))) %>%
      unlist() %>%
      unique()
    for (i in seq_along(x)) {
      attr(x[[i]], "class") <- class(x[[i]])[1]
      for (name in x_attr_names) {
        if (! name %in% names(attributes(x[[i]]))) {
          attr(x[[i]], name) <- NA
        }
      }
    }
    atts <- x %>%
      purrr::map(~ attributes(.)) %>%
      dplyr::bind_rows()
    atts$class <- NULL
    x %<>% purrr::reduce(c)
    attributes(x) <- atts
  }
  stopifnot(length(x) == l)
  x
}

#' Number and brightness time series image classes.
#'
#' The `number_ts_img` and `brightness_ts_img` classes are designed to hold
#' objects which are images calculated from the *number and brightness*
#' technique.
#'
#' An object of class `number_ts_img` or `brightness_ts_img` is a 3- or
#' 4-dimensional array of real numbers with 4 attributes:
#' \describe{\item{`def`}{Are we using the `"N"` or `"n"` definition of number,
#' or the `"B"` or `"epsilon"` definition of brightness?} \item{`thresh`}{A
#' positive integer, possibly an object of class [autothresholdr::th] detailing
#' which threshold and thresholding method was used in preprocessing (in the
#' multi-channel case, one threshold per channel is given).} \item{`tau`}{A
#' positive number indicating the tau parameter used for detrending with an
#' attribute `auto` which is a logical indicating whether or not the tau
#' parameter was chosen automatically (in the multi-channel case, one `tau` per
#' channel is given).} \item{`frames_per_set`}{A positive integer detailing how
#' many frames were used in the calculation of each point in the number or
#' brightness time series.}}
#'
#' @param img The calculated number or brightness time series image series.
#' @inheritParams number_img
#' @param frames_per_set The number of frames used in the calculation of each
#'   point in the number or brightness time series.
#'
#' @return An object of class `number`.
#'
#' @seealso [number_time_series()], [brightness_time_series()].
#'
#' @name nb-ts-img-classes
NULL

#' @rdname nb-ts-img-classes
#' @export
number_ts_img <- function(img, def, frames_per_set,
                          thresh = NA, tau = NA, filt = NA) {
  checkmate::assert_array(img, min.d = 3, max.d = 4)
  n_ch <- dim(img) %>% {dplyr::if_else(length(.) == 4, .[3], 1L)}
  img %<>% number_img_common(n_ch = n_ch, def = def,
                         thresh = thresh, tau = tau, filt = filt)
  class(img) %<>% c("number_ts_img", .)
  attr(img, "frames_per_set") <- frames_per_set
  img
}

#' @rdname nb-ts-img-classes
#' @export
brightness_ts_img <- function(img, def, frames_per_set,
                              thresh = NA, tau = NA, filt = NA) {
  checkmate::assert_string(def)
  def %<>% tolower()
  if (def == "b") {
    def <- "B"
  } else if (startsWith("epsilon", def)) {
    def <- "epsilon"
  } else {
    stop("'def' must be one of 'B' or 'epsilon'.")
  }
  if (!isTRUE(def %in% c("B", "epsilon")))
    stop("def must be one of 'B' or 'epsilon'.")
  out <- number_ts_img(img, "n", frames_per_set = frames_per_set,
                       thresh = thresh, tau = tau, filt = filt)
  attr(out, "def") <- def
  class(out)[class(out) == "number_ts_img"] <- "brightness_ts_img"
  out
}