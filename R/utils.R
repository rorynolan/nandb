extend_for_all_chs <- function(x, n_ch) {
  if (is.null(x)) x <- list(NA)[rep(1, n_ch)]
  if (length(x) == 1) {
    if (is.list(x)) {
      x <- x[rep(1, n_ch)]
    } else {
      x <- list(x)[rep(1, n_ch)]
    }
  }
  sq <- seq_len(n_ch)
  x <- x[sq]
  for (i in sq) if (is.null(x[[i]])) x[[i]] <- NA
  x
}

make_nb_filename_ending <- function(nb_img) {
  checkmate::assert_array(nb_img, d = 4)
  def <- attr(nb_img, "def")
  tau <- attr(nb_img, "tau")
  auto <- attr(tau, "auto")
  thresh <- attr(nb_img, "thresh")
  pasted_class <- paste(class(nb_img), collapse = "")
  nb <- dplyr::if_else(stringr::str_detect(pasted_class, "number"),
                       "number", "brightness")
  is_ts <- stringr::str_detect(pasted_class, "_ts_")
  d <- dim(nb_img)
  n_ch <- d[3]
  if ("autothresh_method" %in% names(attributes(thresh))) {
    thresh_method <- attr(thresh, "autothresh_method")
  } else {
    thresh_method <- rep(NA, n_ch)
  }
  filt <- attr(nb_img, "filt")
  stopifnot(filesstrings::all_equal(c(length(tau), length(auto), length(thresh),
                                      length(thresh_method), length(filt))))
  tau_part <- ""
  for (i in seq_len(n_ch)) {
    tau_part %<>% paste0(dplyr::if_else(auto[i], "auto=", ""), tau[i], ",")
  }
  tau_part %<>% stringr::str_sub(1, -2)
  thresh_part <- ""
  for (i in seq_len(n_ch)) {
    thresh_part %<>% paste0(dplyr::if_else(is.na(thresh_method[i]), "",
                                           paste0(thresh_method[i], "=")),
                            thresh[i], ",")
  }
  thresh_part %<>% stringr::str_sub(1, -2)
  paste0("_", nb, "_", def, "_",
         dplyr::if_else(is_ts, paste0("timeseries", "_", "frames=",
                                      attr(nb_img, "frames_per_set"), "_"), ""),
         "tau=", tau_part, "_",
         "thresh=", thresh_part, "_",
         "filt=", paste(filt, collapse = ","))
}

fix_filt <- function(filt) {
  if (is.null(filt)) filt <- NA_character_
  checkmate::assert_character(filt)
  filt %<>% tolower()
  filt[startsWith("mean", filt)] <- "mean"
  filt[startsWith("median", filt)] <- "median"
  if (!all((filt %in% c("mean", "median")) | is.na(filt)))
    stop("All elements of 'filt' must be either 'mean', 'median' or 'NA'.")
  filt
}

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

deduplicate_nb_filename <- function(path) {
  checkmate::assert_string(path)
  has_been_through_detrendr <- path %>%
    stringr::str_detect(paste0("detrended_",
                               c("boxcar", "exponential", "polynomial"))) %>%
    any()
  if (any(stringr::str_detect(path, paste0("detrended_",
                                           c("boxcar",
                                             "exponential",
                                             "polynomial"))))) {
    thresh_pattern <- "_thresh=.*?_"
    n_threshs <- filesstrings::count_matches(path, thresh_pattern)
    if (n_threshs == 2) {
      second_thresh_indices <- stringr::str_locate_all(path,
                                                       thresh_pattern)[[1]][2, ]
      if (stringr::str_detect(stringr::str_sub(path,
                                               second_thresh_indices["start"],
                                               second_thresh_indices["end"]),
                              "_thresh=NA,*(NA,)*(NA)*_")) {
        if (second_thresh_indices["end"] == nchar(path)) {
          path %<>% stringr::str_sub(1, second_thresh_indices["start"] - 1)
        } else {
          path <- paste0(
                    stringr::str_sub(path,
                                     1, second_thresh_indices["start"] - 1),
                    stringr::str_sub(path,
                                     second_thresh_indices["end"], -1))
        }
      }
    }
    tau_pattern <- "_tau=.*?_"
    n_taus <- filesstrings::count_matches(path, tau_pattern)
    if (n_taus == 2) {
      second_tau_indices <- stringr::str_locate_all(path,
                                                    tau_pattern)[[1]][2, ]
      if (stringr::str_detect(stringr::str_sub(path,
                                               second_tau_indices["start"],
                                               second_tau_indices["end"]),
                              "_tau=NA,*(NA,)*(NA)*_")) {
        if (second_tau_indices["end"] == nchar(path)) {
          path %<>% stringr::str_sub(1, second_tau_indices["start"] - 1)
        } else {
          path <- paste0(stringr::str_sub(path,
                                          1, second_tau_indices["start"] - 1),
                         stringr::str_sub(path,
                                          second_tau_indices["end"], -1))
        }
      }
    }
  }
  path
}
