#' Calculate brightness from image series.
#'
#' Given a time stack of images, `brightness()` performs a calculation of the
#' brightness for each pixel.
#'
#' @param def A character. Which definition of brightness do you want to use,
#'   `"B"` or `"epsilon"`?
#' @inheritParams detrendr::img_detrend_exp
#' @inheritParams number
#'
#' @return A matrix, the brightness image.
#'
#' @references Digman MA, Dalal R, Horwitz AF, Gratton E. Mapping the Number of
#'   Molecules and Brightness in the Laser Scanning Microscope. Biophysical
#'   Journal. 2008;94(6):2320-2332. \doi{10.1529/biophysj.107.114645}.
#'
#'   Dalal, RB, Digman, MA, Horwitz, AF, Vetri, V, Gratton, E (2008).
#'   Determination of particle number and brightness using a laser scanning
#'   confocal microscope operating in the analog mode. Microsc. Res. Tech., 71,
#'   1:69-81. \doi{10.1002/jemt.20526}.
#'
#' @examples
#' img <- ijtiff::read_tif(system.file('extdata', '50.tif', package = 'nandb'))
#' ijtiff::display(img[, , , 1])
#' b <- brightness(img, "e", tau = NA, thresh = "Huang")
#' b <- brightness(img, "B", tau = 10, thresh = "tri")
#' @export
brightness <- function(img, def, tau = NULL,
                       thresh = NULL, fail = NA, filt = NULL,
                       s = 1, offset = 0, readout_noise = 0, gamma = 1,
                       seed = NULL, parallel = FALSE) {
  checkmate::assert_string(def)
  if (startsWith("epsilon", tolower(def))) def <- "epsilon"
  if (def == "b") def <- "B"
  if (! def %in% c("epsilon", "B"))
    stop("'def' must be one of 'B' or 'epsilon'.")
  img %<>% nb_get_img()
  d <- dim(img)
  n_ch <- dplyr::if_else(length(d) == 3, 1L, d[3])
  if (n_ch == 1 && length(d) == 4) img %<>% {.[, , 1, ]}
  thresh %<>% extend_for_all_chs(n_ch)
  tau %<>% extend_for_all_chs(n_ch)
  if (!is.null(filt)) filt %<>% fix_filt()
  filt %<>% extend_for_all_chs(n_ch) %>% unlist() %>% as.character()
  tau_atts <- extend_for_all_chs(rlang::set_attrs(NA, auto = FALSE),
                                 n_ch)
  thresh_atts <- extend_for_all_chs(NA, n_ch)
  if (n_ch == 1) {
    if (!is.na(thresh)) {
      img %<>% autothresholdr::mean_stack_thresh(method = thresh, fail = fail)
      thresh_atts <- attr(img, "thresh")
      img <- img[, , 1, ]
    }
    if (!is.na(tau)) {
      img %<>% detrendr::img_detrend_exp(tau, seed = seed, parallel = parallel)
      tau_atts <- attr(img, "parameter")
      attr(tau_atts, "auto") <- attr(img, "auto")
      img <- img[, , 1, ]
    }
    out <- (detrendr::var_pillars(img, parallel = parallel) - readout_noise) /
      (detrendr::mean_pillars(img, parallel = parallel) - offset)
    if (def == "epsilon") out <- out - s
    if (!is.na(filt)) {
      checkmate::assert_string(filt)
      if (filt == "median") {
        out %<>% median_filter(na_count = TRUE)
      } else {
        out %<>% smooth_filter(na_count = TRUE)
      }
    }
  } else {
    out <- img[, , , 1]
    thresh_atts <- list()
    tau_atts <- list()
    for (i in seq_len(n_ch)) {
      if (!is.null(seed)) seed <- seed + i
      out_i <- brightness(img[, , i, ], def = def, tau = tau[[i]],
                          thresh = thresh[[i]], filt = filt[[i]],
                          s = s, offset = offset, readout_noise = readout_noise,
                          gamma = gamma, seed = seed, parallel = parallel)
      out[, , i] <- out_i
      thresh_atts[[i]] <- attr(out_i, "thresh")
      tau_atts[[i]] <- attr(out_i, "tau")
    }
  }
  brightness_img(out, def, thresh_atts, tau_atts, filt)
}

#' Create a brightness time-series.
#'
#' Given a stack of images `img`, use the first `frames_per_set` of them to
#' create one brightness image, the next `frames_per_set` of them to create the
#' next brightness image and so on to get a time-series of brightness images.
#'
#' @param frames_per_set The number of frames with which to calculate the
#'   successive brightnesses.
#'
#' This may discard some images, for example if 175 frames are in the input and
#' `frames_per_set = 50`, then the last 25 are discarded. If bleaching
#' correction is selected, it is performed on the whole image stack before the
#' sectioning is done for calculation of numbers.
#'
#' @inheritParams brightness
#' @inheritParams number
#'
#' @return An object of class [brightness_ts_img].
#'
#'   \itemize{\item If `img` is 3-dimensional (i.e. 1-channel), a 3-dimensional
#'   array `arr` is returned with `arr[y, x, t]` being pixel \eqn{(x, y)} of the
#'   \eqn{t}th brightness image in the brightness time series. \item If  `img`
#'   is 4-dimensional (i.e. 2-channel), a 4-dimensional array `arr` is returned
#'   with `arr[y, x, c, t]` being pixel \eqn{(x, y)} of the \eqn{c}th channel of
#'   the \eqn{t}th brightness image in the brightness time series.}
#'
#' @seealso [brightness()].
#'
#' @examples
#' img <- ijtiff::read_tif(system.file('extdata', '50.tif', package = 'nandb'))
#' bts <- brightness_time_series(img, "e", frames_per_set = 20,
#'                               tau = NA, thresh = "Huang", parallel = 2)
#' @export
brightness_time_series <- function(img, def, frames_per_set,
                                   tau = NULL, thresh = NULL, fail = NA,
                                   filt = NULL, s = 1, offset = 0,
                                   readout_noise = 0, gamma = 1,
                                   parallel = FALSE, seed = NULL) {
  if (startsWith("epsilon", tolower(def))) def <- "epsilon"
  if (def == "b") def <- "B"
  img %<>% nb_get_img()
  d <- dim(img)
  n_ch <- dplyr::if_else(length(d) == 3, 1L, d[3])
  if (n_ch == 1 && length(d) == 4) img %<>% {.[, , 1, ]}
  thresh %<>% extend_for_all_chs(n_ch)
  tau %<>% extend_for_all_chs(n_ch)
  if (!is.null(filt)) filt %<>% fix_filt()
  filt %<>% extend_for_all_chs(n_ch) %>% unlist() %>% as.character()
  tau_atts <- extend_for_all_chs(rlang::set_attrs(NA, auto = FALSE),
                                 n_ch)
  thresh_atts <- extend_for_all_chs(NA, n_ch)
  if (n_ch == 1) {
    frames <- dim(img)[3]
    if (frames < frames_per_set) {
      stop("You have selected ", frames_per_set, " frames per set, ",
           "but there are only ", frames, " frames in total.")
    }
    sets <- frames %/% frames_per_set
    if (!is.na(thresh)) {
      img %<>% autothresholdr::mean_stack_thresh(method = thresh, fail = fail)
      thresh_atts <- attr(img, "thresh")
      img <- img[, , 1, ]
    }
    if (!is.na(tau)) {
      img %<>% detrendr::img_detrend_exp(tau = tau,
                                         seed = seed, parallel = parallel)
      tau_atts <- attr(img, "parameter")
      attr(tau_atts, "auto") <- attr(img, "auto")
      img <- img[, , 1, ]
    }
    out <- img[, , seq_len(sets)]
    if (sets == 1) dim(out) %<>% c(1)
    for (i in seq_len(sets)) {
      if (!is.null(seed)) seed <- seed + i
      indices_i <- seq((i - 1) * frames_per_set + 1, i * frames_per_set)
      out[, , i] <- brightness(img[, , indices_i], def = def, filt = filt,
                               s = s, offset = offset,
                               readout_noise = readout_noise, gamma = gamma,
                               seed = seed, parallel = parallel)
    }
  } else {
    frames <- d[4]
    sets <- frames %/% frames_per_set
    out <- img[, , , seq_len(sets)]
    if (sets == 1) dim(out) %<>% c(1)
    thresh_atts <- list()
    tau_atts <- list()
    for (i in seq_len(n_ch)) {
      if (!is.null(seed)) seed <- seed + i
      out_i <- brightness_time_series(img[, , i, ], def = def,
                                     frames_per_set = frames_per_set,
                                      tau = tau[[i]], thresh = thresh[[i]],
                                      filt = filt[[i]], offset = offset,
                                      readout_noise = readout_noise,
                                      seed = seed + i, parallel = parallel)
      out[, , i, ] <- out_i
      thresh_atts[[i]] <- attr(out_i, "thresh")
      tau_atts[[i]] <- attr(out_i, "tau")
    }
  }
  brightness_ts_img(out, def = def, frames_per_set = frames_per_set,
                    thresh = thresh_atts, tau = tau_atts, filt = filt)
}

brightness_file <- function(path, def, tau = NULL,
                            thresh = NULL, fail = NA, filt = NULL,
                            s = 1, offset = 0, readout_noise = 0, gamma = 1,
                            seed = NULL, parallel = FALSE) {
  checkmate::assert_string(path)
  b <- brightness(path, def, tau = tau,
                  thresh = thresh, fail = fail, filt = filt,
                  s = s, offset = offset, readout_noise = readout_noise,
                  gamma = gamma, seed = seed, parallel = parallel)
  path %<>% filesstrings::before_last_dot()
  ijtiff::write_tif(b, paste0(path, make_nb_filename_ending(b)))
}

brightness_time_series_file <- function(path, def, frames_per_set,
                                        tau = NULL, thresh = NULL, fail = NA,
                                        filt = NULL, s = 1, offset = 0,
                                        readout_noise = 0, gamma = 1,
                                        parallel = FALSE, seed = NULL) {
  if (startsWith("epsilon", tolower(def))) def <- "epsilon"
  if (def == "b") def <- "B"
  checkmate::assert_string(path)
  bts <- brightness_time_series(path, def, frames_per_set = frames_per_set,
                                tau = tau,
                                thresh = thresh, fail = fail, filt = filt,
                                s = s, offset = offset,
                                readout_noise = readout_noise, gamma = gamma,
                                seed = seed, parallel = parallel)
  path %<>% filesstrings::before_last_dot()
  ijtiff::write_tif(bts, paste0(path, make_nb_filename_ending(bts)))
}


#' Brightness calculations for every image in a folder.
#'
#' Perform [brightness()] calculations on all tif images in a folder and save the
#' resulting brightness images to disk.
#'
#' @inheritParams brightness
#' @inheritParams number
#' @inheritParams number_folder
#'
#' @seealso [number()]
#'
#' @examples
#' setwd(tempdir())
#' img <- ijtiff::read_tif(system.file('extdata', '50.tif', package = 'nandb'))
#' ijtiff::write_tif(img, 'img1.tif')
#' ijtiff::write_tif(img, 'img2.tif')
#' brightness_folder(def = "B", tau = NA, thresh = "Huang", parallel = 2)
#' suppressWarnings(file.remove(list.files()))  # cleanup
#' @export
brightness_folder <- function(folder_path = ".", def,
                              tau = NULL, thresh = NULL, fail = NA, filt = NULL,
                              s = 1, offset = 0, readout_noise = 0, gamma = 1,
                              seed = NULL, parallel = FALSE) {
  if (startsWith("epsilon", tolower(def))) def <- "epsilon"
  if (def == "b") def <- "B"
  init_dir <- getwd()
  on.exit(setwd(init_dir))
  setwd(folder_path)
  file_names <- list.files(pattern = "\\.tiff*$")
  purrr::map(file_names, brightness_file, def = def, tau = tau,
             thresh = thresh, fail = fail, filt = filt, s = s, offset = offset,
             readout_noise = readout_noise, gamma = gamma,
             seed = seed, parallel = parallel) %>%
    magrittr::set_names(filesstrings::before_last_dot(file_names)) %>%
    invisible()
}

#' Brightness time-series calculations for every image in a folder.
#'
#' Perform [brightness_time_series()] calculations on all tif images in a folder
#' and save the resulting number images to disk.
#'
#' @inheritParams brightness
#' @inheritParams brightness_time_series
#' @inheritParams brightness_folder
#' @inheritParams number
#' @inheritParams number_time_series
#' @inheritParams number_folder
#'
#' @seealso [brightness_time_series()]
#'
#' @examples
#' setwd(tempdir())
#' img <- ijtiff::read_tif(system.file('extdata', '50.tif', package = 'nandb'))
#' ijtiff::write_tif(img, 'img1.tif')
#' ijtiff::write_tif(img, 'img2.tif')
#' brightness_time_series_folder(def = "e", tau = NA, thresh = "Huang",
#'                               frames_per_set = 20, parallel = 2)
#' suppressWarnings(file.remove(list.files()))  # cleanup
#' @export
brightness_time_series_folder <- function(folder_path = ".", def,
                                          frames_per_set, tau = NULL,
                                          thresh = NULL, fail = NA, filt = NULL,
                                          s = 1, offset = 0, readout_noise = 0,
                                          gamma = 1, seed = NULL,
                                          parallel = FALSE) {
  if (startsWith("epsilon", tolower(def))) def <- "epsilon"
  if (def == "b") def <- "B"
  init_dir <- getwd()
  on.exit(setwd(init_dir))
  setwd(folder_path)
  file_names <- list.files(pattern = "\\.tif")
  purrr::map(file_names, brightness_time_series_file, def = def,
             frames_per_set = frames_per_set,
             tau = tau, thresh = thresh, fail = fail, filt = filt,
             s = s, offset = offset, readout_noise = readout_noise,
             gamma = gamma, seed = seed, parallel = parallel) %>%
    magrittr::set_names(filesstrings::before_last_dot(file_names)) %>%
    invisible()
}

