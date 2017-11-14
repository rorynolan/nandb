#' Calculate number from image series.
#'
#' Given a time stack of images, `number()` performs a calculation of the number
#' for each pixel. `number_txt_folder()` does this calculation for an entire
#' folder, writing the results as text files via [write_txt_img()].
#'
#' @param img The image to perform the calculation on. This can be a
#'   single-channel image (which would be a 3-d array `img[y, x, frame]`) or a
#'   multi channel (`img[y, x, channel, frame]`). To perform this on a file that
#'   has not yet been read in, set this argument to the path to that file (a
#'   string).
#' @param def A character. Which definition of number do you want to use, `"n"`
#'   or `"N"`?
#' @param n_ch The number of channels in the image (default 1).
#' @param tau The exponential parameter to be passed to
#'   [detrendr::img_detrend_exp()]. This can be a positive number or `"auto"`.
#' @param thresh The threshold or thresholding method (see
#'   [autothresholdr::mean_stack_thresh()]) to use on the image prior to
#'   detrending and brightness calculations.
#' @param fail If thresholding is done, to which value should pixels not
#'   exceeding the threshold be set?
#' @param filt Do you want to smooth (`filt = 'mean'`) or median (`filt =
#'   'median'`) filter the number image using [smooth_filter()] or
#'   [median_filter()] respectively? If selected, these are invoked here with a
#'   filter radius of 1 (with corners included, so each median is the median of
#'   9 elements) and with the option `na_count = TRUE`. If you want to
#'   smooth/median filter the number image in a different way, first calculate
#'   the numbers without filtering (`filt = NULL`) using this function and then
#'   perform your desired filtering routine on the result.
#' @param offset,readout_noise Microscope acquisition parameters. See reference
#'   Dalal et al.
#' @param s A number. The \eqn{S}-factor of microscope acquisition.
#' @param gamma Factor for correction of number \eqn{n} due to the illumination
#'   profile. The default (`gamma = 1`) has no effect. Changing gamma will have
#'   the effect of dividing the result by `gamma`, so the result with `gamma =
#'   0.5` is two times the result with `gamma = 1`. For a Gaussian illumination
#'   profile, use `gamma = 0.3536`; for a Gaussian-Lorentzian illumination
#'   profile, use `gamma = 0.0760`.
#' @inheritParams detrendr::img_detrend_exp
#'
#' @return A matrix, the number image.
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
#' img <- read_tif(system.file('extdata', '50.tif', package = 'nandb'))
#' display(img[, , 1])
#' num <- number(img, "N", tau = NA, thresh = "Huang")
#' num <- number(img, "n", tau = 10, thresh = "tri")
#' @export
number <- function(img, def, n_ch = 1, tau = NULL,
                   thresh = NULL, fail = NA, filt = NULL,
                   s = 1, offset = 0, readout_noise = 0, gamma = 1,
                   seed = NULL, parallel = FALSE) {
  checkmate::assert_string(def)
  if (! def %in% c("n", "N")) stop("'def' must be one of 'n' or 'N'.")
  img %<>% nb_get_img(n_ch = n_ch, min_d = 3, max_d = 4)
  thresh %<>% extend_for_all_chs(n_ch)
  tau %<>% extend_for_all_chs(n_ch)
  if (!is.null(filt)) filt %<>% fix_filt()
  filt %<>% extend_for_all_chs(n_ch)
  tau_atts <- radiant.data::set_attr(NA, "auto", FALSE)
  thresh_atts <- NA
  if (n_ch == 1) {
    if (!is.na(thresh)) {
      img %<>% autothresholdr::mean_stack_thresh(method = thresh, fail = fail)
      thresh_atts <- attr(img, "thresh")
    }
    if (!is.na(tau)) {
      img %<>% detrendr::img_detrend_exp(tau, seed = seed, parallel = parallel)
      tau_atts <- attr(img, "parameter")
      attr(tau_atts, "auto") <- attr(img, "auto")
    }
    if (def == "N") {
      out <- (detrendr::mean_pillars(img, parallel = parallel) - offset) ^ 2 /
        (detrendr::var_pillars(img, parallel = parallel) - readout_noise)
    } else {
      out <- (detrendr::mean_pillars(img, parallel = parallel) - offset) %>% {
        . ^ 2 / (detrendr::var_pillars(img, parallel = parallel) -
                   readout_noise - s * .)
      } / gamma
    }
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
    if (!is.null(seed)) seed <- seed + i
    for (i in seq_len(n_ch)) {
      out_i <- number(img[, , i, ], def = def, tau = tau[[i]],
                      thresh = thresh[[i]], filt = filt[[i]],
                      s = s, offset = offset, readout_noise = readout_noise,
                      gamma = gamma, seed = seed, parallel = parallel)
      out[, , i] <- out_i
      thresh_atts[[i]] <- attr(out_i, "thresh")
      tau_atts[[i]] <- attr(out_i, "tau")
    }
  }
  number_img(out, def, thresh_atts, tau_atts, filt)
}


#' Create a number time-series.
#'
#' Given a stack of images `img`, use the first `frames_per_set` of them to
#' create one number image, the next `frames_per_set` of them to create the next
#' number image and so on to get a time-series of number images.
#'
#' This may discard some images, for example if 175 frames are in the input and
#' `frames_per_set = 50`, then the last 25 are discarded. If bleaching
#' correction is selected, it is performed on the whole image stack before the
#' sectioning is done for calculation of numbers.
#'
#' @inheritParams number
#' @param frames_per_set The number of frames with which to calculate the
#'   successive numbers.
#'
#' @return An object of class [number_ts_img].
#'
#'   \itemize{\item If `img` is 3-dimensional (i.e. 1-channel), a 3-dimensional
#'   array `arr` is returned with `arr[y, x, t]` being pixel \eqn{(x, y)} of the
#'   \eqn{t}th number image in the number time series. \item If  `img` is
#'   4-dimensional (i.e. 2-channel), a 4-dimensional array `arr` is returned
#'   with `arr[y, x, c, t]` being pixel \eqn{(x, y)} of the \eqn{c}th channel of
#'   the \eqn{t}th number image in the number time series.}
#'
#' @seealso [number()].
#'
#' @examples
#' img <- read_tif(system.file('extdata', '50.tif', package = 'nandb'))
#' nts <- number_time_series(img, "n", frames_per_set = 20,
#'                           tau = NA, thresh = "Huang", parallel = 2)
#' @export
number_time_series <- function(img, def, frames_per_set, n_ch = 1,
                               tau = NULL, thresh = NULL, fail = NA,
                               filt = NULL, s = 1, offset = 0,
                               readout_noise = 0, gamma = 1,
                               parallel = FALSE, seed = NULL) {
  if (! def %in% c("n", "N")) stop("'def' must be one of 'n' or 'N'.")
  img %<>% nb_get_img(n_ch = n_ch, min_d = 3, max_d = 4)
  d <- dim(img)
  thresh %<>% extend_for_all_chs(n_ch)
  tau %<>% extend_for_all_chs(n_ch)
  if (!is.null(filt)) filt %<>% fix_filt()
  filt %<>% extend_for_all_chs(n_ch)
  thresh_atts <- NA
  tau_atts <- radiant.data::set_attr(NA, "auto", FALSE)
  if (n_ch == 1) {
    frames <- d[3]
    if (frames < frames_per_set) {
      stop("You have selected ", frames_per_set, " frames per set, ",
           "but there are only ", frames, " frames in total.")
    }
    sets <- frames %/% frames_per_set
    if (!is.na(thresh)) {
      img %<>% autothresholdr::mean_stack_thresh(method = thresh, fail = fail)
      thresh_atts <- attr(img, "thresh")
    }
    if (!is.na(tau)) {
      img %<>% detrendr::img_detrend_exp(tau = tau,
                                         seed = seed, parallel = parallel)
      tau_atts <- attr(img, "parameter")
      attr(tau_atts, "auto") <- attr(img, "auto")
    }
    out <- img[, , seq_len(sets)]
    if (sets == 1) dim(out) %<>% c(1)
    for (i in seq_len(sets)) {
      if (!is.null(seed)) seed <- seed + i
      indices_i <- seq((i - 1) * frames_per_set + 1, i * frames_per_set)
      out[, , i] <- number(img[, , indices_i], def = def, filt = filt,
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
      out_i <- number_time_series(img[, , i, ], def = def,
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
  number_ts_img(out, def = def, frames_per_set = frames_per_set,
                thresh = thresh_atts, tau = tau_atts, filt = filt)
}

number_file <- function(path, def, n_ch = 1, tau = NULL,
                        thresh = NULL, fail = NA, filt = NULL,
                        s = 1, offset = 0, readout_noise = 0, gamma = 1,
                        seed = NULL, parallel = FALSE, rds = FALSE) {
  if (! def %in% c("n", "N")) stop("'def' must be one of 'n' or 'N'.")
  checkmate::assert_string(path)
  num <- number(path, def, n_ch = n_ch, tau = tau,
                thresh = thresh, fail = fail, filt = filt,
                s = s, offset = offset, readout_noise = readout_noise,
                gamma = gamma, seed = seed, parallel = parallel)
  path %<>% filesstrings::before_last_dot()
  write_txt_img(num, paste0(path, make_nb_filename_ending(num)), rds = rds)
}

number_time_series_file <- function(path, def, frames_per_set, n_ch = 1,
                                    tau = NULL, thresh = NULL, fail = NA,
                                    filt = NULL, s = 1, offset = 0,
                                    readout_noise = 0, gamma = 1,
                                    parallel = FALSE, seed = NULL,
                                    rds = FALSE) {
  if (! def %in% c("n", "N")) stop("'def' must be one of 'n' or 'N'.")
  checkmate::assert_string(path)
  nts <- number_time_series(path, def, frames_per_set = frames_per_set,
                            n_ch = n_ch, tau = tau,
                            thresh = thresh, fail = fail, filt = filt,
                            s = s, offset = offset,
                            readout_noise = readout_noise, gamma = gamma,
                            seed = seed, parallel = parallel)
  path %<>% filesstrings::before_last_dot()
  write_txt_img(nts, paste0(path, make_nb_filename_ending(nts)), rds = rds)
}


#' Number calculations for every image in a folder.
#'
#' Perform [number()] calculations on all tif images in a folder and save the
#' resulting number images to disk as text images (and optionally also as RDS
#' files).
#'
#' @param folder_path The path (relative or absolute) to the folder you wish to
#'   process.
#'
#' @inheritParams number
#' @inheritParams write_txt_img
#'
#' @seealso [number()]
#'
#' @examples
#' setwd(tempdir())
#' img <- read_tif(system.file('extdata', '50.tif', package = 'nandb'))
#' write_tif(img, 'img1.tif')
#' write_tif(img, 'img2.tif')
#' number_folder(def = "n", tau = NA, thresh = "Huang", parallel = 2, n_ch = 1)
#' suppressWarnings(file.remove(list.files()))  # cleanup
#' @export
number_folder <- function(folder_path = ".", def, n_ch = 1,
                          tau = NULL, thresh = NULL, fail = NA, filt = NULL,
                          s = 1, offset = 0, readout_noise = 0, gamma = 1,
                          seed = NULL, parallel = FALSE, rds = FALSE) {
  if (! def %in% c("n", "N")) stop("'def' must be one of 'n' or 'N'.")
  init_dir <- getwd()
  on.exit(setwd(init_dir))
  setwd(folder_path)
  file_names <- list.files(pattern = "\\.tif")
  purrr::map(file_names, number_file, def = def, n_ch = n_ch, tau = tau,
             thresh = thresh, fail = fail, filt = filt, s = s, offset = offset,
             readout_noise = readout_noise, gamma = gamma,
             seed = seed, parallel = parallel, rds = rds) %>%
    magrittr::set_names(filesstrings::before_last_dot(file_names)) %>%
    invisible()
}

#' Number time-series calculations for every image in a folder.
#'
#' Perform [number_time_series()] calculations on all tif images in a folder and
#' save the resulting number images to disk as text images (and optionally also
#' as RDS files).
#'
#' @inheritParams number
#' @inheritParams number_time_series
#' @inheritParams number_folder
#' @inheritParams write_txt_img
#'
#' @seealso [number_time_series()]
#'
#' @examples
#' setwd(tempdir())
#' img <- read_tif(system.file('extdata', '50.tif', package = 'nandb'))
#' write_tif(img, 'img1.tif')
#' write_tif(img, 'img2.tif')
#' number_time_series_folder(def = "n", tau = NA, thresh = "Huang",
#'                           frames_per_set = 20, parallel = 2, n_ch = 1)
#' suppressWarnings(file.remove(list.files()))  # cleanup
#' @export
number_time_series_folder <- function(folder_path = ".", def, frames_per_set,
                                      n_ch = 1, tau = NULL,
                                      thresh = NULL, fail = NA, filt = NULL,
                                      s = 1, offset = 0, readout_noise = 0,
                                      gamma = 1, seed = NULL, parallel = FALSE,
                                      rds = FALSE) {
  if (! def %in% c("n", "N")) stop("'def' must be one of 'n' or 'N'.")
  init_dir <- getwd()
  on.exit(setwd(init_dir))
  setwd(folder_path)
  file_names <- list.files(pattern = "\\.tif")
  purrr::map(file_names, number_time_series_file, def = def, n_ch = n_ch,
             frames_per_set = frames_per_set,
             tau = tau, thresh = thresh, fail = fail, filt = filt,
             s = s, offset = offset, readout_noise = readout_noise,
             gamma = gamma, seed = seed, parallel = parallel, rds = rds) %>%
    magrittr::set_names(filesstrings::before_last_dot(file_names)) %>%
    invisible()
}

