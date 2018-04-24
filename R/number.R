#' Calculate number from image series.
#'
#' Given a time stack of images, `number()` performs a calculation of the number
#' for each pixel.
#'
#' @param img A 4-dimensional array of images indexed by `img[y, x, channel,
#'   frame]` (an object of class [ijtiff::ijtiff_img]). The image to perform the
#'   calculation on. To perform this on a file that has not yet been read in,
#'   set this argument to the path to that file (a string).
#' @param def A character. Which definition of number do you want to use, `"n"`
#'   or `"N"`?
#' @param thresh The threshold or thresholding method (see
#'   [autothresholdr::mean_stack_thresh()]) to use on the image prior to
#'   detrending and number calculations.
#' @param tau The exponential detrending parameter to be passed to
#'   [detrendr::img_detrend_exp()]. This can be a positive number or `"auto"`.
#'   Default is no detrending.
#' @param filt Do you want to smooth (`filt = 'mean'`) or median (`filt =
#'   'median'`) filter the number image using [smooth_filter()] or
#'   [median_filter()] respectively? If selected, these are invoked here with a
#'   filter radius of 1 (with corners included, so each median is the median of
#'   9 elements) and with the option `na_count = TRUE`. If you want to
#'   smooth/median filter the number image in a different way, first calculate
#'   the numbers without filtering (`filt = NULL`) using this function and then
#'   perform your desired filtering routine on the result.
#' @param correct Apply the number/brightness correction detailed in equation 7
#'   of Hur et al. (2014). This is another correction for the effects of
#'   bleaching and is needed in addition to the more conventional correction
#'   controlled by the `tau` parameter.
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
#'   Hur K-H, Macdonald PJ, Berk S, Angert CI, Chen Y, Mueller JD (2014)
#'   Quantitative Measurement of Brightness from Living Cells in the Presence of
#'   Photodepletion. PLoS ONE 9(5): e97440. \doi{10.1371/journal.pone.0097440}.
#'
#' @examples
#' img <- ijtiff::read_tif(system.file('extdata', '50.tif', package = 'nandb'))
#' ijtiff::display(img[, , 1, 1])
#' num <- number(img, "N", tau = NA, thresh = "Huang")
#' num <- number(img, "n", tau = 10, thresh = "tri")
#' @export
number <- function(img, def, tau = NULL, thresh = NULL, filt = NULL,
                   correct = FALSE, s = 1, offset = 0, readout_noise = 0,
                   gamma = 1, parallel = FALSE) {
  checkmate::assert_string(def)
  if (! def %in% c("n", "N")) stop("'def' must be one of 'n' or 'N'.")
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
      img %<>% autothresholdr::mean_stack_thresh(method = thresh,
                                                 ignore_na = TRUE)
      thresh_atts <- attr(img, "thresh")
      img <- img[, , 1, ]
    }
    if (!is.na(tau)) {
      img %<>% detrendr::img_detrend_exp(tau, parallel = parallel,
                                         purpose = "FFS")
      tau_atts <- attr(img, "parameter")
      attr(tau_atts, "auto") <- attr(img, "auto")
      img <- img[, , 1, ]
    }
    epsilon <- brightness(img, def = "e", s = s, offset = offset,
                          readout_noise = readout_noise,
                          correct = correct)
    pillar_means <- detrendr::mean_pillars(img, parallel = parallel)
    out <- (pillar_means - offset) / (s * epsilon)
    if (def == "N") {
      out %<>% {(epsilon * .) / (1 + epsilon)}
    } else {
      out %<>% {. / gamma}
    }
    if (length(dim(out)) > 2) dim(out) %<>% {.[1:2]}
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
      for (i in seq_len(n_ch)) {
        out_i <- number(img[, , i, ], def = def, tau = tau[[i]],
                        thresh = thresh[[i]], filt = filt[[i]],
                        correct = correct,
                        s = s, offset = offset, readout_noise = readout_noise,
                        gamma = gamma, parallel = parallel)
        out[, , i] <- out_i
        thresh_atts[[i]] <- attr(out_i, "thresh")
        tau_atts[[i]] <- attr(out_i, "tau")
      }
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
#' img <- ijtiff::read_tif(system.file('extdata', '50.tif', package = 'nandb'))
#' nts <- number_timeseries(img, "n", frames_per_set = 20,
#'                           tau = NA, thresh = "Huang", parallel = 2)
#' @export
number_timeseries <- function(img, def, frames_per_set,
                               tau = NULL, thresh = NULL,
                               filt = NULL, correct = FALSE, s = 1, offset = 0,
                               readout_noise = 0, gamma = 1,
                               parallel = FALSE) {
  if (! def %in% c("n", "N")) stop("'def' must be one of 'n' or 'N'.")
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
      img %<>% autothresholdr::mean_stack_thresh(method = thresh,
                                                 ignore_na = TRUE)
      thresh_atts <- attr(img, "thresh")
      img <- img[, , 1, ]
    }
    if (!is.na(tau)) {
      img %<>% detrendr::img_detrend_exp(tau = tau, purpose = "FFS",
                                         parallel = parallel)
      tau_atts <- attr(img, "parameter")
      attr(tau_atts, "auto") <- attr(img, "auto")
      img <- img[, , 1, ]
    }
    out <- img[, , seq_len(sets)]
    if (sets == 1) dim(out) %<>% c(1)
    for (i in seq_len(sets)) {
      indices_i <- seq((i - 1) * frames_per_set + 1, i * frames_per_set)
      out[, , i] <- number(img[, , indices_i], def = def, filt = filt,
                           s = s, offset = offset,
                           readout_noise = readout_noise, gamma = gamma,
                           parallel = parallel)
    }
  } else {
    frames <- d[4]
    sets <- frames %/% frames_per_set
    out <- img[, , , seq_len(sets)]
    if (sets == 1) dim(out) %<>% c(1)
    thresh_atts <- list()
    tau_atts <- list()
    for (i in seq_len(n_ch)) {
      out_i <- number_timeseries(img[, , i, ], def = def,
                                  frames_per_set = frames_per_set,
                                  tau = tau[[i]], thresh = thresh[[i]],
                                  filt = filt[[i]], offset = offset,
                                  readout_noise = readout_noise,
                                  parallel = parallel)
      out[, , i, ] <- out_i
      thresh_atts[[i]] <- attr(out_i, "thresh")
      tau_atts[[i]] <- attr(out_i, "tau")
    }
  }
  number_ts_img(out, def = def, frames_per_set = frames_per_set,
                thresh = thresh_atts, tau = tau_atts, filt = filt)
}

number_file <- function(path, def, tau = NULL,
                        thresh = NULL, filt = NULL,
                        s = 1, offset = 0, readout_noise = 0, gamma = 1,
                        parallel = FALSE) {
  if (! def %in% c("n", "N")) stop("'def' must be one of 'n' or 'N'.")
  checkmate::assert_file_exists(path)
  need_to_change_dir <- stringr::str_detect(path, "/")
  if (need_to_change_dir) {
    dir <- filesstrings::str_before_last(path, "/")
    cwd <- getwd()
    on.exit(setwd(cwd))
    setwd(dir)
    path %<>% filesstrings::str_after_last("/")
  }
  num <- number(path, def, tau = tau,
                thresh = thresh, filt = filt,
                s = s, offset = offset, readout_noise = readout_noise,
                gamma = gamma, parallel = parallel)
  num[abs(num) >= float_max()] <- NA
  suppressMessages(filesstrings::create_dir("number"))
  path %<>% filesstrings::before_last_dot() %>%
    paste0("number", "/", ., make_nb_filename_ending(num)) %>%
    deduplicate_nb_filename()
  ijtiff::write_tif(num, path)
}

number_timeseries_file <- function(path, def, frames_per_set,
                                    tau = NULL, thresh = NULL,
                                    filt = NULL, s = 1, offset = 0,
                                    readout_noise = 0, gamma = 1,
                                    parallel = FALSE) {
  if (! def %in% c("n", "N")) stop("'def' must be one of 'n' or 'N'.")
  checkmate::assert_file_exists(path)
  need_to_change_dir <- stringr::str_detect(path, "/")
  if (need_to_change_dir) {
    dir <- filesstrings::str_before_last(path, "/")
    cwd <- getwd()
    on.exit(setwd(cwd))
    setwd(dir)
    path %<>% filesstrings::str_after_last("/")
  }
  nts <- number_timeseries(path, def, frames_per_set = frames_per_set,
                            tau = tau,
                            thresh = thresh, filt = filt,
                            s = s, offset = offset,
                            readout_noise = readout_noise, gamma = gamma,
                            parallel = parallel)
  nts[abs(nts) >= float_max()] <- NA
  suppressMessages(filesstrings::create_dir("number_timeseries"))
  path %<>% filesstrings::before_last_dot() %>%
    paste0("number_timeseries", "/", ., make_nb_filename_ending(nts)) %>%
    deduplicate_nb_filename()
  ijtiff::write_tif(nts, path)
}


#' Number calculations for every image in a folder.
#'
#' Perform [number()] calculations on all tif images in a folder and save the
#' resulting number images to disk.
#'
#' @note Extreme number values (of magnitude greater than 3.40282e+38) will be
#'   written to the TIFF file as `NA`, since TIFF files cannot handle such huge
#'   numbers.
#'
#' @param folder_path The path (relative or absolute) to the folder you wish to
#'   process.
#'
#' @inheritParams number
#'
#' @seealso [number()]
#'
#' @examples
#' \dontrun{
#' setwd(tempdir())
#' img <- ijtiff::read_tif(system.file('extdata', '50.tif', package = 'nandb'))
#' ijtiff::write_tif(img, 'img2.tif')
#' number_folder(def = "n", tau = NA, thresh = "Huang", parallel = 2)
#' }
#' @export
number_folder <- function(folder_path = ".", def,
                          tau = NULL, thresh = NULL, filt = NULL,
                          s = 1, offset = 0, readout_noise = 0, gamma = 1,
                          parallel = FALSE) {
  if (! def %in% c("n", "N")) stop("'def' must be one of 'n' or 'N'.")
  init_dir <- getwd()
  on.exit(setwd(init_dir))
  setwd(folder_path)
  file_names <- list.files(pattern = "\\.tif")
  purrr::map(file_names, number_file, def = def, tau = tau,
             thresh = thresh, filt = filt, s = s, offset = offset,
             readout_noise = readout_noise, gamma = gamma,
             parallel = parallel) %>%
    magrittr::set_names(filesstrings::before_last_dot(file_names)) %>%
    invisible()
}

#' Number time-series calculations for every image in a folder.
#'
#' Perform [number_timeseries()] calculations on all tif images in a folder and
#' save the resulting number images to disk.
#'
#' @note Extreme number values (of magnitude greater than 3.40282e+38) will be
#'   written to the TIFF file as `NA`, since TIFF files cannot handle such huge
#'   numbers.
#'
#' @inheritParams number
#' @inheritParams number_timeseries
#' @inheritParams number_folder
#'
#' @seealso [number_timeseries()]
#'
#' @examples
#' \dontrun{
#' setwd(tempdir())
#' img <- ijtiff::read_tif(system.file('extdata', '50.tif', package = 'nandb'))
#' ijtiff::write_tif(img, 'img1.tif')
#' ijtiff::write_tif(img, 'img2.tif')
#' number_timeseries_folder(def = "n", tau = NA, thresh = "Huang",
#'                           frames_per_set = 20, parallel = 2)
#' }
#' @export
number_timeseries_folder <- function(folder_path = ".", def, frames_per_set,
                                      tau = NULL, thresh = NULL,
                                      filt = NULL, s = 1, offset = 0,
                                      readout_noise = 0, gamma = 1,
                                      parallel = FALSE) {
  if (! def %in% c("n", "N")) stop("'def' must be one of 'n' or 'N'.")
  init_dir <- getwd()
  on.exit(setwd(init_dir))
  setwd(folder_path)
  file_names <- list.files(pattern = "\\.tif")
  purrr::map(file_names, number_timeseries_file, def = def,
             frames_per_set = frames_per_set,
             tau = tau, thresh = thresh, filt = filt,
             s = s, offset = offset, readout_noise = readout_noise,
             gamma = gamma, parallel = parallel) %>%
    magrittr::set_names(filesstrings::before_last_dot(file_names)) %>%
    invisible()
}

