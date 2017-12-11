test_that("number works", {
  context("number()")
  img <- system.file('extdata', '50.tif', package = 'nandb')
  number <- number(img, "n", tau = 'auto', thresh = 'Huang',
                           filt = 'median', seed = 0, parallel = 2)
  expect_equal(round(median(number, na.rm = TRUE), 3), -31.781)
  number <- number(img, "N", tau = 1000, thresh = 'Huang',
                   filt = 'mean', seed = 9)
  expect_equal(round(mean(number, na.rm = TRUE), 1), 17.5)
  expect_error(number(img, "n", tau = "abc"), "If tau is a string")
  expect_error(number(img, "n", tau = FALSE),
               "tau must be specified as a positive number or")
  img %<>% ijtiff::read_tif()
  number <- number(img, "n", filt = "median")
  two_channel_img <- abind::abind(img, img, along = 3)
  number_2ch <- number(two_channel_img, "n", filt = "median")
  expect_equal(number_2ch %>% {list(dim(.), as.vector(.))},
               abind::abind(number, number, along = 3) %>%
               {list(dim(.), as.vector(.))},
               check.attributes = FALSE)
  expect_error(number(matrix(1:4, nrow = 2), "n"), "has dimension 2")
  expect_error(number(matrix(1:4, nrow = 2), "n"),
               "but has dimension 2")
})

test_that("number_time_series works", {
  context("number_time_series()")
  library(magrittr)
  img <- system.file('extdata', '50.tif', package = 'nandb')
  nts <- number_time_series(img, "N", 20, tau = 100, thresh = 'Huang',
                            filt = 'median', parallel = 2, seed = 9)
  expect_equal(round(mean(nts, na.rm = TRUE)), 18)
  nts <- number_time_series(img, "n", 30, tau = NA, thresh = 'tri',
                            filt = 'median', parallel = 2, seed = 99)
  expect_equal(round(median(nts, na.rm = TRUE), 3), -3.693)
  expect_error(number_time_series(img, "n", 51),
               paste("You have selected 51 frames per set,",
                     "but there are only 50 frames in total"))
  img <- ijtiff::read_tif(img)
  nts <- number_time_series(img, "N", 10)
  two_channel_img <- abind::abind(img, img, along = 3)
  nts_2ch <- number_time_series(two_channel_img, "N", 10)
  expect_equal(nts_2ch %>% {list(dim(.), as.vector(.))},
               abind::abind(nts, nts, along = 3) %>%
               {list(dim(.), as.vector(.))})
  expect_error(number_time_series(matrix(1:4, nrow = 2)),
               "argument.*def.*is missing, with no default")
  expect_error(number_time_series(matrix(1:4, nrow = 2), "N"),
               "Assertion.*failed.*Must have >=3 dimensions.*has dimension 2")
  setwd(tempdir())
  ijtiff::write_tif(img, '50.tif')
  ijtiff::write_tif(img, '50again.tif')
  number_time_series_folder(def = "n", tau = 333, thresh = 'tri',
                                frames_per_set = 20, seed = 0)
  expect_true(all(
    paste0("50", c("_number_n_timeseries_",
                   "again_number_n_timeseries_"),
           c("frames=20_tau=333_thresh=Triangle=0.68_filt=NA.tif",
             "frames=20_tau=333_thresh=Triangle=0.68_filt=NA.tif")) %in%
      list.files()))
  suppressWarnings(file.remove(list.files()))  # cleanup
})

test_that("number_folder works", {
  context("number_folder()")
  img <- ijtiff::read_tif(system.file('extdata', '50.tif', package = 'nandb'))
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(tempdir())
  ijtiff::write_tif(img, '50.tif')
  ijtiff::write_tif(img, '50again.tif')
  ijtiff::write_tif(array(4, dim = rep(3, 4)), "const.tif")
  number_folder(def = "N", tau = NA, seed = 3)
  expect_true(
    all(c("50_number_N_tau=NA_thresh=NA_filt=NA.tif",
          "50again_number_N_tau=NA_thresh=NA_filt=NA.tif",
          paste0("const_number_N_tau=NA,NA,NA_",
                 "thresh=NA,NA,NA_filt=NA,NA,NA.tif")) %in%
          list.files())
  )
  suppressWarnings(file.remove(list.files()))
})
