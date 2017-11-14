context("Brightness")

test_that("brightness works", {
  img <- system.file('extdata', '50.tif', package = 'nandb')
  brightness <- brightness(img, "e", tau = 'auto', thresh = 'Huang',
                           filt = 'median', seed = 0, parallel = 2)
  if (not_windows()) {
    expect_equal(round(mean(brightness, na.rm = TRUE), 4), 0.0242)
  }
  brightness <- brightness(img, "B", tau = 1000, thresh = 'Huang',
                           filt = 'mean', seed = 9)
  if (not_windows()) {
    expect_equal(round(mean(brightness, na.rm = TRUE), 4), 1.0433)
  }
  expect_error(brightness(img, "e", tau = "abc"), "If tau is a string")
  expect_error(brightness(img, "B", tau = FALSE),
               "tau must be specified as a positive number or")
  img %<>% read_tif()
  brightness <- brightness(img, "B", filt = "median")
  two_channel_img <- abind::abind(img, img, along = 4) %>% aperm(c(1, 2, 4, 3))
  brightness_2ch <- brightness(two_channel_img, "B", n_ch = 2, filt = "median")
  expect_equal(brightness_2ch %>% {list(dim(.), as.vector(.))},
               abind::abind(brightness, brightness, along = 3) %>%
                 {list(dim(.), as.vector(.))},
               check.attributes = FALSE)
  expect_error(brightness(two_channel_img, "e"),
               "as having 1 channels, but it looks like it has 2")
  expect_error(brightness(matrix(1:4, nrow = 2), "B"), "has dimension 2")
  expect_error(brightness(matrix(1:4, nrow = 2), "B", n_ch = 2),
               "but has dimension 2")
})

 test_that("brightness_time_series works", {
  library(magrittr)
  img <- system.file('extdata', '50.tif', package = 'nandb')
  bts <- brightness_time_series(img, "e", 20, tau = 100, thresh = 'Huang',
                                filt = 'median', parallel = 2, seed = 9)
  if (not_windows()) {
    expect_equal(round(mean(bts, na.rm = TRUE), 3), -0.018)
  }
  bts <- brightness_time_series(img, "B", 30, tau = NA, thresh = 'tri',
                                filt = 'median', parallel = 2, seed = 99)
  if (not_windows()) {
    expect_equal(round(mean(bts, na.rm = TRUE), 3), 1.008)
  }
  expect_error(brightness_time_series(img, "b", 51),
               paste("You have selected 51 frames per set,",
                     "but there are only 50 frames in total"))
  img <- read_tif(img)
  bts <- brightness_time_series(img, "b", 10)
  two_channel_img <- abind::abind(img, img, along = 4) %>% aperm(c(1, 2, 4, 3))
  bts_2ch <- brightness_time_series(two_channel_img, "b", 10, n_ch = 2)
  expect_equal(bts_2ch %>% {list(dim(.), as.vector(.))},
               abind::abind(bts, bts, along = 4) %>% aperm(c(1, 2, 4, 3)) %>%
                 {list(dim(.), as.vector(.))})
  expect_error(brightness_time_series(two_channel_img, "e"),
               "as having 1 channels, but it looks like it has 2")
  expect_error(brightness_time_series(matrix(1:4, nrow = 2)),
               "argument.*def.*is missing, with no default")
  expect_error(brightness_time_series(matrix(1:4, nrow = 2), "b", n_ch = 2),
               "Assertion.*failed.*Must have >=3 dimensions.*has dimension 2")
  setwd(tempdir())
  write_tif(img, '50.tif')
  write_tif(img, '50again.tif')
  brightness_time_series_folder(def = "B", tau = 333, thresh = 'tri',
                                frames_per_set = 20, seed = 0)
  if (not_windows()) {
    expect_true(all(
      paste0("50", c("_brightness_B_timeseries_",
                     "again_brightness_B_timeseries_"),
             c("frames=20_tau=333_thresh=Triangle=0.68_filt=NA_1.txt",
               "frames=20_tau=333_thresh=Triangle=0.68_filt=NA_2.txt")) %in%
        list.files()))
  }
  suppressWarnings(file.remove(list.files()))  # cleanup
})

test_that("brightness_folder works", {
  img <- read_tif(system.file('extdata', '50.tif', package = 'nandb'))
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(tempdir())
  write_tif(img, '50.tif')
  write_tif(img, '50again.tif')
  write_tif(array(4, dim = rep(3, 3)), "const.tif")
  brightness_folder(def = "B", tau = NA, seed = 3)
  expect_true(
    all(c("50_brightness_B_tau=NA_thresh=NA_filt=NA.txt",
          "50again_brightness_B_tau=NA_thresh=NA_filt=NA.txt") %in%
          list.files())
  )
  suppressWarnings(file.remove(list.files()))
})
