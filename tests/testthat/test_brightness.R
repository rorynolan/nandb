test_that("Brightness works", {
  library(magrittr)
  img <- system.file('extdata', '50.tif', package = 'nandb')
  set.seed(33)
  brightness <- Brightness(img, tau = 'auto', mst = 'Huang', filt = 'median',
                           verbose = TRUE)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(round(mean(brightness, na.rm = TRUE), 4), 1.0239)
  }
  brightness <- Brightness(img, tau = 1000, mst = 'Huang', filt = 'smooth',
                           verbose = TRUE)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(round(mean(brightness, na.rm = TRUE), 4), 1.0442)
  }
  expect_error(Brightness(img, "abc"),
               "If tau is a string, it must be 'auto'.")
  expect_error(Brightness(img, FALSE),
               "If tau is not numeric, then it must be NA or 'auto'.")
  img <- ReadImageData(img)
  brightness <- Brightness(img, filt = "smooth")
  two.channel.img <- abind::abind(img, img, along = 4) %>% aperm(c(1, 2, 4, 3))
  brightness.2ch <- Brightness(two.channel.img, n.ch = 2, filt = "smooth")
  expect_equal(brightness.2ch, abind::abind(brightness, brightness, along = 3))
  expect_error(Brightness(two.channel.img), "dimensional one")
  expect_error(Brightness(matrix(1:4, nrow = 2)), "dimensional one")
  expect_error(Brightness(matrix(1:4, nrow = 2), n.ch = 2),
               "nothing to be done")
})

 test_that("BrightnessTimeSeries works", {
  library(magrittr)
  img <- system.file('extdata', '50.tif', package = 'nandb')
  set.seed(333)
  bts <- BrightnessTimeSeries(img, 20, tau = 100, mst = 'Huang',
                              filt = 'median', mcc = 2, seed = 9,
                              verbose = TRUE)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(round(mean(bts, na.rm = TRUE), 3), 0.977)
  }
  bts <- BrightnessTimeSeries(img, 30, tau = NA, mst = 'tri',
                              filt = 'median', mcc = 2, seed = 99,
                              verbose = TRUE)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(round(mean(bts, na.rm = TRUE), 3), 1.008)
  }
  expect_error(BrightnessTimeSeries(img, 51),
               "frames.per.set must not be greater than the depth of arr3d")
  img <- ReadImageData(img)
  bts <- BrightnessTimeSeries(img, 10)
  two.channel.img <- abind::abind(img, img, along = 4) %>% aperm(c(1, 2, 4, 3))
  bts.2ch <- BrightnessTimeSeries(two.channel.img, 10, n.ch = 2)
  expect_equal(bts.2ch,
               abind::abind(bts, bts, along = 4) %>% aperm(c(1, 2, 4, 3)))
  expect_error(BrightnessTimeSeries(two.channel.img), "dimensional one")
  expect_error(BrightnessTimeSeries(matrix(1:4, nrow = 2)), "dimensional one")
  expect_error(BrightnessTimeSeries(matrix(1:4, nrow = 2), n.ch = 2),
               "nothing to be done")
  setwd(tempdir())
  WriteIntImage(img, '50.tif')
  WriteIntImage(img, '50again.tif')
  BrightnessTimeSeriesTxtFolder(tau = 333, mst = 'tri', mcc = 2,
                                frames.per.set = 20, seed = 0)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_true(all(
      c("50_brightness_frames=20_tau=333_mst=tri_filter=NA_1.csv",
        "50_brightness_frames=20_tau=333_mst=tri_filter=NA_2.csv") %in%
        list.files()))
  }
  suppressWarnings(file.remove(list.files()))  # cleanup
})

test_that("BrightnessTxtFolder works", {
  img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(tempdir())
  WriteIntImage(img, '50.tif')
  WriteIntImage(img, '50again.tif')
  WriteIntImage(array(4, dim = rep(3, 3)), "const.tif")
  BrightnessTxtFolder(tau = NA, mst = NULL, mcc = 2, seed = 3)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_true(
      all(c("50_brightness_frames=50_tau=NA_mst=NA_filter=NA.csv",
            "50again_brightness_frames=50_tau=NA_mst=NA_filter=NA.csv") %in%
        list.files())
    )
    BrightnessPlotFolder()
    expect_true(
      all(c("50_brightness_frames=50_tau=NA_mst=NA_filter=NA.pdf",
            "50again_brightness_frames=50_tau=NA_mst=NA_filter=NA.pdf") %in%
            list.files())
    )
  }
  suppressWarnings(file.remove(list.files()))
})

test_that("Brightnesses works", {
  img.paths <- rep(system.file('extdata', '50.tif', package = 'nandb'), 2)
  brightnesses <- Brightnesses(img.paths, mst = 'Huang', tau = 250, mcc = 2,
                               seed = 7)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(round(mean(unlist(brightnesses), na.rm = TRUE), 3), 1.039)
  }
  expect_error(Brightnesses(1:2, mst = "tri"),
               "must either be a list of 3d arrays")
})

test_that("BrightnessTimeSeriess works", {
  img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
  arr3d.list <- list(img, img)
  btss <- nandb:::BrightnessTimeSeriess(arr3d.list, 20)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(round(purrr::map_dbl(btss, mean, na.rm = TRUE), 4),
                 rep(1.0124, 2))
  }
  btss <- nandb:::BrightnessTimeSeriess(arr3d.list, 20, mst = "otsu")
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(round(purrr::map_dbl(btss, mean, na.rm = TRUE), 4),
                 rep(1.0105, 2))
  }
  expect_error(nandb:::BrightnessTimeSeriess(1:2, 10, mst = "tri"),
               "either.*vector of path")
})
