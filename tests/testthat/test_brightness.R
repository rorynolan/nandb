test_that("Brightness works", {
  img <- system.file('extdata', '50.tif', package = 'nandb')
  set.seed(33)
  brightness <- Brightness(img, tau = 'auto', mst = 'Huang', filt = 'median',
                           verbose = TRUE)
  expect_equal(round(mean(brightness, na.rm = TRUE), 4), 1.0314)
  brightness <- Brightness(img, tau = 'auto', mst = 'Huang', filt = 'smooth',
                           verbose = TRUE)
  expect_equal(round(mean(brightness, na.rm = TRUE), 4), 1.0448)
  expect_error(Brightness(img, "abc"),
               "If tau is a string, it must be 'auto'.")
  expect_error(Brightness(img, FALSE),
               "If tau is not numeric, then it must be NA or 'auto'.")
})

test_that("BrightnessTimeSeries works", {
  img <- system.file('extdata', '50.tif', package = 'nandb')
  set.seed(333)
  bts <- BrightnessTimeSeries(img, 10, tau = 'auto', mst = 'tri',
                              filt = 'median', mcc = 2, verbose = TRUE)
  expect_equal(round(mean(bts, na.rm = TRUE), 3), 0.959)
  bts <- BrightnessTimeSeries(img, 30, tau = 'auto', mst = 'tri',
                              filt = 'median', mcc = 2, verbose = TRUE)
  expect_equal(round(mean(bts, na.rm = TRUE), 3), 1.016)
  expect_error(BrightnessTimeSeries(img, 51),
               "frames.per.set must not be greater than the depth of arr3d")
})

test_that("BrightnessTxtFolder works", {
  img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(tempdir())
  WriteIntImage(img, '50.tif')
  WriteIntImage(img, '50again.tif')
  WriteIntImage(array(4, dim = rep(3, 3)), "const.tif")
  BrightnessTxtFolder(tau = NA, mst = 'tri', mcc = 2, seed = 3,
                      skip.consts = TRUE)
  expect_true(
    all(c("50_brightness_frames=50_tau=NA_mst=tri_filter=NA.csv",
          "50again_brightness_frames=50_tau=NA_mst=tri_filter=NA.csv") %in%
      list.files())
  )
  BrightnessPlotFolder()
  expect_true(
    all(c("50_brightness_frames=50_tau=NA_mst=tri_filter=NA.pdf",
          "50again_brightness_frames=50_tau=NA_mst=tri_filter=NA.pdf") %in%
          list.files())
  )
  suppressWarnings(file.remove(list.files()))
})

test_that("Brightnesses works", {
  img.paths <- rep(system.file('extdata', '50.tif', package = 'nandb'), 2)
  brightnesses <- Brightnesses(img.paths, mst = 'Huang', tau = 'auto', mcc = 2,
                               seed = 7)
  expect_equal(round(mean(unlist(brightnesses), na.rm = TRUE), 3), 1.045)
  expect_error(Brightnesses(1:2, mst = "tri"),
               "must either be a list of 3d arrays")
})
