# Graphics are difficult to test well, and I don't manage it here.
# They still require a manual check, even if they are kind of tested here.

test_that("ArrArrHexPlot is ggplot", {
  img <- ReadImageData(system.file("extdata", "50.tif",
    package = "nandb"))
  brightness <- Brightness(img)
  mean.intensity <- MeanPillars(img)
  expect_is(ArrArrHexPlot(mean.intensity, brightness), "ggplot")
})

test_that("KmerPlot is ggplot", {
  img <- ReadImageData(system.file("extdata", "50.tif",
    package = "nandb"))
  brightness <- Brightness(img)
  expect_is(KmerPlot(brightness, 1.02), "ggplot")
  expect_is(KmerPlot(brightness, 100), "ggplot")
  expect_is(KmerPlot(brightness, 1.02, log.trans = TRUE), "ggplot")
  expect_is(KmerPlot(brightness, 1.52), "ggplot")
})

test_that("MatrixRasterPlot is ggplot", {
  img <- ReadImageData(system.file("extdata", "50.tif",
    package = "nandb"))
  brightness <- Brightness(img)
  expect_is(MatrixRasterPlot(brightness, scale.name = "brightness"), "ggplot")
  expect_is(MatrixRasterPlot(brightness, scale.name = "brightness",
    log.trans = TRUE), "ggplot")
  expect_is(MatrixRasterPlot(brightness, scale.name = "brightness",
    log.trans = TRUE, breaks = 1:3), "ggplot")
  expect_is(MatrixRasterPlot(brightness, scale.name = "brightness",
    ranges = seq(0.5, 3, length.out = 6), range.names = paste0(1:5,
      "mer")), "ggplot")
  expect_is(MatrixRasterPlot(brightness, scale.name = "brightness",
    ranges = seq(0.5, 3, length.out = 6),
    range.names = paste0(1:5, "mer"), log.trans = TRUE), "ggplot")
  expect_error(MatrixRasterPlot(brightness, scale.name = "brightness",
    ranges = seq(0.5, 3, length.out = 6),
    range.names = paste0(1:59, "mer"), log.trans = TRUE))
  expect_error(MatrixRasterPlot(brightness, scale.name = "brightness",
    ranges = seq(0.5, 3, length.out = 6), colours = viridis::viridis(999),
    range.names = paste0(1:59, "mer"), log.trans = TRUE))
  expect_is(MatrixRasterPlot(brightness, scale.name = "brightness",
    ranges = seq(0.5, 3, length.out = 6), include.breaks = 1.25,
    range.names = NULL, log.trans = TRUE), "ggplot")
  expect_is(MatrixRasterPlot(brightness, scale.name = "brightness",
    ranges = seq(0.5, 3, length.out = 6), include.breaks = 1.25,
    range.names = NULL, log.trans = FALSE), "ggplot")
  expect_is(MatrixRasterPlot(brightness, scale.name = "brightness",
                             limits = c(1, 1.25), clip = TRUE), "ggplot")
  expect_is(MatrixRasterPlot(brightness, scale.name = "brightness",
                             include.breaks = 1.25), "ggplot")
  expect_is(MatrixRasterPlot(brightness, scale.name = "brightness",
                             include.breaks = 1.25, log.trans = TRUE), "ggplot")
  expect_error(MatrixRasterPlot(brightness, scale.name = "brightness",
                                limits = c(-1, 1), log.trans = TRUE))
})
