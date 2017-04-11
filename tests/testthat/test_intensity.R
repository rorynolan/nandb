test_that("MeanIntensity works", {
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(tempdir())
  img <- ReadImageData(system.file("extdata", "50.tif",
    package = "nandb"))
  mean.intensity <- MeanIntensity(img, mst = "Huang",
    filt = "median")
  expect_equal(round(mean(mean.intensity, na.rm = TRUE), 3), 23.571)
  mean.intensity <- MeanIntensity(img, mst = "Huang",
                                  filt = "smooth")
  expect_equal(round(mean(mean.intensity, na.rm = TRUE), 3), 23.604)
  WriteIntImage(img, "50.tif")
  WriteIntImage(img, "50again.tif")
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    # these fail on windows due to an issue with readTIFF(..., as.is = TRUE)
    expect_equal(round(mean(unlist(MeanIntensityTxtFolder(mcc = 2))), 1),
      round(mean(unlist(MeanIntensities(list(img)[rep(1, 2)], mcc = 2))), 1))
    expect_equal(mean(ApplyOnPillars(img, mean)),
                 mean(unlist(MeanIntensities(list(img)[rep(1, 2)], mcc = 2))))
    expect_error(MeanIntensities(1:2, mst = "tri"))
    mean.intensities <- MeanIntensities(list(img, 2 * img), mst = "otsu")
    expect_equal(round(mean(unlist(mean.intensities), na.rm = TRUE), 3), 37.96)
    img <- system.file("extdata", "50.tif", package = "nandb")
    mean.intensities <- MeanIntensities(rep(img, 2), mst = "otsu")
    expect_equal(round(mean(unlist(mean.intensities), na.rm = TRUE), 3), 25.289)
  }
  suppressWarnings(file.remove(list.files()))
})
