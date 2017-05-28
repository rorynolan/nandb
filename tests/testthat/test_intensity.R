test_that("MeanIntensity works", {
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(tempdir())
  img <- ReadImageData(system.file("extdata", "50.tif",
    package = "nandb"))
  mean.intensity <- MeanIntensity(img)
  two.channel.img <- abind::abind(img, img, along = 4) %>% aperm(c(1, 2, 4, 3))
  mint.2ch <- MeanIntensity(two.channel.img)
  expect_equal(mint.2ch,
               abind::abind(mean.intensity, mean.intensity, along = 3))
  mean.intensity <- MeanIntensity(img, mst = "Huang", filt = "median")
  expect_equal(round(mean(mean.intensity, na.rm = TRUE), 3), 17.237)
  mean.intensity <- MeanIntensity(img, filt = "smooth")
  expect_equal(round(mean(mean.intensity, na.rm = TRUE), 3), 13.858)
  mean.intensity <- MeanIntensity(img, mst = "Huang", filt = "smooth")
  mint.2ch <- MeanIntensity(two.channel.img, mst = "h", filt = "s")
  expect_equal(mint.2ch,
               abind::abind(mean.intensity, mean.intensity, along = 3))
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
    expect_equal(round(mean(unlist(mean.intensities), na.rm = TRUE), 3), 26.659)
    img <- system.file("extdata", "50.tif", package = "nandb")
    mean.intensities <- MeanIntensities(rep(img, 2), mst = "otsu")
    expect_equal(round(mean(unlist(mean.intensities), na.rm = TRUE), 3), 17.752)
  }
  suppressWarnings(file.remove(list.files()))
})
