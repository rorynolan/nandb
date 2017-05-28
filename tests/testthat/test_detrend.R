test_that("CorrectForBleaching works", {
  img <- ReadImageData(system.file("extdata", "50.tif", package = "nandb"))
  set.seed(9)
  autotau <- CorrectForBleaching(img, "auto")
  expect_equal(round(mean(autotau), 4), 15.0092)
  expect_error(CorrectForBleaching(img, "abc"),
               "If tau is a string, it must be 'auto'.")
  expect_error(CorrectForBleaching(img, FALSE),
               "If tau is not numeric, then it must be NA or 'auto'.")
  library(magrittr)
  two.channel.img <- abind::abind(img, img, along = 4) %>% aperm(c(1, 2, 4, 3))
  twoch <- CorrectForBleaching(two.channel.img, 100)
  expect_equal(mean(twoch), mean(CorrectForBleaching(img, 100)))
})

test_that("CorrectForBleachingFolder works", {
  cwd <- getwd()
  setwd(tempdir())
  on.exit(setwd(cwd))
  img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
  WriteIntImage(img, '50.tif')
  WriteIntImage(img, '50again.tif')
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    set.seed(19)
    CorrectForBleachingFolder(tau = "auto", mst = 'tri', mcc = 1, na = "s",
                              seed = 8)
    expect_true(all(c("50_tau=auto=155_mst=tri.tif",
                      "50again_tau=auto=144_mst=tri.tif") %in%
                      list.files()))
  }
  expect_error(nandb:::CorrectForBleachingFile("50.tif", mst = "h"),
               "If you select thresholding")
  suppressWarnings(file.remove(list.files()))
})

test_that("BestTau works", {
  set.seed(5)
  img <- system.file('extdata', '50.tif', package = 'nandb')
  expect_true(round(BestTau(img, mst = 'tri', tol = 3)) == 100)
  expect_error(BestTau(array(1, dim = rep(3, 3))),
               "Your raw brightness mean is below 1,")
  savage <- abind::abind(matrix(0, nrow = 2, ncol = 2),
                         matrix(100, nrow = 2, ncol = 2),
                         matrix(0, nrow = 2, ncol = 2), along = 3)
  expect_error(BestTau(savage), "savage")
})

test_that("ExpSmoothNaive works", {
  expect_equal(round(nandb:::ExpSmoothNaive(1:3, 2), 2), c(1.68, 2, 2.32))
})
