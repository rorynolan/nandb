test_that("CorrectForBleaching works", {
  img <- ReadImageData(system.file("extdata", "50.tif", package = "nandb"))
  tau10 <- CorrectForBleaching(img, 10)
  expect_equal(round(mean(tau10), 4), 21.1089)
  set.seed(9)
  autotau <- CorrectForBleaching(img, "auto")
  expect_equal(round(mean(autotau), 4), 21.1122)
  expect_error(CorrectForBleaching(img, "abc"),
               "If tau is a string, it must be 'auto'.")
  expect_error(CorrectForBleaching(img, FALSE),
               "If tau is not numeric, then it must be NA or 'auto'.")
})

test_that("CorrectForBleachingFolder works", {
  cwd <- getwd()
  setwd(tempdir())
  on.exit(setwd(cwd))
  img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
  WriteIntImage(img, '50.tif')
  WriteIntImage(img, '50again.tif')
  CorrectForBleachingFolder(tau = 'auto', mst = 'tri', mcc = 2, na = "s",
                            seed = 8)
  expect_true(all(c("50_tau=auto=268.tif", "50again_tau=auto=NA.tif") %in%
                    list.files()))
  suppressWarnings(file.remove(list.files()))
})

test_that("BestTau works", {
  set.seed(3)
  img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
  expect_true(round(BestTau(img, mst = 'tri', tol = 3)) == 32480)
  set.seed(3)
  img <- system.file('extdata', '50.tif', package = 'nandb')
  expect_true(round(BestTau(img, mst = 'tri', tol = 3)) == 32480)
  expect_error(BestTau(array(1, dim = rep(3, 3))),
               "Your raw brightness mean is below 1,")
  savage <- abind::abind(matrix(0, nrow = 2, ncol = 2),
                         matrix(100, nrow = 2, ncol = 2),
                         matrix(0, nrow = 2, ncol = 2), along = 3)
  expect_error(BestTau(savage), "savage")
})
