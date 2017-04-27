test_that("Number works", {
  img <- system.file('extdata', '50.tif', package = 'nandb')
  set.seed(33)
  number <- Number(img, tau = 'auto', mst = 'Huang', filt = 'median',
                           verbose = TRUE)
  expect_equal(round(mean(number, na.rm = TRUE), 4), 22.9993)
  number <- Number(img, tau = 1000, mst = 'Huang', filt = 'smooth',
                   verbose = TRUE)
  expect_equal(round(mean(number, na.rm = TRUE), 4), 23.5807)
  expect_error(Number(img, "abc"),
               "If tau is a string, it must be 'auto'.")
  expect_error(Number(img, FALSE),
               "If tau is not numeric, then it must be NA or 'auto'.")
  library(magrittr)
  img <- ReadImageData(img)
  number <- Number(img)
  two.channel.img <- abind::abind(img, img, along = 4) %>% aperm(c(1, 2, 4, 3))
  number.2ch <- Number(two.channel.img)
  expect_equal(number.2ch, abind::abind(number, number, along = 3))
})

test_that("NumberTimeSeries works", {
  img <- system.file('extdata', '50.tif', package = 'nandb')
  set.seed(333)
  nts <- NumberTimeSeries(img, 30, tau = 10, mst = 'tri', filt = 'median',
                          mcc = 2, verbose = TRUE, seed = 4)
  expect_equal(round(mean(nts, na.rm = TRUE), 3), 22.551)
  nts <- NumberTimeSeries(img, 10, tau = 100, mst = 'tri',
                          filt = 'median', mcc = 2, verbose = TRUE, seed = 4)
  expect_equal(round(mean(nts, na.rm = TRUE), 3), 23.594)
  expect_error(NumberTimeSeries(img, 51),
               "frames.per.set must not be greater than the depth of arr3d")
  library(magrittr)
  img <- ReadImageData(img)
  two.channel.img <- abind::abind(img, img, along = 4) %>% aperm(c(1, 2, 4, 3))
  nts.2ch <- NumberTimeSeries(two.channel.img, 10, tau = 100, mst = 'tri',
                              filt = 'median', mcc = 2, verbose = TRUE,
                              seed = 4)
  expect_equal(nts.2ch,
               abind::abind(nts, nts, along = 4) %>% aperm(c(1, 2, 4, 3)))
  nts.2ch <- NumberTimeSeries(two.channel.img, 10)
  expect_equal(round(mean(nts.2ch, na.rm = TRUE), 3), 27.382)
})

test_that("NumberTxtFolder works", {
  img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(tempdir())
  WriteIntImage(img, '50.tif')
  WriteIntImage(img, '50again.tif')
  WriteIntImage(array(4, dim = rep(3, 3)), "const.tif")
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    NumberTxtFolder(tau = NA, mst = NULL, mcc = 2, seed = 3)
    expect_true(
      all(c("50_number_frames=50_tau=NA_mst=NA_filter=NA.csv",
            "50again_number_frames=50_tau=NA_mst=NA_filter=NA.csv") %in%
            list.files())
    )
  }
  suppressWarnings(file.remove(list.files()))
})

test_that("Numbers works", {
  img.paths <- rep(system.file('extdata', '50.tif', package = 'nandb'), 2)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    numbers <- Numbers(img.paths, mst = 'Huang', tau = 7, mcc = 2,
                       seed = 7)
    expect_equal(round(mean(unlist(numbers), na.rm = TRUE), 3), 26.634)
    imgs <- lapply(img.paths, ReadImageData)
    numbers <- Numbers(imgs, mst = 'Huang', tau = 8, mcc = 2, seed = 7)
    expect_equal(round(mean(unlist(numbers), na.rm = TRUE), 3), 26.227)
  }
  expect_error(Numbers(1:2, mst = "tri"),
               "must either be a list of 3d arrays")
})
