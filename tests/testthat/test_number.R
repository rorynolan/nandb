test_that("Number works", {
  img <- system.file('extdata', '50.tif', package = 'nandb')
  set.seed(33)
  number <- Number(img, tau = 'auto', mst = 'Huang', filt = 'median',
                           verbose = TRUE)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(round(mean(number, na.rm = TRUE), 4), 17.1119)
  }
  expect_error(Number(img, "abc"),
               "If tau is a string, it must be 'auto'.")
  expect_error(Number(img, FALSE),
               "If tau is not numeric, then it must be NA or 'auto'.")
  library(magrittr)
  img <- ReadImageData(img)
  number <- Number(img, filt = "smooth")
  two.channel.img <- abind::abind(img, img, along = 4) %>% aperm(c(1, 2, 4, 3))
  number.2ch <- Number(two.channel.img, n.ch = 2, filt = "smooth")
  expect_equal(number.2ch, abind::abind(number, number, along = 3))
  expect_error(Number(two.channel.img), "dimensional one")
  expect_error(Number(matrix(1:4, nrow = 2)), "dimensional one")
  expect_error(Number(matrix(1:4, nrow = 2), n.ch = 2),
               "nothing to be done")
})

test_that("NumberTimeSeries works", {
  img <- system.file('extdata', '50.tif', package = 'nandb')
  set.seed(333)
  nts <- NumberTimeSeries(img, 30, tau = 10, mst = 'Huang', filt = 'smooth',
                          mcc = 2, verbose = TRUE, seed = 4)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(round(mean(nts, na.rm = TRUE), 3), 11.98)
  }
  nts <- NumberTimeSeries(img, 20, tau = NA, mst = 'Huang',
                          mcc = 2, verbose = TRUE, seed = 4)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(round(mean(nts, na.rm = TRUE), 3), 12.504)
  }
  expect_error(NumberTimeSeries(img, 51),
               "frames.per.set must not be greater than the depth of arr3d")
  library(magrittr)
  img <- ReadImageData(img)
  two.channel.img <- abind::abind(img, img, along = 4) %>% aperm(c(1, 2, 4, 3))
  nts.2ch <- NumberTimeSeries(two.channel.img, 20, tau = NA, mst = 'Huang',
                              mcc = 2, verbose = TRUE,
                              seed = 4, n.ch = 2)
  expect_equal(nts.2ch,
               abind::abind(nts, nts, along = 4) %>% aperm(c(1, 2, 4, 3)))
  nts.2ch <- NumberTimeSeries(two.channel.img, 20, n.ch = 2)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(round(mean(nts.2ch, na.rm = TRUE), 3), 12.504)
  }
  expect_error(NumberTimeSeries(two.channel.img), "dimensional one")
  expect_error(NumberTimeSeries(matrix(1:4, nrow = 2)), "dimensional one")
  expect_error(NumberTimeSeries(matrix(1:4, nrow = 2), n.ch = 2),
               "nothing to be done")
})

test_that("NumberTxtFolder works", {
  img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(tempdir())
  WriteIntImage(img, '50.tif')
  WriteIntImage(img, '50again.tif')
  WriteIntImage(array(4, dim = rep(3, 3)), "const.tif")
  NumberTxtFolder(tau = NA, mst = NULL, mcc = 2, seed = 3)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
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
  numbers <- Numbers(img.paths, mst = 'Huang', tau = NA, mcc = 2,
                     seed = 7)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(round(mean(unlist(numbers), na.rm = TRUE), 3), 17.457)
  }
  imgs <- lapply(img.paths, ReadImageData)
  numbers <- Numbers(imgs, mst = 'Huang', tau = NA, mcc = 2, seed = 7)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(round(mean(unlist(numbers), na.rm = TRUE), 3), 17.457)
  }
  expect_error(Numbers(1:2, mst = "tri"),
               "must either be a list of 3d arrays")
})
