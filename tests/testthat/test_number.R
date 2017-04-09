test_that("Number works", {
  img <- system.file('extdata', '50.tif', package = 'nandb')
  set.seed(33)
  number <- Number(img, tau = 'auto', mst = 'Huang', filt = 'median',
                           verbose = TRUE)
  expect_equal(round(mean(number, na.rm = TRUE), 4), 22.9993)
  number <- Number(img, tau = 'auto', mst = 'Huang', filt = 'smooth',
                           verbose = TRUE)
  expect_equal(round(mean(number, na.rm = TRUE), 4), 23.5807)
  expect_error(Number(img, "abc"),
               "If tau is a string, it must be 'auto'.")
  expect_error(Number(img, FALSE),
               "If tau is not numeric, then it must be NA or 'auto'.")
})

test_that("NumberTimeSeries works", {
  img <- system.file('extdata', '50.tif', package = 'nandb')
  set.seed(333)
  bts <- NumberTimeSeries(img, 10, tau = 'auto', mst = 'tri',
                              filt = 'median', mcc = 2, verbose = TRUE)
  expect_equal(round(mean(bts, na.rm = TRUE), 3), 23.545)
  bts <- NumberTimeSeries(img, 30, tau = 'auto', mst = 'tri',
                              filt = 'median', mcc = 2, verbose = TRUE)
  expect_equal(round(mean(bts, na.rm = TRUE), 3), 21.043)
  expect_error(NumberTimeSeries(img, 51),
               "frames.per.set must not be greater than the depth of arr3d")
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
    NumberTxtFolder(tau = NA, mst = 'tri', mcc = 2, seed = 3,
                        skip.consts = TRUE)
    expect_true(
      all(c("50_number_frames=50_tau=NA_mst=tri_filter=NA.csv",
            "50again_number_frames=50_tau=NA_mst=tri_filter=NA.csv") %in%
            list.files())
    )
  }
  suppressWarnings(file.remove(list.files()))
})

test_that("Numbers works", {
  img.paths <- rep(system.file('extdata', '50.tif', package = 'nandb'), 2)
  numbers <- Numbers(img.paths, mst = 'Huang', tau = 'auto', mcc = 2,
                               seed = 7)
  expect_equal(round(mean(unlist(numbers), na.rm = TRUE), 3), 23.716)
  imgs <- lapply(img.paths, ReadImageData)
  numbers <- Numbers(imgs, mst = 'Huang', tau = 'auto', mcc = 2, seed = 7)
  expect_equal(round(mean(unlist(numbers), na.rm = TRUE), 3), 23.716)
  expect_error(Numbers(1:2, mst = "tri"),
               "must either be a list of 3d arrays")
})
