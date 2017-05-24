test_that("ForceChannels works", {
  library(magrittr)
  x <- lapply(1:300, function(x) matrix(runif(4), nrow = 2)) %>%
    Reduce(function(x, y) abind::abind(x, y, along = 3),
      .)
  expect_equal(ForceChannels(x, 6) %>% dim, c(2, 2, 6, 50))
  expect_error(ForceChannels(x, 7), "multiple")
})

test_that("Stack2DTifs works", {
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(tempdir())
  img <- array(0:8, dim = rep(2, 3))
  WriteIntImage(img[, , 1], "50_1.tif")
  WriteIntImage(img[, , 2], "50_2.tif")
  Stack2DTifs(c("50_1.tif", "50_2.tif"), "50_1_2")
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    # these fail on windows due to an issue with readTIFF(..., as.is = TRUE)
    expect_equal(ReadImageData("50_1_2.tif"),
                 abind::abind(img[, , 1], img[, , 2], along = 3),
                 check.attributes = FALSE)
  }
  suppressWarnings(file.remove(list.files()))
  img <- ReadImageData(system.file("extdata", "50.tif",
                                   package = "nandb"))
  WriteIntImage(img[, , 1], "50_1.tif")
  WriteIntImage(img[, , 1:2], "50_2.tif")
  expect_error(Stack2DTifs(c("50_1.tif", "50_2.tif"), "50_1_2"), "same dim")
  suppressWarnings(file.remove(list.files()))
  img <- ReadImageData(system.file("extdata", "50.tif",
                                   package = "nandb"))
  WriteIntImage(img[, , 1:2], "50_1.tif")
  WriteIntImage(img[, , 1:2], "50_2.tif")
  expect_error(Stack2DTifs(c("50_1.tif", "50_2.tif"), "50_1_2"), "2-dim")
  suppressWarnings(file.remove(list.files()))
})

test_that("WriteIntImage works", {
  img <- ReadImageData(system.file("extdata", "50.tif",
    package = "nandb"))
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(tempdir())
  WriteIntImage(img, "50.tif")
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    # these fail on windows due to an issue with readTIFF(..., as.is = TRUE)
    expect_equal(ReadImageData("50.tif"), img, check.attributes = FALSE)
  }
  suppressWarnings(file.remove(list.files()))
})

test_that("ReadImageData works", {
  file.path <- system.file('extdata', '50.tif', package = 'nandb')
  expect_equal(EBImage::imageData(EBImage::readImage(file.path, as.is = TRUE)),
               ReadImageData(file.path, 3))
  expect_error(ReadImageData(file.path, TRUE), "integer")

})

test_that("WriteImageTxt works", {
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(tempdir())
  img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
  expect_equal(mean(unlist(WriteImageTxt(img, 'temp'))), mean(img))
  expect_error(WriteImageTxt(1:3, "abc"), "dimension")
  expect_error(WriteIntImage(matrix(0.5), "a"), "integer")
  expect_error(WriteIntImage(matrix(-1), "a"), "negative")
  WriteIntImage(matrix(2^9), "16bit")
  expect_error(WriteIntImage(matrix(2 ^ 17)),
               "The maximum value")
  expect_equal(mean(unlist(WriteImageTxt(img, 'temp'))), mean(img))
  img_01 <- ReadImageTxt("temp_01.csv")
  expect_equal(img_01, img[, , 1], check.attributes = FALSE)
  four.d <- array(1:(2^4), dim = rep(2, 4))
  WriteImageTxt(four.d, "fourD")
  for (i in 1:2) {
    for (j in 1:2) {
      expect_equal(ReadImageTxt(paste0("fourD_", i, "_", j, ".csv")),
                   four.d[, , i, j], check.attributes = FALSE)
    }
  }
  suppressWarnings(file.remove(list.files()))
})

test_that("Bin2Tiff works", {
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(tempdir())
  dir.create("temp_dir")
  expect_true(file.copy(system.file("extdata", "b.bin", package = "nandb"),
                         "temp_dir"))
  Bin2Tiff("temp_dir/b.bin", height = 2, width = 2, bits = 8)
  Bin2TiffFolder("temp_dir", height = 2, width = 2, bits = 8)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    # these fail on windows due to an issue with readTIFF(..., as.is = TRUE)
    expect_equal(list.files("temp_dir"), c("b.bin", "b.tif"))
    setwd("temp_dir")
    expect_equal(readBin("b.bin", "int", size = 1, n = 4),
                 as.vector(ReadImageData("b.tif")))
  }
  setwd("..")
  filesstrings::RemoveDirs("temp_dir")
})
