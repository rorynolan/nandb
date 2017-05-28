test_that("KmersFromBrightnesses works", {
  set.seed(3)
  brightnesses <- runif(100, 1, 3)
  monomer.med <- 1.2
  expect_equal(KmersFromBrightnesses(brightnesses, monomer.med),
               c(23, 21, 19, 21, 12), check.names = FALSE)
  expect_equal(KmersFromBrightnesses(t(0.5), 2), c(`1mers` = 0))
})


test_that("KmersFromImage works", {
  cwd <- getwd()
  on.exit(setwd(cwd))
  img <- system.file("extdata", "50.tif", package = "nandb")
  library(magrittr)
  expected <- c(50, 1) %>%
    set_names(c("1mers", "2mers")) %T>% {
      attr(., "mean.intensity") <- 15.0075
    }
  set.seed(3)
  expect_equal(KmersFromImage(img, 2.1, tau = NA, mst = "tri"),
               expected, tolerance = 1e-4, check.attributes = TRUE)
  expected <- list(expected %T>% {
                     .[1] <- 50
                     .[2] <- 1})[c(1, 1)]
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(KmersFromImages(c(img, img), 2.1, tau = NA, mst = "tri",
                                 seed = 8), expected, tolerance = 1e-4,
                 check.attributes = FALSE)
  }
  setwd(tempdir())
  file.copy(img, ".")
  KmersFromImagesFolder(monomer = 2.1, seed = 3)
  expect_equal(as.matrix(readr::read_csv("results.csv")),
               as.matrix(tibble::tibble(ImageName = "50.tif",
                                        MeanIntensity = 15.0075,
                                        `1mers` = 71L, `2mers` = 2L)))
  suppressWarnings(file.remove(list.files()))
})

test_that("KmerTIFFsFromBrightnessCSVs works", {
  cwd <- getwd()
  on.exit(setwd(cwd))
  img <- system.file("extdata", "50.tif", package = "nandb")
  setwd(tempdir())
  file.copy(img, ".")
  set.seed(6)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    BrightnessTxtFolder(tau = NA, mst = "tri", seed = 8)
    KmerTIFFsFromBrightnessCSVs(1.111)
    expect_equal(round(mean(ReadImageData(list.files(pattern = "tau.*tif"))), 2),
               0.8, tolerance = 0.01)
    expect_error(KmerTIFFsFromBrightnessCSVs(1.111, csv.paths = "a",
                                           out.names = c("a", "b")))
  }
  suppressWarnings(file.remove(list.files()))
})

test_that("KmerArray works", {
  set.seed(0)
  img <- ReadImageData(system.file("extdata", "50.tif",
                                   package = "nandb"))
  brightness <- Brightness(img, tau = NA, mst = "tri")
  ka <- KmerArray(brightness, 1.1)
  expect_equal(round(mean(ka), 3), 0.892)
  ka0 <- KmerArray(brightness, max(brightness + 1, na.rm = TRUE))
  expect_true(all(unique(ka0) %in% c(NA, 0)))
})
