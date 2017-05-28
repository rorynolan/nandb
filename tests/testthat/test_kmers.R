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
  expected <- c(66, 2) %>%
    set_names(c("1mers", "2mers")) %T>% {
      attr(., "mean.intensity") <- 15.8987
    }
  set.seed(3)
  expect_equal(KmersFromImage(img, 1.3, tau = NA, mst = "Huang"),
               expected, tolerance = 1e-4, check.attributes = TRUE)
  expected <- list(expected %T>% {
                     .[1] <- 66
                     .[2] <- 2})[c(1, 1)]
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(KmersFromImages(c(img, img), 1.3, tau = NA, mst = "Huang",
                                 seed = 8), expected, tolerance = 1e-4,
                 check.attributes = FALSE)
  }
  setwd(tempdir())
  file.copy(img, ".")
  KmersFromImagesFolder(monomer = 1.3, seed = 3)
  expect_equal(as.matrix(readr::read_csv("results.csv")),
               as.matrix(tibble::tibble(ImageName = "50.tif",
                                        MeanIntensity = 15.8987,
                                        `1mers` = 92L, `2mers` = 2L)))
  suppressWarnings(file.remove(list.files()))
})

test_that("KmerTIFFsFromBrightnessCSVs works", {
  cwd <- getwd()
  on.exit(setwd(cwd))
  img <- system.file("extdata", "50.tif", package = "nandb")
  setwd(tempdir())
  file.copy(img, ".")
  set.seed(6)
  BrightnessTxtFolder(tau = NA, mst = "Huang", seed = 8)
  KmerTIFFsFromBrightnessCSVs(1.111)
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(round(mean(ReadImageData(list.files(pattern = "tau.*tif"))), 2),
               0.57, tolerance = 0.01)
    expect_error(KmerTIFFsFromBrightnessCSVs(1.111, csv.paths = "a",
                                           out.names = c("a", "b")))
  }
  suppressWarnings(file.remove(list.files()))
})

test_that("KmerArray works", {
  set.seed(0)
  img <- ReadImageData(system.file("extdata", "50.tif",
                                   package = "nandb"))
  brightness <- Brightness(img, tau = NA, mst = "Huang")
  ka <- KmerArray(brightness, 1.1)
  expect_equal(round(mean(ka), 3), 0.625)
  ka0 <- KmerArray(brightness, max(brightness + 1, na.rm = TRUE))
  expect_true(all(unique(ka0) %in% c(NA, 0)))
})
