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
  expected <- c(192, 1) %>%
    set_names(c("1mers", "2mers")) %T>% {
      attr(., "mean.intensity") <- 21.11217
    }
  set.seed(3)
  expect_equal(KmersFromImage(img, 2.1, tau = "auto", mst = "huang"),
               expected, tolerance = 1e-5, check.attributes = TRUE)
  expected <- list(expected %T>% {.[1] <- 198},
                   expected %T>% {
                     . <- .[-2]
                     .[1] <- 175
                   }) %>% rev
  if (!stringr::str_detect(tolower(Sys.info()['sysname']), "windows")) {
    expect_equal(KmersFromImages(c(img, img), 2.1, tau = "auto", mst = "huang",
                                 seed = 8), expected, tolerance = 1e-5,
                 check.attributes = FALSE)
  }
  setwd(tempdir())
  file.copy(img, ".")
  KmersFromImagesFolder(monomer = 2.1, seed = 3)
  expect_equal(readr::read_csv("results.csv"),
               tibble::tibble(ImageName = "50.tif", MeanIntensity = 21.112174,
                              `1mers` = 236L, `2mers` = 2L),
               tolerance = 1e-4)
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
    BrightnessTxtFolder(tau = "auto", mst = "tri", seed = 8)
    KmerTIFFsFromBrightnessCSVs(1.111)
    expect_equal(round(mean(ReadImageData(list.files(pattern = "tau.*tif"))), 2),
               0.85, tolerance = 0.1)
    expect_error(KmerTIFFsFromBrightnessCSVs(1.111, csv.paths = "a",
                                           out.names = c("a", "b")))
  }
  suppressWarnings(file.remove(list.files()))
})

test_that("KmerArray works", {
  img <- ReadImageData(system.file("extdata", "50.tif",
                                   package = "nandb"))
  brightness <- Brightness(img, tau = "auto", mst = "Huang",
                           filt = "median")
  ka <- KmerArray(brightness, 1.1)
  expect_equal(round(mean(ka), 4), 0.4725)
  ka0 <- KmerArray(brightness, max(brightness + 1, na.rm = TRUE))
  expect_true(all(unique(ka0) %in% c(NA, 0)))
})
