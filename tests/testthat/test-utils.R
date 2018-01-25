context("Utils")

test_that("nb_get_img() errors correctly", {
  expect_error(nb_get_img(c("a", "b")), "char.*vec.*length 1")
  img <- array(runif(8), dim = rep(2, 3))
  expect_error(nb_get_img(img), "pos.*int.*NA.*only")
})

test_that("extend_for_all_chs() edge case works", {
  expect_equal(extend_for_all_chs(list(NULL), 4), as.list(rep(NA, 4)))
})

test_that("fix_filt() edge cases and exceptions work correctly", {
  expect_equal(fix_filt(NULL), NA_character_)
  expect_error(fix_filt("abc"), "must be either")
})

test_that("deduplicate_nb_filename() works correctly", {
  path <- paste0("detrended_exponential_thresh=4,3_tau=1,2",
                 "_brightness_thresh=NA,NA_tau=NA,NA_filt=NA,NA.tif")
  expect_equal(deduplicate_nb_filename(path),
               paste0("detrended_exponential_thresh=4,3_tau=1,2",
                      "_brightness_filt=NA,NA.tif"))
  path <- paste0("detrended_exponential_thresh=4,3_tau=1,2",
                 "_brightness_thresh=NA,NA_tau=NA,NA_")
  expect_equal(deduplicate_nb_filename(path),
               paste0("detrended_exponential_thresh=4,3_tau=1,2",
                      "_brightness"))
  path <- paste0("detrended_exponential_thresh=4,3_tau=1,2",
                 "_brightness_thresh=NA,NA_")
  expect_equal(deduplicate_nb_filename(path),
               paste0("detrended_exponential_thresh=4,3_tau=1,2",
                      "_brightness"))
  expect_equal(deduplicate_nb_filename("abc"), "abc")
})
