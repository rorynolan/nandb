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
