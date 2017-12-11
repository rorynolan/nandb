context("Class constructors")

test_that("Class construction edge cases function correctly", {
  img <- array(runif(8), dim = rep(2, 3))
  eg <- number_img(img, "N", 4, radiant.data::set_attr(8, "auto", TRUE), NA)
  expect_equal(attr(eg, "thresh"), rep(4, 2))
  expect_equal(attr(eg, "tau"),
               radiant.data::set_attr(c(8, 8), "auto", c(T, T)))
  expect_error(number_img(img, "n", 4, 5, NA), "tau.*spec.*must.*attr.*auto")
  expect_error(brightness_img(img, "a", 4, 9, NA), "def.*must.*one of")
  expect_error(brightness_ts_img(img, "a", 10, 4, 9, NA), "def.*must.*one of")
  x <- list(radiant.data::set_attr(4, "a", 0),
            radiant.data::set_attr(5, "b", 0))
  ans <- 4:5
  attributes(ans) <- list(a = c(0, NA), b = c(NA, 0))
  expect_equal(c_list_attr_na(x), ans)
  img <- array(runif(3 ^ 3), dim = rep(3, 3))
  expect_error(number_img(img, "n", 5,
                          radiant.data::set_attr(4:6, "auto", rep(FALSE, 2)),
                          NA), "auto.*attribute.*tau.*same length.*tau.*itself")
  expect_equal(number_img(img, "n", 5,
                          radiant.data::set_attr(4:6, "auto", FALSE), NA),
               number_img(img, "n", 5, radiant.data::set_attr(4:6, "auto",
                                                              rep(FALSE, 3)),
                          NA))
  expect_error(number_img(img, "n", 5,
                          radiant.data::set_attr(4:6, "auto", NA), NA),
               "tau.*att.*not NA")
  expect_error(number_img(img, "n", 5,
                          radiant.data::set_attr(4:5, "auto", rep(FALSE, 2)),
                          NA), "thresh.*tau.*filt.*same.*channels")
})
