test_that("ExpSmooth works", {
  expect_equal(round(ExpSmooth(1:5, 1), 1),
               c(1.5, 2.2, 3, 3.8, 4.5))
})

test_that("ExpSmoothPillars works", {
  m3 <- array(1:12, dim = c(2, 2, 3))
  expect_equal(round(mean(ExpSmoothPillars(m3, 7)), 4), 6.5)
})

test_that("MedianFilterB works", {
  m <- matrix(1:9, nrow = 3)
  m[2:3, 2:3] <- NA
  expect_equal(round(mean(MedianFilterB(m), na.rm = TRUE), 2), NaN)
  expect_equal(round(mean(MedianFilterB(m, na_rm = TRUE), na.rm = TRUE), 3),
               3.812)
  expect_equal(round(mean(MedianFilterB(m, na_count = TRUE), na.rm = TRUE), 2),
               3.42)
  expect_equal(round(mean(SmoothFilterB(m), na.rm = TRUE), 2), NaN)
  expect_equal(round(mean(SmoothFilterB(m, na_rm = TRUE), na.rm = TRUE), 3),
               3.592)
  expect_equal(round(mean(SmoothFilterB(m, na_count = TRUE), na.rm = TRUE), 2),
               3.34)
})

test_that("MeanPillars works", {
  m3 <- array(1:16, dim = c(2, 2, 4))
  expect_equal(MeanPillars(m3), matrix(7:10, nrow = 2))
  expect_equal(MedianPillars(m3), matrix(7:10, nrow = 2))
  expect_equal(round(VarPillars(m3)), matrix(27, nrow = 2, ncol = 2))
})
