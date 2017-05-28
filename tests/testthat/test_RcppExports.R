test_that("ExpSmooth works", {
  expect_equal(round(ExpSmooth(1:5, 1), 2),
               c(0.97, 1.99, 3, 4.01, 5.03))
  expect_equal(round(ExpSmooth(c(0, 0, 3, 0, 0), 10), 2),
               c(0.48, 0.49, 0.52, 0.49, 0.48))
  expect_equal(ExpSmooth(1:2, 10, extended = T), rep(NA_real_, 2))
})

test_that("ExpSmoothPillars works", {
  m3 <- array(1:12, dim = c(2, 2, 3))
  expect_equal(round(mean(ExpSmoothPillars(m3, 7)), 4), 6.5)
  expect_equal(ExpSmoothPillars(m3, NA), m3)
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

test_that("MedReflectExtend works", {
  expect_equal(nandb:::MedReflectExtend(1:10), (-8):19)
  expect_equal(nandb:::MedReflectExtend(8), 8)
  expect_equal(nandb:::ReflectIndexMed(1, 0, ""), NA_real_)
})

test_that("Smooth works", {
  expect_equal(nandb:::Smooth(3), 3)
})

test_that("ExpSmoothRows works with extended = FALSE", {
  rows <- rbind(c(1, 3, 4),
                c(6, 7, 8))
  library(magrittr)
  smoothed.rounded <- ExpSmoothRows(rows, 22) %>% round(2)
  expect_equal(smoothed.rounded, rbind(c(2.55, 2.67, 2.78),
                                       c(6.92, 7, 7.08)))
})

