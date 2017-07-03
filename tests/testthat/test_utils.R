test_that("PillarsListToArr works", {
  arr <- array(1:27, dim = rep(3, 3))
  expect_equal(ListPillars(arr), split(1:27, rep(1:9, 3)), check.names = FALSE)
  expect_equal(PillarsListToArr(ListPillars(arr), dim(arr)), arr)
})

test_that("PillarsDF works", {
  library(magrittr)
  arr <- array(1:27, dim = rep(3, 3))
  expect_equal(PillarsDF(arr),
               ListPillars(arr) %>%
                 set_names(apply(expand.grid(1:3, 1:3), 1,
                                 paste, collapse = "_")) %>%
                 tibble::as_tibble())
})

test_that("CollapseRanges works", {
  set.seed(0)
  ranges <- sort(sample(1:100, 11))
  expect_equal(CollapseRanges(ranges, 6, c(3, 7)), c(6, 27, 37, 56, 61, 85, 96))
  expect_equal(CollapseRanges(ranges, 6, c(3, 7), prefer.low = TRUE),
               c(6, 20, 27, 37, 61, 85, 96))
  expect_equal(CollapseRanges(ranges, 6, c(3, 7), prefer.high = TRUE),
               c(6, 27, 37, 61, 85, 90, 96))
  expect_error(CollapseRanges(ranges, 6, c(3, 7), prefer.high = TRUE,
                              prefer.low = TRUE),
               "One cannot select both prefer.high and prefer.low.")
  expect_error(CollapseRanges(c(1, 2, 0.5), 6, c(3, 7), prefer.high = TRUE),
               "strictly increasing")
  expect_error(CollapseRanges(0:11, 6, 2:9, prefer.high = TRUE),
               "One cannot try to preserve more ranges than one wants overall")
  expect_error(CollapseRanges(0:11, 5, c(3, 6, 9)),
               "The way in which you've chosen n.out and preserve")
  expect_error(CollapseRanges(1:2, 2), "must reduce")
})

test_that("GroupClose works", {
  expect_equal(GroupClose(1:10, 1), list(1:10))
  expect_equal(GroupClose(1:10, 0.5), as.list(1:10))
  expect_equal(GroupClose(c(1, 2, 4, 10, 11, 14, 20, 25, 27), 3),
               list(c(1, 2, 4), c(10, 11, 14), 20, c(25, 27)))
  expect_equal(GroupClose(c(0:4, 6:10)), list(0:4, 6:10))
  expect_error(GroupClose(numeric(0)), "greater than zero")
  expect_equal(GroupClose(1), list(1))
  expect_error(GroupClose(10:1), "strictly increasing")
})

test_that("SpreadSpecific works", {
  expect_equal(SpreadSpecific(c(0, 10), 1, 3), c(1, 5.5, 10))
  expect_equal(SpreadSpecific(c(1, 10), 2, 5, TRUE),
               c(1, 2, 3.42, 5.848, 10), tolerance = 0.001)
  expect_error(SpreadSpecific(c(0, 10), 2, 5, TRUE))
  expect_error(SpreadSpecific(c(0, 10), 11, 5), "members")
})

test_that("Closest works", {
  expect_equal(Closest(pi, 1:10), 3)
})

test_that("FixLUTError works", {
  has.lut.error <- abind::abind(matrix(0, nrow = 2, ncol = 2),
                                matrix(1:4, nrow = 2), matrix(0, nrow = 2, ncol = 2),
                                along = 3)
  expect_equal(FixLUTError(has.lut.error, 2), matrix(1:4, nrow = 2))
  has.lut.error3d <- abind::abind(has.lut.error, has.lut.error,
                                  along = 4)
  expect_equal(FixLUTError(has.lut.error3d, 3),
               abind::abind(matrix(1:4, nrow = 2, ncol = 2),
                            matrix(1:4, nrow = 2, ncol = 2),
                            along = 3),
               check.attributes = FALSE)
  has.no.error <- abind::abind(matrix(runif(4), nrow = 2, ncol = 2),
                               matrix(runif(4), nrow = 2),
                               matrix(runif(4), nrow = 2, ncol = 2),
                               along = 3)
  expect_equal(FixLUTError(has.no.error, 3), has.no.error)
  expect_error(FixLUTError(has.lut.error3d, 2), "dimensionality by 1")
  expect_error(FixLUTError(abind::abind(has.lut.error, has.lut.error,
                                        along = 3), 2),
                           "the third dimension has value 3, however")
  will.lut.error <- abind::abind(matrix(3, nrow = 2, ncol = 2),
                                 matrix(1:4, nrow = 2),
                                 matrix(0, nrow = 2, ncol = 2),
                                 along = 3)
  expect_error(FixLUTError(will.lut.error, 2),
               "all but one of these color channels is zero")
})

test_that("ScaleBreaksFun works", {
  expect_equal(nandb:::ScaleBreaksFun(3:4, FALSE)(c(0, 10)), c(0, 3, 4, 7, 10))
  expect_equal(nandb:::ScaleBreaksFun(3:4, FALSE)(c(2, 10)), c(3, 4, 6, 8, 10))
  expect_equal(round(nandb:::ScaleBreaksFun(3:4, TRUE)(c(2, 10)), 1),
               c(2, 3, 4, 6.3, 10))
  expect_equal(round(nandb:::ScaleBreaksFun(log = TRUE)(c(2, 10)), 1),
               c(2, 3, 4.5, 6.7, 10))
  expect_equal(nandb:::ScaleBreaksFun(), ggplot2::waiver())

})

test_that("Mat2ColList works", {
  expect_equal(nandb:::Mat2ColList(matrix(1:6, nrow = 3)), list(1:3, 4:6))
})

test_that("ListChannels works", {
  arr4d <- array(1:(3 ^ 4), dim = rep(3, 4))
  expect_equal(nandb:::ListChannels(arr4d, 3), list(arr4d[, , 1, ],
                                                 arr4d[, , 2, ],
                                                 arr4d[, , 3, ]))
})

test_that("Slices errors correctly", {
  arr3d <- array(1:27, dim = rep(3, 3))
  expect_error(nandb:::Slices(99, arr3d))
  arr4d <- array(1:(3 ^ 4), dim = rep(3, 4))
  expect_error(nandb:::Slices(1, arr4d))
})

test_that("ApplyOnPillars works in special case", {
  arr3d <- array(1:27, dim = rep(3, 3))
  expect_equal(ApplyOnPillars(arr3d, identity), arr3d)
})

test_that("WhichInterval errors correctly", {
  expect_error(nandb:::WhichInterval(1:3, t(c(5, 5))))
})

test_that("BrightnessVec works", {
  img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
  expect_equal(Brightness(img), ApplyOnPillars(img, nandb:::BrightnessVec),
               check.attributes = FALSE)
  expect_equal(nandb:::BrightnessVec(rep(0, 10)), 0)
})

test_that("ListChannels errors correctly", {
  library(magrittr)
  arr <- array(seq_len(2^3), dim = rep(2, 3)) %>%
    abind::abind(., ., ., ., ., along = 4) %>%
    aperm(c(1, 2, 4, 3))
  expect_error(nandb:::ListChannels(arr, 2), "but you have indicated")
})

test_that("ChannelList2Arr errors correctly", {
  library(magrittr)
  expect_error(nandb:::ChannelList2Arr(list(matrix(1:4, nrow = 2),
                                            matrix(1:6, nrow = 2))),
               "dimensions.*same")
  expect_error(nandb:::ChannelList2Arr(list(array(seq_len(2^5),
                                                  dim = rep(2, 5)))[c(1, 1)]),
               "dimensional")

})

