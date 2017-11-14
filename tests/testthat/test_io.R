context("I/O")

test_that("text-image-io works", {
  setwd(tempdir())
  mm <- matrix(1:60, nrow = 4)
  mmm <- abind::abind(mm, mm, along = 3)
  mmmm <- abind::abind(mmm, mmm, along = 4)
  write_txt_img(mmmm, "mmmm")
  expect_equal(as.vector(mmmm),
               unlist(lapply(dir(pattern = "txt$"), read_txt_img)))
  suppressWarnings(file.remove(list.files()))
})
