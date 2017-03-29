test_that("MeanStackThresh works", {
  img <- ReadImageData(system.file("extdata", "50.tif",
    package = "nandb"))
  img_thresh_mask <- MeanStackThresh(img, "Otsu")
  expect_equal(round(mean(img_thresh_mask, na.rm = TRUE), 4), 25.2892)
  img_thresh_mask <- MeanStackThresh(img, "Triangle")
  expect_equal(round(mean(img_thresh_mask, na.rm = TRUE), 4), 23.5941)
})
