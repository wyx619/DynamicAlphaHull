test_that("getDynamicAlphaHull returns an un-clipped polygon", {
  points <- data.frame(
    Longitude = c(116.0, 116.5, 116.4, 116.1),
    Latitude = c(39.0, 39.0, 39.4, 39.3)
  )

  result <- getDynamicAlphaHull(points, buff = 1000, initialAlpha = 2, clipToCoast = "no")

  expect_s4_class(result$hull, "SpatVector")
  expect_false(any(terra::is.empty(result$hull)))
  expect_identical(result$alpha, "alpha2")
})

test_that("coastline clipping works from the bundled offline data", {
  points <- data.frame(
    Longitude = c(116.0, 116.5, 116.4, 116.1),
    Latitude = c(39.0, 39.0, 39.4, 39.3)
  )

  result <- getDynamicAlphaHull(points, buff = 1000, initialAlpha = 2, clipToCoast = "terrestrial")

  expect_s4_class(result$hull, "SpatVector")
  expect_false(any(terra::is.empty(result$hull)))
})

test_that("invalid input and collinear points are handled deterministically", {
  expect_error(
    getDynamicAlphaHull(data.frame(Longitude = c(116, 116), Latitude = c(39, 40))),
    "minimum of 3"
  )

  collinear <- data.frame(
    Longitude = c(116, 116, 116, 116),
    Latitude = c(39, 39.1, 39.2, 39.3)
  )
  result <- getDynamicAlphaHull(collinear, buff = 1000, clipToCoast = "no")
  expect_identical(result$alpha, "alphaMCH")
  expect_s4_class(result$hull, "SpatVector")
})
