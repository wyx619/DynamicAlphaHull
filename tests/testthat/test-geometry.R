test_that("pointInPolygon handles closed and degenerate polygon boundaries", {
  polygonX <- c(0, 2, 2, 0)
  polygonY <- c(0, 0, 2, 2)

  expect_true(pointInPolygon(1, 1, polygonX, polygonY))
  expect_false(pointInPolygon(3, 1, polygonX, polygonY))
  expect_false(pointInPolygon(1, 1, numeric(), numeric()))
  expect_false(pointInPolygon(1, 1, 0, 0))
})
