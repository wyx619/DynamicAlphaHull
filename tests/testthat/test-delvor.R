test_that("optimized Delaunay extraction preserves deldir triangle order", {
  set.seed(20260811)
  points <- cbind(runif(30), runif(30))
  triangulation <- deldir::deldir(points[, 1], points[, 2], plot = FALSE)

  expectedCandidates <- deldir:::prelimtlist(triangulation)
  expect_identical(delaunayCandidates(triangulation), expectedCandidates)

  expectedTriangles <- as.matrix(data.table::rbindlist(purrr::map(
    deldir::triang.list(triangulation),
    function(triangle) as.list(triangle[["ptNum"]])
  )))
  expect_identical(
    unname(delaunayTriangles(triangulation)),
    unname(expectedTriangles)
  )
})
