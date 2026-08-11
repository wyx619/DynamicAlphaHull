#' Enumerate Delaunay triangle candidates from an edge table
#'
#' Reproduces the candidate order of `deldir`'s internal triangle enumeration
#' with one adjacency-table pass instead of repeatedly scanning all edges.
#'
#' @param triangulation A `deldir` triangulation.
#'
#' @return An integer matrix with one three-vertex candidate per row. Some
#'   candidates can be non-faces and are filtered by [delaunayTriangles()].
#' @examples
#' triangulation <- deldir::deldir(
#'   x = c(0, 1, 1, 0, 0.5), y = c(0, 0, 1, 1, 0.5), plot = FALSE
#' )
#' delaunayCandidates(triangulation)
#' @export
delaunayCandidates <- function(triangulation) {
  edges <- triangulation$delsgs
  pointCount <- nrow(triangulation$summary)
  adjacency <- data.table(
    vertex = c(edges$ind1, edges$ind2),
    neighbour = c(edges$ind2, edges$ind1)
  )[, .(neighbours = list(neighbour)), by = vertex]
  neighbours <- vector("list", pointCount)
  neighbours[adjacency$vertex] <- adjacency$neighbours
  candidates <- vector("list", 0L)

  for (first in seq_len(pointCount)) {
    secondCandidates <- sort(unique(neighbours[[first]]))
    secondCandidates <- secondCandidates[secondCandidates > first]
    if (length(secondCandidates) == 0L) {
      next
    }
    for (second in secondCandidates) {
      thirdCandidates <- neighbours[[second]]
      thirdCandidates <- thirdCandidates[
        thirdCandidates %in% secondCandidates & thirdCandidates > second
      ]
      if (length(thirdCandidates) > 0L) {
        for (third in thirdCandidates) {
          candidates[[length(candidates) + 1L]] <- c(first, second, third)
        }
      }
    }
  }

  if (length(candidates) == 0L) {
    return(matrix(integer(), ncol = 3L))
  }
  matrix(unlist(candidates, use.names = FALSE), ncol = 3L, byrow = TRUE)
}

#' Extract Delaunay triangle indices without materialising a triangle list
#'
#' Reproduces the candidate order, orientation correction, and `intri` filter
#' used by [deldir::triang.list()], while returning one integer matrix instead
#' of thousands of small data frames.
#'
#' @param triangulation A `deldir` triangulation.
#'
#' @return An integer matrix with one counter-clockwise Delaunay triangle per
#'   row. Indices refer to the original coordinate order.
#' @examples
#' triangulation <- deldir::deldir(
#'   x = c(0, 1, 1, 0, 0.5), y = c(0, 0, 1, 1, 0.5), plot = FALSE
#' )
#' delaunayTriangles(triangulation)
#' @export
delaunayTriangles <- function(triangulation) {
  indices <- triangulation$ind.orig
  candidates <- delaunayCandidates(triangulation)
  if (nrow(candidates) == 0L) {
    return(matrix(integer(), ncol = 3L))
  }

  x <- triangulation$summary[, "x"]
  y <- triangulation$summary[, "y"]
  triangleX <- matrix(x[candidates], nrow(candidates), 3L)
  triangleY <- matrix(y[candidates], nrow(candidates), 3L)
  centredY <- triangleY - min(y)
  deltaX <- cbind(
    triangleX[, 2] - triangleX[, 1],
    triangleX[, 3] - triangleX[, 2],
    triangleX[, 1] - triangleX[, 3]
  )
  clockwise <- rowSums(deltaX * cbind(
    centredY[, 1] + centredY[, 2],
    centredY[, 2] + centredY[, 3],
    centredY[, 3] + centredY[, 1]
  )) > 0
  if (any(clockwise)) {
    triangleX[clockwise, ] <- triangleX[clockwise, c(1, 3, 2), drop = FALSE]
    triangleY[clockwise, ] <- triangleY[clockwise, c(1, 3, 2), drop = FALSE]
    candidates[clockwise, ] <- candidates[clockwise, c(1, 3, 2), drop = FALSE]
  }

  keep <- vapply(seq_len(nrow(candidates)), function(index) {
    as.logical(.Fortran(
      "intri",
      x = as.double(triangleX[index, ]),
      y = as.double(triangleY[index, ]),
      u = as.double(x),
      v = as.double(y),
      n = as.integer(length(x)),
      okay = integer(1),
      PACKAGE = "deldir"
    )[["okay"]])
  }, logical(1))

  matrix(indices[candidates[keep, , drop = FALSE]], ncol = 3L)
}

#' Build the Delaunay mesh used by alpha-shape algorithms
#'
#' @param x,y Coordinates accepted by [grDevices::xy.coords()].
#'
#' @return An object of class `delvor` containing the Delaunay mesh.
#' @examples
#' coordinates <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0.5, 0.5))
#' triangulation <- delvor(coordinates)
#' head(triangulation$mesh)
#' @export
delvor <- function(x, y = NULL) {
  coordinates <- xy.coords(x, y)
  coordinates <- cbind(coordinates$x, coordinates$y)
  if (nrow(coordinates) <= 2 || qr(scale(coordinates, scale = FALSE))$rank < 2) {
    stop("At least three non-collinear points are required.")
  }

  triangulation <- deldir(coordinates[, 1], coordinates[, 2], plot = FALSE)
  triangles <- delaunayTriangles(triangulation)
  if (nrow(triangles) == 0L) {
    stop("Delaunay triangulation did not produce any triangles.")
  }

  circumcentres <- t(vapply(seq_len(nrow(triangles)), function(index) {
    vertices <- coordinates[triangles[index, ], , drop = FALSE]
    denominator <- 2 * (vertices[1, 1] * (vertices[2, 2] - vertices[3, 2]) +
      vertices[2, 1] * (vertices[3, 2] - vertices[1, 2]) +
      vertices[3, 1] * (vertices[1, 2] - vertices[2, 2]))
    if (abs(denominator) < .Machine$double.eps) {
      stop("Delaunay triangulation contains a degenerate triangle.")
    }
    squaredNorm <- rowSums(vertices^2)
    c(
      sum(squaredNorm * c(vertices[2, 2] - vertices[3, 2], vertices[3, 2] - vertices[1, 2], vertices[1, 2] - vertices[2, 2])) / denominator,
      sum(squaredNorm * c(vertices[3, 1] - vertices[2, 1], vertices[1, 1] - vertices[3, 1], vertices[2, 1] - vertices[1, 1])) / denominator
    )
  }, numeric(2)))

  triangleCount <- nrow(triangles)
  edgeFirst <- as.vector(t(cbind(
    triangles[, 1], triangles[, 2], triangles[, 1]
  )))
  edgeSecond <- as.vector(t(cbind(
    triangles[, 2], triangles[, 3], triangles[, 3]
  )))
  edgeLower <- pmin(edgeFirst, edgeSecond)
  edgeUpper <- pmax(edgeFirst, edgeSecond)
  edgeTriangle <- rep(seq_len(triangleCount), each = 3L)
  edgeKey <- paste(edgeLower, edgeUpper, sep = "-")
  orderedRows <- order(edgeKey, seq_along(edgeKey), method = "radix")
  orderedKeys <- edgeKey[orderedRows]
  groupStarts <- c(TRUE, orderedKeys[-1] != orderedKeys[-length(orderedKeys)])
  firstPositions <- which(groupStarts)
  groupCounts <- diff(c(firstPositions, length(orderedRows) + 1L))
  if (any(groupCounts > 2L)) {
    stop("Delaunay mesh contains an edge with more than two adjacent triangles.")
  }
  firstRows <- orderedRows[firstPositions]
  secondRows <- orderedRows[firstPositions + (groupCounts == 2L)]
  edgeCount <- length(firstRows)
  away <- max(diff(range(coordinates[, 1])), diff(range(coordinates[, 2])))
  boundaryHull <- triangulation$summary[
    chull(triangulation$summary[, c("x", "y")]), c("x", "y")
  ]

  mesh <- matrix(0, nrow = edgeCount, ncol = 12L)
  mesh[, 1] <- edgeLower[firstRows]
  mesh[, 2] <- edgeUpper[firstRows]
  mesh[, 3:4] <- coordinates[mesh[, 1], , drop = FALSE]
  mesh[, 5:6] <- coordinates[mesh[, 2], , drop = FALSE]
  mesh[, 7:8] <- circumcentres[edgeTriangle[firstRows], , drop = FALSE]
  interior <- groupCounts == 2L
  if (any(interior)) {
    mesh[interior, 9:10] <- circumcentres[edgeTriangle[secondRows[interior]], , drop = FALSE]
  }
  if (any(!interior)) {
    for (index in which(!interior)) {
      mesh[index, 9:10] <- dummycoor(
        triangulation,
        mesh[index, 3:4],
        mesh[index, 5:6],
        mesh[index, 7:8],
        away,
        boundaryHull
      )
    }
    mesh[!interior, 12] <- 1
  }
  colnames(mesh) <- c(
    "ind1", "ind2", "x1", "y1", "x2", "y2", "mx1", "my1",
    "mx2", "my2", "bp1", "bp2"
  )

  result <- list(mesh = mesh, x = coordinates, triangulation = triangulation)
  class(result) <- "delvor"
  invisible(result)
}
