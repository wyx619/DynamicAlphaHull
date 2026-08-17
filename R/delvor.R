#' Build the Delaunay mesh used by alpha-shape algorithms
#'
#' Reproduces `alphahull::delvor`, which triangulates with `interp::tri.mesh`.
#' Using the same triangulation as the reference implementation keeps the edge
#' order (and therefore the alpha-hull arc order) identical, so that downstream
#' polygon construction matches the original rangeBuilder output exactly.
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
  X <- xy.coords(x, y)
  coordinates <- cbind(X$x, X$y)
  if (nrow(coordinates) <= 2) {
    stop("At least three non-collinear points are required.")
  }

  tri.obj <- interp::tri.mesh(X)
  tri <- interp::triangles(tri.obj)

  # Circumcentres, vectorised from interp::circum (src/circum.cpp). The centre
  # depends only on the double-precision side lengths and barycentric weights;
  # the float fields circum also computes (radius, aspect ratio) are unused
  # here. The operations mirror circum.cpp exactly for bit-identical results.
  x1 <- coordinates[tri[, 1], 1]
  y1 <- coordinates[tri[, 1], 2]
  x2 <- coordinates[tri[, 2], 1]
  y2 <- coordinates[tri[, 2], 2]
  x3 <- coordinates[tri[, 3], 1]
  y3 <- coordinates[tri[, 3], 2]
  sideA <- sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))
  sideB <- sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2))
  sideC <- sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3))
  weightA <- sideA * sideA * (-sideA * sideA + sideB * sideB + sideC * sideC)
  weightB <- sideB * sideB * (sideA * sideA - sideB * sideB + sideC * sideC)
  weightC <- sideC * sideC * (sideA * sideA + sideB * sideB - sideC * sideC)
  weightSum <- weightA + weightB + weightC
  weightA <- weightA / weightSum
  weightB <- weightB / weightSum
  weightC <- weightC / weightSum
  circenter <- cbind(
    circumx = weightB * x1 + weightC * x2 + weightA * x3,
    circumy = weightB * y1 + weightC * y2 + weightA * y3
  )
  tri.info <- cbind(tri, circenter)

  n.tri <- nrow(tri.info)
  if (n.tri == 1) {
    aux1 <- cbind(matrix(tri.info[, c("arc1", "node2", "node3")], ncol = 3, nrow = 1),
                  1:n.tri, tri.info[, "tr1"])
    aux2 <- cbind(matrix(tri.info[, c("arc2", "node1", "node3")], ncol = 3, nrow = 1),
                  1:n.tri, tri.info[, "tr2"])
    aux3 <- cbind(matrix(tri.info[, c("arc3", "node1", "node2")], ncol = 3, nrow = 1),
                  1:n.tri, tri.info[, "tr3"])
  } else {
    aux1 <- cbind(tri.info[, c("arc1", "node2", "node3")], 1:n.tri, tri.info[, "tr1"])
    aux2 <- cbind(tri.info[, c("arc2", "node1", "node3")], 1:n.tri, tri.info[, "tr2"])
    aux3 <- cbind(tri.info[, c("arc3", "node1", "node2")], 1:n.tri, tri.info[, "tr3"])
  }
  aux <- rbind(aux1, aux2, aux3)
  repeated <- duplicated(aux[, 1])
  aux <- aux[!repeated, ]
  colnames(aux) <- c("arc", "ind1", "ind2", "indm1", "indm2")

  bp1 <- (aux[, "indm1"] == 0)
  bp2 <- (aux[, "indm2"] == 0)
  is.dummy <- which(bp2)
  # Preallocate the exact number of exterior circumcentres (one per boundary
  # edge) so the fill loop never re-copies the growing matrix. rbind() names
  # each appended row after its argument ("dum"); reproduce that row-name
  # attribute so the final mesh is unchanged from the reference output.
  circumcentres <- matrix(nrow = nrow(tri.info) + length(is.dummy), ncol = 2)
  circumcentres[seq_len(nrow(tri.info)), ] <- tri.info[, c("circumx", "circumy")]
  rownames(circumcentres) <- c(rep("", nrow(tri.info)), rep("dum", length(is.dummy)))
  away <- max(diff(range(coordinates[, 1])), diff(range(coordinates[, 2])))
  for (i in is.dummy) {
    n.tri <- n.tri + 1
    circumcentres[n.tri, ] <- dummycoor(
      tri.obj,
      coordinates[aux[i, "ind1"], ],
      coordinates[aux[i, "ind2"], ],
      tri.info[aux[i, "indm1"], c("circumx", "circumy")],
      away
    )
    aux[i, "indm2"] <- n.tri
  }

  mesh <- cbind(
    aux[, c("ind1", "ind2")],
    coordinates[aux[, "ind1"], ],
    coordinates[aux[, "ind2"], ],
    circumcentres[aux[, "indm1"], ],
    circumcentres[aux[, "indm2"], ],
    bp1, bp2
  )
  colnames(mesh) <- c(
    "ind1", "ind2", "x1", "y1", "x2", "y2", "mx1", "my1",
    "mx2", "my2", "bp1", "bp2"
  )

  result <- list(mesh = mesh, x = coordinates, tri.obj = tri.obj)
  class(result) <- "delvor"
  invisible(result)
}
