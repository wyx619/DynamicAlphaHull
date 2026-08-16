#' Compute the angular limits of an arc
#'
#' @param vector Unit direction vector from an arc centre.
#' @param theta Half-angle of the arc in radians.
#'
#' @return A numeric vector of start and end angles in radians.
#' @examples
#' anglesArc(c(1, 0), pi / 4)
#' @export
anglesArc <- function(vector, theta) {
  angle <- if (vector[2] >= 0) acos(vector[1]) else 2 * pi - acos(vector[1])
  c(angle - theta, angle + theta)
}

#' Rotate a two-dimensional vector
#'
#' @param vector Numeric vector of length two.
#' @param theta Rotation angle in radians.
#'
#' @return A rotated numeric vector of length two.
#' @examples
#' rotation(c(1, 0), pi / 2)
#' @export
rotation <- function(vector, theta) {
  c(
    cos(theta) * vector[1] + sin(theta) * vector[2],
    -sin(theta) * vector[1] + cos(theta) * vector[2]
  )
}

#' Calculate the intersection geometry of two circles
#'
#' @param c11,c12 Centre coordinates of the first circle.
#' @param r1 Radius of the first circle.
#' @param c21,c22 Centre coordinates of the second circle.
#' @param r2 Radius of the second circle.
#'
#' @return A list with the number of intersections, two direction vectors, and
#'   the corresponding half-angles.
#' @examples
#' inter(0, 0, 1, 1, 0, 1)
#' @export
inter <- function(c11, c12, r1, c21, c22, r2) {
  distance <- sqrt((c21 - c11)^2 + (c22 - c12)^2)
  vector1 <- c(0, 0)
  theta1 <- 0
  theta2 <- 0
  if (distance == 0 && r1 == r2) {
    cuts <- Inf
  } else if (distance > r1 + r2 || distance < abs(r1 - r2)) {
    cuts <- 0
  } else if (distance == r1 + r2 || distance == abs(r1 - r2)) {
    cuts <- 1
  } else {
    cuts <- 2
    distance1 <- (distance^2 - r2^2 + r1^2) / (2 * distance)
    distance2 <- distance - distance1
    vector1 <- c((c21 - c11) / distance, (c22 - c12) / distance)
    theta1 <- acos(distance1 / r1)
    theta2 <- acos(distance2 / r2)
  }
  list(n.cut = cuts, v1 = vector1, theta1 = theta1, v2 = -vector1, theta2 = theta2)
}

#' Calculate alpha-hull arc length
#'
#' @param arcs Matrix or data frame containing `theta` and `r` columns.
#'
#' @return Total arc length.
#' @examples
#' arcs <- data.frame(theta = c(pi / 4, pi / 2), r = c(2, 1))
#' lengthahull(arcs)
#' @export
lengthahull <- function(arcs) {
  sum(2 * arcs[, "theta"] * arcs[, "r"])
}

#' Test whether a point lies in a polygon
#'
#' @param x,y Coordinates of the query point.
#' @param polygonX,polygonY Coordinates of polygon vertices.
#'
#' @return `TRUE` when the point is inside the polygon, otherwise `FALSE`.
#' @examples
#' pointInPolygon(0.5, 0.5, c(0, 1, 1, 0), c(0, 0, 1, 1))
#' @export
pointInPolygon <- function(x, y, polygonX, polygonY) {
  vertexCount <- length(polygonX)
  if (vertexCount == 0) return(FALSE)
  nextVertex <- if (vertexCount == 1) 1L else c(2:vertexCount, 1L)
  intersects <- (polygonY > y) != (polygonY[nextVertex] > y) &
    x < (polygonX[nextVertex] - polygonX) * (y - polygonY) /
      (polygonY[nextVertex] - polygonY) + polygonX
  sum(intersects) %% 2L == 1L
}

#' Construct an exterior circumcentre for a boundary edge
#'
#' @param triangulation A `deldir` triangulation returned by [deldir::deldir()].
#' @param firstPoint,secondPoint Numeric coordinate vectors defining an edge.
#' @param centre Circumcentre of the adjacent Delaunay triangle.
#' @param away Distance used to place the exterior point.
#' @param hull Optional precomputed convex-hull coordinate matrix.
#'
#' @return Numeric coordinate vector for the exterior circumcentre.
#' @examples
#' triangulation <- deldir::deldir(c(0, 1, 0), c(0, 0, 1), plot = FALSE)
#' dummycoor(
#'   triangulation,
#'   firstPoint = c(0, 0), secondPoint = c(1, 0),
#'   centre = c(0.5, 0.5), away = 1
#' )
#' @export
dummycoor <- function(triangulation, firstPoint, secondPoint, centre, away,
                      hull = NULL) {
  normal <- c(secondPoint[2] - firstPoint[2], firstPoint[1] - secondPoint[1])
  norm <- sum(normal^2)
  if (norm == 0) return(centre)
  normal <- normal / norm
  midpoint <- (firstPoint + secondPoint) / 2
  testPoint <- midpoint + normal * 1e-5
  if (is.null(hull)) {
    hull <- triangulation$summary[
      chull(triangulation$summary[, c("x", "y")]), c("x", "y")
    ]
  }
  if (pointInPolygon(testPoint[1], testPoint[2], hull[, 1], hull[, 2])) {
    centre - away * normal
  } else {
    centre + away * normal
  }
}
