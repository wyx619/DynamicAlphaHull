#' Convert an alpha hull to a terra polygon vector
#'
#' @param x An object of class `ahull`.
#' @param increment Number of points used to approximate a full circle.
#' @param rnd Decimal places retained while joining arc endpoints.
#' @param crs Output coordinate reference system.
#' @param tol Coordinate tolerance for closing rings.
#'
#' @return A `terra::SpatVector` containing polygon features.
#' @examples
#' coordinates <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0.5, 0.5))
#' hull <- ahull(coordinates, alpha = 0.8)
#' polygon <- ah2terra(hull)
#' plot(polygon)
#' @export
ah2terra <- function(x, increment = 360, rnd = 10, crs = "EPSG:4326", tol = 1e-4) {
  if (!inherits(x, "ahull")) {
    stop("x must be an ahull object.")
  }

  arcs <- as.data.frame(x$arcs)
  arcs <- arcs[arcs$r > 0, , drop = FALSE]
  if (nrow(arcs) == 0) {
    return(vect(matrix(numeric(), ncol = 2), type = "polygons", crs = crs))
  }

  arcs$flip <- FALSE
  unused <- seq_len(nrow(arcs))
  rings <- list()

  while (length(unused) > 0) {
    first <- unused[1]
    chain <- first
    unused <- unused[-1]
    start <- arcs$end1[first]
    endpoint <- arcs$end2[first]

    while (endpoint != start && length(unused) > 0) {
      forward <- unused[arcs$end1[unused] == endpoint]
      backward <- unused[arcs$end2[unused] == endpoint]
      if (length(forward) > 0) {
        nextArc <- forward[1]
      } else if (length(backward) > 0) {
        nextArc <- backward[1]
        arcs$flip[nextArc] <- TRUE
      } else {
        break
      }
      chain <- c(chain, nextArc)
      unused <- unused[unused != nextArc]
      endpoint <- if (arcs$flip[nextArc]) arcs$end1[nextArc] else arcs$end2[nextArc]
    }

    if (endpoint != start) {
      next
    }

    coordinates <- NULL
    for (index in chain) {
      row <- arcs[index, ]
      angles <- anglesArc(c(row$v.x, row$v.y), row$theta)
      if (row$flip) {
        angles <- rev(angles)
      }
      arcAngles <- seq(angles[1], angles[2], length.out = 2 + round(increment * row$theta / 2))
      arcCoordinates <- cbind(
        round(row$c1 + row$r * cos(arcAngles), rnd),
        round(row$c2 + row$r * sin(arcAngles), rnd)
      )
      coordinates <- if (is.null(coordinates)) arcCoordinates else rbind(coordinates, arcCoordinates[-1, , drop = FALSE])
    }

    if (nrow(coordinates) >= 4 &&
        max(abs(coordinates[1, ] - coordinates[nrow(coordinates), ])) <= tol) {
      coordinates[nrow(coordinates), ] <- coordinates[1, ]
      rings[[length(rings) + 1]] <- coordinates
    }
  }

  if (length(rings) == 0) {
    return(vect(matrix(numeric(), ncol = 2), type = "polygons", crs = crs))
  }

  geometry <- rbindlist(imap(rings, function(ring, index) {
    as.data.table(cbind(as.integer(index), 1L, ring))
  }))
  vect(as.matrix(geometry), type = "polygons", crs = crs)
}
