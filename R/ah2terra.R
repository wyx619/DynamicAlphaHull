#' Convert an alpha hull to a terra polygon vector
#'
#' Reproduces the arc ordering and ring construction of the reference
#' implementation (`rangeBuilder::ah2sf`, itself derived from `alphahull::ah2sp`)
#' so that the resulting polygons match the original rangeBuilder output exactly.
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
  arcs <- cbind(arcs, flip = rep(FALSE, nrow(arcs)))

  # Reorder arcs so that consecutive arcs share an endpoint. Ported verbatim from
  # rangeBuilder::ah2sf.
  if (nrow(arcs) > 0) {
    k <- 1
    repeat {
      if (is.na(arcs[k + 1, "end1"])) {
        break
      }
      if (arcs[k, "end2"] == arcs[k + 1, "end1"]) {
        k <- k + 1
      } else if (arcs[k, "end2"] != arcs[k + 1, "end1"] &&
                 !arcs[k, "end2"] %in% arcs$end1[k + 1:nrow(arcs)] &&
                 !arcs[k, "end2"] %in% arcs$end2[k + 1:nrow(arcs)]) {
        k <- k + 1
      } else if (arcs[k, "end2"] != arcs[k + 1, "end1"] &&
                 arcs[k, "end2"] %in% arcs$end1[k + 1:nrow(arcs)] &&
                 !arcs[k, "end2"] %in% arcs$end2[k + 1:nrow(arcs)]) {
        m <- which(arcs$end1[k + 1:nrow(arcs)] == arcs[k, "end2"]) + k
        arcs <- rbind(arcs[1:k, ], arcs[m, ], arcs[setdiff((k + 1):nrow(arcs), m), ])
      } else if (arcs[k, "end2"] != arcs[k + 1, "end1"] &&
                 !arcs[k, "end2"] %in% arcs$end1[k + 1:nrow(arcs)] &&
                 arcs[k, "end2"] %in% arcs$end2[k + 1:nrow(arcs)]) {
        m <- which(arcs$end2[k + 1:nrow(arcs)] == arcs[k, "end2"]) + k
        tmp1 <- arcs[m, "end1"]
        tmp2 <- arcs[m, "end2"]
        arcs[m, "end1"] <- tmp2
        arcs[m, "end2"] <- tmp1
        arcs[m, "flip"] <- TRUE
        arcs <- rbind(arcs[1:k, ], arcs[m, ], arcs[setdiff((k + 1):nrow(arcs), m), ])
      } else {
        k <- k + 1
      }
    }
  }

  arcs <- arcs[arcs$r > 0, , drop = FALSE]

  if (nrow(arcs) == 0) {
    return(vect(matrix(numeric(), ncol = 2), type = "polygons", crs = crs))
  }

  # Build rings from the ordered arcs, ported from rangeBuilder::ah2sf.
  rings <- list()
  prevx <- NULL
  prevy <- NULL

  for (i in seq_len(nrow(arcs))) {
    rowi <- arcs[i, ]
    v <- c(rowi$v.x, rowi$v.y)
    theta <- rowi$theta
    r <- rowi$r
    cc <- c(rowi$c1, rowi$c2)
    ipoints <- 2 + round(increment * (rowi$theta / 2), 0)
    angles <- anglesArc(v, theta)
    if (isTRUE(rowi$flip)) {
      angles <- rev(angles)
    }
    seqang <- seq(angles[1], angles[2], length = ipoints)
    xx <- round(cc[1] + r * cos(seqang), rnd)
    yy <- round(cc[2] + r * sin(seqang), rnd)

    if (is.null(prevx)) {
      prevx <- xx
      prevy <- yy
    } else if ((xx[1] == round(prevx[length(prevx)], rnd) ||
                abs(xx[1] - prevx[length(prevx)]) < tol) &&
               (yy[1] == round(prevy[length(prevy)], rnd) ||
                abs(yy[1] - prevy[length(prevy)]) < tol)) {
      if (i == nrow(arcs)) {
        prevx <- append(prevx, xx[2:ipoints])
        prevy <- append(prevy, yy[2:ipoints])
        prevx[length(prevx)] <- prevx[1]
        prevy[length(prevy)] <- prevy[1]
        rings[[length(rings) + 1]] <- cbind(prevx, prevy)
      } else {
        prevx <- append(prevx, xx[2:ipoints])
        prevy <- append(prevy, yy[2:ipoints])
      }
    } else {
      prevx[length(prevx)] <- prevx[1]
      prevy[length(prevy)] <- prevy[1]
      rings[[length(rings) + 1]] <- cbind(prevx, prevy)
      prevx <- NULL
      prevy <- NULL
    }
  }

  # Drop rings with fewer than four points (matches ah2sf's badLines filter).
  rings <- Filter(function(m) nrow(m) >= 4, rings)

  if (length(rings) == 0) {
    return(vect(matrix(numeric(), ncol = 2), type = "polygons", crs = crs))
  }

  geometry <- rbindlist(imap(rings, function(ring, index) {
    as.data.table(cbind(as.integer(index), 1L, ring))
  }))
  vect(as.matrix(geometry), type = "polygons", crs = crs)
}
