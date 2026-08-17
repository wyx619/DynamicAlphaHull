#' Construct an adaptive alpha-hull range polygon
#'
#' Starting at `initialAlpha`, the function increases alpha by
#' `alphaIncrement` until the hull contains at least `fraction` of the input
#' occurrences and has no more than `partCount` polygon parts. If no suitable
#' alpha can be found before `alphaCap`, it returns a buffered minimum convex
#' hull instead.
#'
#' Input coordinates are interpreted as WGS 84 longitude and latitude.
#' Buffering is performed in the Equal Earth projection, so `buff` is in
#' metres. Coastline clipping uses the Natural Earth 1:50m land polygon bundled
#' with this package and never downloads or writes data.
#'
#' @param x A data frame or matrix containing longitude and latitude columns.
#' @param fraction Minimum fraction of occurrences included in the polygon, in
#'   `(0, 1]`.
#' @param partCount Maximum number of disjoint polygon parts.
#' @param buff Buffer distance in metres.
#' @param initialAlpha Starting alpha value.
#' @param coordHeaders Names or positions of the longitude and latitude
#'   columns. Ignored when `x` has exactly two columns.
#' @param clipToCoast One of `"no"`, `"terrestrial"`, or `"aquatic"`.
#'   `FALSE` remains a backwards-compatible alias for `"no"`.
#' @param alphaIncrement Positive amount by which alpha increases per iteration.
#' @param verbose Whether to report attempted alpha values.
#' @param alphaCap Maximum alpha to try before the convex-hull fallback.
#'
#' @return A list containing `hull`, a `terra::SpatVector` polygon, and
#'   `alpha`, the selected alpha label. The fallback label is `"alphaMCH"`.
#' @examples
#' points <- data.frame(
#'   Longitude = c(116.0, 116.5, 116.4, 116.1),
#'   Latitude = c(39.0, 39.0, 39.4, 39.3)
#' )
#' range <- getDynamicAlphaHull(points, buff = 1000, clipToCoast = "no")
#' plot(range$hull)
#' points(points$Longitude, points$Latitude, pch = 16)
#' @export
getDynamicAlphaHull <- function(x, fraction = 0.95, partCount = 3,
                                buff = 10000, initialAlpha = 3,
                                coordHeaders = c("Longitude", "Latitude"),
                                clipToCoast = "terrestrial",
                                alphaIncrement = 1, verbose = FALSE,
                                alphaCap = 400) {
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop("x must be a data frame or matrix containing longitude and latitude columns.")
  }
  if (ncol(x) < 2) {
    stop("x must contain at least two columns.")
  }
  if (ncol(x) == 2) {
    coordHeaders <- c(1L, 2L)
  }
  if (length(coordHeaders) != 2 ||
      !all(coordHeaders %in% seq_len(ncol(x)) | coordHeaders %in% colnames(x))) {
    stop("coordHeaders must identify exactly two columns in x.")
  }
  if (!is.numeric(fraction) || length(fraction) != 1 || !is.finite(fraction) || fraction <= 0 || fraction > 1) {
    stop("fraction must be a single number in (0, 1].")
  }
  if (!is.numeric(partCount) || length(partCount) != 1 || !is.finite(partCount) || partCount < 1 || partCount != as.integer(partCount)) {
    stop("partCount must be a positive integer.")
  }
  if (!is.numeric(buff) || length(buff) != 1 || !is.finite(buff) || buff < 0) {
    stop("buff must be a single non-negative distance in metres.")
  }
  if (!is.numeric(initialAlpha) || length(initialAlpha) != 1 || !is.finite(initialAlpha) || initialAlpha < 0) {
    stop("initialAlpha must be a single non-negative number.")
  }
  if (!is.numeric(alphaIncrement) || length(alphaIncrement) != 1 || !is.finite(alphaIncrement) || alphaIncrement <= 0) {
    stop("alphaIncrement must be a single positive number.")
  }
  if (!is.numeric(alphaCap) || length(alphaCap) != 1 || !is.finite(alphaCap) || alphaCap < initialAlpha) {
    stop("alphaCap must be a finite number greater than or equal to initialAlpha.")
  }

  if (isFALSE(clipToCoast)) {
    clipToCoast <- "no"
  }
  clipToCoast <- as.character(clipToCoast)
  clipToCoast <- match.arg(clipToCoast, c("no", "terrestrial", "aquatic"))

  coordinates <- x[, coordHeaders, drop = FALSE] |> as.data.table()
  if (!all(vapply(coordinates, is.numeric, logical(1)))) {
    stop("The coordinate columns must be numeric.")
  }
  setnames(coordinates, c("Longitude", "Latitude"))
  coordinates <- coordinates[
    complete.cases(coordinates) &
      is.finite(coordinates$Longitude) & is.finite(coordinates$Latitude)
  ]
  coordinates <- unique(coordinates)
  if (nrow(coordinates) < 3) {
    stop("This function requires a minimum of 3 unique, finite coordinates.")
  }
  if (any(coordinates$Longitude < -180 | coordinates$Longitude > 180) ||
      any(coordinates$Latitude < -90 | coordinates$Latitude > 90)) {
    stop("Coordinates must be WGS 84 longitude and latitude values.")
  }

  # Shuffle the leading rows while the first three points are collinear,
  # mirroring rangeBuilder's pre-triangulation guard. This consumes RNG via
  # sample() and, under a fixed seed, keeps the point order identical to the
  # reference implementation.
  coordinates <- as.matrix(coordinates)
  allSameX <- all(coordinates[, 1] == coordinates[1, 1])
  allSameY <- all(coordinates[, 2] == coordinates[1, 2])
  while (!allSameX && !allSameY &&
         ((coordinates[1, 1] == coordinates[2, 1] &&
           coordinates[2, 1] == coordinates[3, 1]) ||
          (coordinates[1, 2] == coordinates[2, 2] &&
           coordinates[2, 2] == coordinates[3, 2]))) {
    coordinates <- coordinates[sample(1:nrow(coordinates), size = nrow(coordinates)), , drop = FALSE]
  }

  points <- vect(coordinates, type = "points", crs = "EPSG:4326")
  pointCoordinates <- crds(points)
  alpha <- initialAlpha
  problem <- FALSE

  if (verbose) {
    message("\talpha: ", alpha)
  }

  hull <- try(ahull(pointCoordinates, alpha = alpha), silent = TRUE)

  # Drop near-duplicate points by great-circle distance, mirroring
  # rangeBuilder's sf::st_distance closest-pair search. The full distance
  # matrix is built once (in row chunks to bound peak memory, element-wise
  # identical to the row loop) and then maintained incrementally: removing a
  # point does not change the pairwise distances among the survivors, so each
  # pass only discards the dropped row and column instead of rebuilding the
  # whole matrix.
  if (inherits(hull, "try-error") && any(grepl("duplicate points", hull))) {
    pointCoordinates <- crds(points)
    n <- nrow(pointCoordinates)
    lon <- pointCoordinates[, 1] * pi / 180
    lat <- pointCoordinates[, 2] * pi / 180
    distance <- matrix(Inf, n, n)
    CHUNK <- 512L
    for (start in seq.int(1L, n, by = CHUNK)) {
      end <- min(start + CHUNK - 1L, n)
      idx <- start:end
      deltaLat <- outer(lat[idx], lat, "-") / 2
      deltaLon <- outer(lon[idx], lon, "-") / 2
      haversine <- sin(deltaLat)^2 +
        outer(cos(lat[idx]), cos(lat)) * sin(deltaLon)^2
      d <- 2 * asin(pmin(1, sqrt(haversine)))
      dim(d) <- c(length(idx), n)
      distance[idx, ] <- d
    }
    diag(distance) <- Inf
  }
  while (inherits(hull, "try-error") && any(grepl("duplicate points", hull))) {
    if (nrow(points) <= 3) {
      stop("Fewer than 3 usable coordinates remain after duplicate-point removal.")
    }
    closest <- which(distance == min(distance), arr.ind = TRUE)
    points <- points[-closest[1, 1], ]
    distance <- distance[-closest[1, 1], -closest[1, 1], drop = FALSE]
    pointCoordinates <- crds(points)
    hull <- try(ahull(pointCoordinates, alpha = alpha), silent = TRUE)
    if (verbose) {
      message("\tdropping a duplicate point")
    }
  }

  # Second attempt at the current alpha, matching rangeBuilder.
  hull <- try(ahull(pointCoordinates, alpha = alpha), silent = TRUE)

  while (inherits(hull, "try-error")) {
    if (verbose) {
      message("\talpha: ", alpha)
    }
    alpha <- alpha + alphaIncrement
    hull <- try(ahull(pointCoordinates, alpha = alpha), silent = TRUE)
    if (alpha > alphaCap) {
      problem <- TRUE
      break
    }
  }

  if (!problem) {
    hull <- try(ah2terra(hull), silent = TRUE)

    validityCheck <- function(h) {
      if (!is.null(h) && !inherits(h, "try-error")) {
        if (!all(is.valid(h))) TRUE else FALSE
      } else {
        FALSE
      }
    }

    while (is.null(hull) || inherits(hull, "try-error") || validityCheck(hull)) {
      alpha <- alpha + alphaIncrement
      if (verbose) {
        message("\talpha: ", alpha)
      }
      hull <- try(ah2terra(ahull(pointCoordinates, alpha = alpha)), silent = TRUE)
      if (alpha > alphaCap) {
        problem <- TRUE
        break
      }
    }
  }

  if (!problem) {
    pointWithin <- relate(points, hull, "intersects")
    alphaVal <- alpha
    buffered <- FALSE

    while (nrow(hull) > partCount ||
           sum(pointWithin) / nrow(points) < fraction ||
           !all(is.valid(hull))) {
      alpha <- alpha + alphaIncrement
      if (verbose) {
        message("\talpha: ", alpha)
      }
      hull <- try(ahull(pointCoordinates, alpha = alpha), silent = TRUE)
      while (inherits(hull, "try-error") && alpha <= alphaCap) {
        alpha <- alpha + alphaIncrement
        hull <- try(ahull(pointCoordinates, alpha = alpha), silent = TRUE)
      }
      if (!inherits(hull, "try-error")) {
        hull <- ah2terra(hull)
        hull <- project(hull, "EPSG:8857")
        if (all(is.valid(hull))) {
          hull <- buffer(hull, width = buff)
          hull <- project(hull, "EPSG:4326")
          hull <- makeValid(hull)
          buffered <- TRUE
          pointWithin <- relate(points, hull, "intersects")
        }
      }
      alphaVal <- alpha
      if (alpha > alphaCap) {
        hull <- convHull(points)
        hull <- project(hull, "EPSG:8857")
        hull <- buffer(hull, width = buff)
        hull <- project(hull, "EPSG:4326")
        buffered <- TRUE
        alphaVal <- "MCH"
        break
      }
    }
  } else {
    hull <- convHull(points)
    hull <- project(hull, "EPSG:8857")
    hull <- buffer(hull, width = buff)
    hull <- project(hull, "EPSG:4326")
    buffered <- TRUE
    alphaVal <- "MCH"
  }

  if (!buffered) {
    hull <- project(hull, "EPSG:8857")
    hull <- buffer(hull, width = buff)
    hull <- project(hull, "EPSG:4326")
    hull <- makeValid(hull)
  }

  if (clipToCoast != "no") {
    world <- loadWorldMap()
    hull <- if (clipToCoast == "terrestrial") {
      intersect(hull, world)
    } else {
      erase(hull, world)
    }
  }

  list(hull = hull, alpha = paste0("alpha", alphaVal))
}

worldMapCache <- new.env(parent = emptyenv())

#' Load the bundled offline land polygon
#'
#' Restores the pre-dissolved Natural Earth 1:50m land layer from the package's
#' lazy internal data. The resulting `SpatVector` is cached for the current R
#' session.
#'
#' @return A `terra::SpatVector` containing one dissolved land polygon.
#' @examples
#' world <- loadWorldMap()
#' terra::nrow(world)
#' @export
loadWorldMap <- function() {
  if (exists("world", envir = worldMapCache, inherits = FALSE)) {
    return(get("world", envir = worldMapCache, inherits = FALSE))
  }
  data("ne_50m_land", package = "DynamicAlphaHull", envir = worldMapCache)
  if (!exists("ne_50m_land", envir = worldMapCache, inherits = FALSE)) {
    stop("The bundled Natural Earth land data are missing.")
  }
  world <- unwrap(get("ne_50m_land", envir = worldMapCache, inherits = FALSE))
  rm("ne_50m_land", envir = worldMapCache)
  assign("world", world, envir = worldMapCache)
  world
}

#' Plot a DynamicAlphaHull range polygon
#'
#' Enables the usual base-R `plot()` call for the `terra::SpatVector` returned
#' in `getDynamicAlphaHull()$hull`, without requiring users to attach terra.
#'
#' @param x A `terra::SpatVector`.
#' @param ... Additional arguments passed to `terra::plot()`.
#' @return The plotted `SpatVector`, invisibly.
#' @examples
#' range <- getDynamicAlphaHull(
#'   data.frame(Longitude = c(116.0, 116.5, 116.4, 116.1),
#'              Latitude = c(39.0, 39.0, 39.4, 39.3)),
#'   buff = 1000,
#'   clipToCoast = "no"
#' )
#' plot(range$hull)
#' @exportS3Method graphics::plot
plot.SpatVector <- function(x, ...) {
  plot(x, ...)
  invisible(x)
}
