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

  points <- vect(as.matrix(coordinates), type = "points", crs = "EPSG:4326")
  pointCoordinates <- crds(points)
  delaunay <- NULL
  alpha <- initialAlpha
  alphaVal <- alpha
  buffered <- FALSE
  problem <- qr(scale(pointCoordinates, scale = FALSE))$rank < 2

  if (!problem) {
    if (verbose) {
      message("\talpha: ", alpha)
    }
    hull <- try(ahull(pointCoordinates, alpha = alpha), silent = TRUE)
    if (!inherits(hull, "try-error")) {
      delaunay <- hull$alphaShape$delaunay
    }

    while (inherits(hull, "try-error") && any(grepl("duplicate points", hull))) {
      if (nrow(points) <= 3) {
        stop("Fewer than 3 usable coordinates remain after duplicate-point removal.")
      }
      pointCoordinates <- crds(points)
      # Drop the closest pair by great-circle distance, matching the original
      # rangeBuilder, which uses sf::st_distance (spherical) rather than the
      # planar degree-space distance. The radius is a constant scale factor and
      # therefore does not affect which pair is closest.
      n <- nrow(pointCoordinates)
      lon <- pointCoordinates[, 1] * pi / 180
      lat <- pointCoordinates[, 2] * pi / 180
      distance <- matrix(Inf, n, n)
      for (i in seq_len(n)) {
        deltaLon <- lon[i] - lon
        haversine <- sin((lat[i] - lat) / 2)^2 +
          cos(lat[i]) * cos(lat) * sin(deltaLon / 2)^2
        distance[i, ] <- 2 * asin(pmin(1, sqrt(haversine)))
      }
      diag(distance) <- Inf
      closest <- which(distance == min(distance), arr.ind = TRUE)
      points <- points[-closest[1, 1], ]
      pointCoordinates <- crds(points)
      delaunay <- NULL
      hull <- try(ahull(pointCoordinates, alpha = alpha), silent = TRUE)
      if (!inherits(hull, "try-error")) {
        delaunay <- hull$alphaShape$delaunay
      }
      if (verbose) {
        message("\tdropping a duplicate point")
      }
    }

    while (inherits(hull, "try-error")) {
      alpha <- alpha + alphaIncrement
      if (alpha > alphaCap) {
        problem <- TRUE
        break
      }
      if (verbose) {
        message("\talpha: ", alpha)
      }
      if (is.null(delaunay)) {
        delaunay <- try(delvor(pointCoordinates), silent = TRUE)
      }
      hull <- if (inherits(delaunay, "try-error")) {
        try(ahull(pointCoordinates, alpha = alpha), silent = TRUE)
      } else {
        try(ahull(delaunay, alpha = alpha), silent = TRUE)
      }
      if (!inherits(hull, "try-error")) {
        delaunay <- hull$alphaShape$delaunay
      }
    }
  }

  if (!problem) {
    hull <- try(ah2terra(hull), silent = TRUE)
    while (inherits(hull, "try-error") || nrow(hull) == 0 ||
           any(is.empty(hull)) || !all(is.valid(hull))) {
      alpha <- alpha + alphaIncrement
      if (alpha > alphaCap) {
        problem <- TRUE
        break
      }
      if (verbose) {
        message("\talpha: ", alpha)
      }
      hull <- try(ah2terra(ahull(delaunay, alpha = alpha)), silent = TRUE)
    }
  }

  if (!problem) {
    pointWithin <- relate(points, hull, "intersects")
    alphaVal <- alpha
    while (nrow(hull) > partCount ||
           sum(pointWithin) / nrow(points) < fraction ||
           !all(is.valid(hull))) {
      alpha <- alpha + alphaIncrement
      if (alpha > alphaCap) {
        problem <- TRUE
        break
      }
      if (verbose) {
        message("\talpha: ", alpha)
      }
      candidate <- try(ahull(delaunay, alpha = alpha), silent = TRUE)
      if (!inherits(candidate, "try-error")) {
        candidate <- try(ah2terra(candidate), silent = TRUE)
      }
      if (!inherits(candidate, "try-error") && nrow(candidate) > 0 &&
          !any(is.empty(candidate)) && all(is.valid(candidate))) {
        hull <- candidate
        hull <- project(hull, "EPSG:8857")
        hull <- buffer(hull, width = buff)
        hull <- project(hull, "EPSG:4326")
        hull <- makeValid(hull)
        buffered <- TRUE
        pointWithin <- relate(points, hull, "intersects")
      }
      alphaVal <- alpha
    }
  }

  if (problem) {
    hull <- convHull(points)
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
