#' Dissolved Natural Earth 1:50m land polygon
#'
#' A pre-dissolved Natural Earth 1:50m land polygon stored as lazy internal
#' package data for offline coastline clipping. It is a packed
#' `terra::SpatVector` in WGS 84.
#'
#' Restore it with `data("ne_50m_land")` and `terra::unwrap()` the result, or
#' simply call `loadWorldMap()`, which performs both steps and caches the object
#' for the current R session.
#'
#' @format A packed `SpatVector` with a single dissolved land multipolygon in
#'   WGS 84.
#' @source Natural Earth 1:50m land,
#'   \url{https://www.naturalearthdata.com/}, dissolved to a single polygon.
#' @name ne_50m_land
#' @docType data
#' @keywords datasets
NULL
