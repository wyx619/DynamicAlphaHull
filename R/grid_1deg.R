#' 1-degree global grid for biogeographic range analysis
#'
#' A wrapped terra SpatVector containing a global 1-degree grid system used for
#' discretizing species geographic ranges. The grid is stored in wrapped format
#' for efficient serialization and must be unwrapped with \code{terra::unwrap()}
#' before use.
#'
#' @format A wrapped terra SpatVector (use \code{terra::unwrap(grid_1deg)} to restore).
#' The unwrapped object is a polygon layer with 19,196 features and 5 attributes:
#' \describe{
#'   \item{grid_id}{Integer. Unique identifier for each grid cell (1 to 19,196)}
#'   \item{latitude}{Numeric. Latitude of grid cell centroid}
#'   \item{longitude}{Numeric. Longitude of grid cell centroid}
#'   \item{纬度段}{Character. Latitudinal band label (Chinese)}
#'   \item{A99}{Numeric. Additional attribute from source data}
#' }
#'
#' @details
#' \strong{Grid Specification:}
#' \itemize{
#'   \item \strong{Resolution}: 1 degree × 1 degree cells
#'   \item \strong{CRS}: WGS 84 (EPSG:4326)
#'   \item \strong{Coverage}: Global, from 56°S to 84°N latitude
#'   \item \strong{Extent}: -180° to 180° longitude, -55.98° to 84.02° latitude
#'   \item \strong{Cell count}: 19,196 polygons
#' }
#'
#' \strong{Usage Pattern:}
#'
#' The grid is stored in wrapped format to enable efficient R package distribution.
#' Always unwrap before use:
#'
#' \code{grid <- terra::unwrap(grid_1deg)}
#'
#' For range-building workflows, use \code{read_range_grid()} and
#' \code{new_range_context()} which handle unwrapping automatically.
#'
#' \strong{Grid Cell Indexing:}
#'
#' Each cell has a unique \code{grid_id} (1–19,196). Range outputs reference
#' cells by this ID, enabling compact storage of species distributions as
#' integer vectors rather than full geometries.
#'
#' \strong{Antimeridian Handling:}
#'
#' The grid properly handles the ±180° meridian. Species ranges that cross the
#' antimeridian (e.g., Potentilla x gorodkovii: -179.76° to +174.33°) are
#' correctly represented across the appropriate grid cells.
#'
#' @source
#' Derived from the reference 1-degree grid system used in rangeBuilder and
#' alphahull validation workflows. Original file: \code{alphahull/R/1d/1.shp}.
#'
#' @seealso
#' \code{\link{read_range_grid}}, \code{\link{new_range_context}},
#' \code{\link{build_family_ranges}}
#'
#' @examples
#' # Load and unwrap the grid
#' data(grid_1deg)
#' grid <- terra::unwrap(grid_1deg)
#'
#' # Inspect grid properties
#' terra::nrow(grid)
#' terra::crs(grid)
#' terra::ext(grid)
#'
#' # Grid cell IDs
#' range(grid$grid_id)
#'
#' # Centroid coordinates
#' head(grid[, c("grid_id", "latitude", "longitude")])
#'
#' \dontrun{
#' # Typical workflow: use helper functions instead of unwrapping manually
#' grid <- read_range_grid(grid_1deg)  # accepts wrapped or file path
#' context <- new_range_context(grid)
#'
#' # Build ranges
#' data(rosales_example)
#' result <- build_family_ranges(rosales_example, context, workers = 4L)
#'
#' # Grid distribution shows which cells each species occupies
#' head(result$grid_distribution)
#' }
"grid_1deg"
