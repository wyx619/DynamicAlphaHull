#' DynamicAlphaHull: adaptive alpha-hull range construction
#'
#' `DynamicAlphaHull` constructs range polygons from longitude-latitude
#' occurrence records. Its primary interface is `getDynamicAlphaHull()`, which
#' returns a `terra::SpatVector`.
#'
#' @import data.table
#' @importFrom grDevices chull xy.coords
#' @importFrom purrr imap
#' @importFrom stats complete.cases
#' @importFrom terra aggregate buffer convHull crds erase intersect is.empty is.valid
#' @importFrom terra plot project relate unwrap vect
#' @importClassesFrom terra SpatVector
#' @importMethodsFrom terra plot
#' @importFrom utils data
#'
#' @keywords internal
"_PACKAGE"
