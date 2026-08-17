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
#' @importFrom terra aggregate buffer convHull crds erase intersect is.empty is.valid makeValid
#' @importFrom terra plot project relate unwrap vect
#' @importClassesFrom terra SpatVector
#' @importMethodsFrom terra plot
#' @importFrom utils data globalVariables
#'
#' @keywords internal
"_PACKAGE"

# data.table non-standard evaluation columns referenced inside `[.data.table`
# expressions in this package.
utils::globalVariables(c(".", "value", "vertex"))
