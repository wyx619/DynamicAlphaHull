#' Construct an alpha shape
#'
#' Selects Delaunay edges admitted by a non-negative alpha radius.
#'
#' @param x A two-column coordinate matrix/data frame, or a `delvor` object.
#' @param y Optional y coordinates when `x` is a numeric vector.
#' @param alpha Non-negative alpha-shape radius in the coordinate units.
#'
#' @return An `ashape` object containing selected edges, alpha extremes, and
#'   its Delaunay mesh.
#' @examples
#' coordinates <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0.5, 0.5))
#' shape <- ashape(coordinates, alpha = 0.8)
#' shape$edges
#' @export
ashape <-
function (x, y = NULL, alpha) 
{
    if (alpha < 0) {
        stop("Parameter alpha must be greater or equal to zero")
    }
    if (!inherits(x, "delvor")) {
        delaunay <- delvor(x, y)
    }
    else {
        delaunay <- x
    }
    xy.data <- delaunay$x
    mesh <- delaunay$mesh
    dm1 <- sqrt((mesh[, "x1"] - mesh[, "mx1"])^2 + (mesh[, "y1"] - 
        mesh[, "my1"])^2)
    dm2 <- sqrt((mesh[, "x1"] - mesh[, "mx2"])^2 + (mesh[, "y1"] - 
        mesh[, "my2"])^2)
    dm1[mesh[, "bp1"] == 1] <- Inf
    dm2[mesh[, "bp2"] == 1] <- Inf
    n <- dim(xy.data)[1]
    ind <- 1:n
    ind.on <- chull(xy.data)
    ind.in <- ind[-ind.on]
    n.on <- length(ind.on)
    n.in <- length(ind.in)
    if (n.in > 0) {
        edgeMaximum <- pmax(dm1, dm2)
        maximumByVertex <- data.table(
            vertex = c(mesh[, "ind1"], mesh[, "ind2"]),
            value = c(edgeMaximum, edgeMaximum)
        )[, .(value = max(value)), by = vertex]
        alphaMaximum <- rep(NA_real_, n)
        alphaMaximum[maximumByVertex$vertex] <- maximumByVertex$value
        alpha.ext <- c(ind.on, ind.in[alpha < alphaMaximum[ind.in]])
    }
    else {
        alpha.ext <- ind.on
    }
    n.edges <- dim(mesh)[1]
    is.edge <- numeric()
    ind.is <- 0
    i1 <- match(mesh[, 1], alpha.ext)
    i2 <- match(mesh[, 2], alpha.ext)
    is.edge <- which(i1 & i2)
    aux <- mesh[is.edge, ]
    n.pos <- dim(aux)[1]
    pm.x <- (aux[, "x1"] + aux[, "x2"]) * 0.5
    pm.y <- (aux[, "y1"] + aux[, "y2"]) * 0.5
    dm <- sqrt((aux[, "x1"] - aux[, "x2"])^2 + (aux[, "y1"] - 
        aux[, "y2"])^2) * 0.5
    betw <- rep(NA_real_, n.pos)
    vertical <- aux[, "mx1"] == aux[, "mx2"]
    betweenVertical <- pm.y > pmin(aux[, "my1"], aux[, "my2"]) &
        pm.y < pmax(aux[, "my1"], aux[, "my2"])
    betweenHorizontal <- pm.x > pmin(aux[, "mx1"], aux[, "mx2"]) &
        pm.x < pmax(aux[, "mx1"], aux[, "mx2"])
    betw[(vertical & betweenVertical) | (!vertical & betweenHorizontal)] <- 1
    l.min <- apply(cbind(dm1[is.edge], dm2[is.edge], dm * betw), 
        1, min, na.rm = TRUE)
    l.max <- apply(cbind(dm1[is.edge], dm2[is.edge], dm * betw), 
        1, max, na.rm = TRUE)
    in.ashape <- (l.min <= alpha & alpha <= l.max)
    edges <- matrix(t(aux[in.ashape, ]), byrow = TRUE, ncol = 12)
    colnames(edges) <- colnames(aux)
    alphaShape <- list(edges = edges, length = sum(2 * dm[in.ashape]), 
        alpha = alpha, alpha.extremes = alpha.ext, delaunay = delaunay, 
        x = xy.data)
    class(alphaShape) <- "ashape"
    invisible(alphaShape)
}
