#' Construct an alpha hull
#'
#' Computes the alpha hull of planar coordinates from their alpha shape and
#' complement. The resulting object can be converted to a `SpatVector` with
#' `ah2terra()`.
#'
#' @param x A two-column coordinate matrix/data frame, or a `delvor` object.
#' @param y Optional y coordinates when `x` is a numeric vector.
#' @param alpha Non-negative alpha-hull radius in the coordinate units.
#'
#' @return An `ahull` object containing hull arcs, retained coordinates,
#'   complement information, and alpha.
#' @examples
#' coordinates <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0.5, 0.5))
#' hull <- ahull(coordinates, alpha = 0.8)
#' hull$arcs
#' @export
ahull <-
function (x, y = NULL, alpha)
{
    alphaShape <- ashape(x, y, alpha)
    compl <- complement(alphaShape$delaunay, alpha = alpha)
    pm.x <- (compl[, "x1"] + compl[, "x2"]) * 0.5
    pm.y <- (compl[, "y1"] + compl[, "y2"]) * 0.5
    dm <- sqrt((compl[, "x1"] - compl[, "x2"])^2 + (compl[, "y1"] -
        compl[, "y2"])^2) * 0.5
    ashape.edges <- matrix(alphaShape$edges[, c("ind1", "ind2")],
        ncol = 2, byrow = FALSE)
    noforget <- alphaShape$alpha.extremes
    ind2 <- integer()
    j <- 0
    nshape <- length(ashape.edges) * 0.5
    arcs <- matrix(0, nrow = nshape, ncol = 6)
    indp <- matrix(0, nrow = nshape, ncol = 2)
    cutp <- alphaShape$x
    if (nshape > 0) {
        # Look up complement rows by their (ind1, ind2) endpoint pair in one
        # pass instead of scanning every complement row for each edge. split()
        # returns each group's row indices in ascending order, exactly as
        # which() would, so the per-edge result is unchanged.
        keyScale <- max(compl[, "ind1"], compl[, "ind2"]) + 1
        compGroups <- split(seq_len(nrow(compl)), compl[, "ind1"] * keyScale + compl[, "ind2"])
        edgeKey <- as.character(ashape.edges[, 1] * keyScale + ashape.edges[, 2])
        for (i in 1:nshape) {
            ind <- compGroups[[edgeKey[i]]]
            if (length(ind) > 0) {
                if (!((1 <= sum((compl[ind, "ind"] == 1))) &
                  (sum((compl[ind, "ind"] == 1)) < length(ind)))) {
                  posIdx <- compl[ind, "r"] > 0
                  if (any(posIdx)) {
                    which <- which(compl[ind, "r"] == min(compl[ind[posIdx], "r"]))
                  j <- j + 1
                  arcs[j, ] <- c(compl[ind[which], 1], compl[ind[which],
                    2], compl[ind[which], 3], compl[ind[which],
                    "v.x"], compl[ind[which], "v.y"], compl[ind[which],
                    "theta"])
                  vaux <- compl[ind[which], c("x1", "y1")] -
                    cbind(pm.x, pm.y)[ind[which], ]
                  theta.aux <- rotation(compl[ind[which], c("v.x",
                    "v.y")], compl[ind[which], "theta"])
                  a2 <- sum(vaux * theta.aux)
                  if (a2 > 0) {
                    indp[j, ] <- compl[ind[which], c("ind1",
                      "ind2")]
                  }
                  else {
                    indp[j, ] <- compl[ind[which], c("ind2",
                      "ind1")]
                  }
                }
                  }
                ind2 <- c(ind2, ind)
            }
        }
    }
    arcs.old <- arcs[arcs[, 3] > 0, ]
    colnames(arcs.old) <- c("c1", "c2", "r", "v.x", "v.y", "theta")
    arcs <- arcs.old
    indp <- indp[indp[, 1] != 0 & indp[, 2] != 0, ]
    n.arc <- dim(arcs)[1]
    watch <- 1
    j <- 1
    if (n.arc > 0) {
        # Working buffers with geometric growth. arcs and indp stay in lockstep
        # (row i of indp describes arc i); cutp keeps its reference-style rbind()
        # appends, which carry deparsed row names that the output must retain.
        # The buffers are trimmed back before downstream use, so the result is
        # identical to the reference implementation which reallocated on every
        # append.
        arcsColnames <- colnames(arcs)
        cap <- max(4L, 2L * n.arc)
        arcsBuffer <- matrix(0, nrow = cap, ncol = 6)
        colnames(arcsBuffer) <- arcsColnames
        arcsBuffer[seq_len(n.arc), ] <- arcs
        indpBuffer <- matrix(0, nrow = cap, ncol = 2)
        indpBuffer[seq_len(n.arc), ] <- indp
        arcs <- arcsBuffer
        indp <- indpBuffer

        ensureCapacity <- function() {
            if (n.arc == cap) {
                cap <<- 2L * cap
                arcs <<- rbind(arcs, matrix(0, nrow = cap - n.arc, ncol = 6))
                indp <<- rbind(indp, matrix(0, nrow = cap - n.arc, ncol = 2))
            }
            invisible(NULL)
        }

        while (watch <= n.arc) {
            initialArcCount <- n.arc
            watchC1 <- arcs[watch, 1]
            watchC2 <- arcs[watch, 2]
            watchR <- arcs[watch, 3]
            deltaX <- watchC1 - arcs[seq_len(initialArcCount), 1]
            deltaY <- watchC2 - arcs[seq_len(initialArcCount), 2]
            centreDistances <- sqrt(deltaX^2 + deltaY^2)
            canIntersect <- centreDistances > abs(watchR - arcs[seq_len(initialArcCount), 3]) &
              centreDistances < watchR + arcs[seq_len(initialArcCount), 3]
            canIntersect[watch] <- FALSE
            while (j <= n.arc) {
                if (j != watch) {
                  if (j <= initialArcCount) {
                    pairCanIntersect <- canIntersect[j]
                  } else {
                    arcjC1 <- arcs[j, 1]
                    arcjC2 <- arcs[j, 2]
                    arcjR <- arcs[j, 3]
                    centreDistance <- sqrt((watchC1 - arcjC1)^2 + (watchC2 - arcjC2)^2)
                    pairCanIntersect <- centreDistance > abs(watchR - arcjR) &&
                      centreDistance < watchR + arcjR
                  }
                  if (pairCanIntersect) {
                    intersection <- inter(watchC1, watchC2, watchR, arcs[j, 1], arcs[j, 2], arcs[j, 3])
                    if (intersection$n.cut == 2) {
                    v.arc <- c(arcs[watch, "v.x"], arcs[watch,
                      "v.y"])
                    if (v.arc[2] >= 0) {
                      ang.OX <- acos(pmax(pmin(v.arc[1], 1), -1))
                    }
                    else {
                      ang.OX <- 2 * pi - acos(pmax(pmin(v.arc[1], 1), -1))
                    }
                    v.arc.rot <- rotation(v.arc, ang.OX)
                    v.int <- intersection$v1
                    v.int.rot <- rotation(v.int, ang.OX)
                    if (v.int.rot[2] >= 0) {
                      ang.v.int.rot.OX <- acos(pmax(pmin(v.int.rot[1], 1), -1))
                      angles <- c(-arcs[watch, "theta"], arcs[watch,
                        "theta"], ang.v.int.rot.OX - intersection$theta1,
                        ang.v.int.rot.OX + intersection$theta1)
                      names(angles) <- c("theta1", "theta2",
                        "beta1", "beta2")
                      order <- names(sort.int(angles))
                      theta1 <- ang.OX - arcs[watch, "theta"]
                      theta2 <- ang.OX + arcs[watch, "theta"]
                      beta1 <- ang.v.int.rot.OX + ang.OX - intersection$theta1
                      beta2 <- ang.v.int.rot.OX + ang.OX + intersection$theta1
                    }
                    else {
                      ang.v.int.rot.OX <- acos(pmax(pmin(v.int.rot[1], 1), -1))
                      angles <- c(-arcs[watch, "theta"], arcs[watch,
                        "theta"], -ang.v.int.rot.OX - intersection$theta1,
                        -ang.v.int.rot.OX + intersection$theta1)
                      names(angles) <- c("theta1", "theta2",
                        "beta1", "beta2")
                      order <- names(sort.int(angles))
                      theta1 <- ang.OX - arcs[watch, "theta"]
                      theta2 <- ang.OX + arcs[watch, "theta"]
                      beta1 <- -ang.v.int.rot.OX + ang.OX - intersection$theta1
                      beta2 <- -ang.v.int.rot.OX + ang.OX + intersection$theta1
                    }
                    if (sum(match(indp[watch, ], indp[j, ], nomatch = 0)) >
                      0) {
                      coinc <- indp[j, sum(match(indp[watch,
                        ], indp[j, ], nomatch = 0))]
                      if (indp[watch, 1] == indp[j, 2]) {
                        if (all(order == c("beta1", "beta2",
                          "theta1", "theta2"))) {
                          case <- 1
                        }
                        else if (all(order == c("theta1", "beta1",
                          "beta2", "theta2"))) {
                          case <- 2
                        }
                        else if (all(order == c("beta1", "theta1",
                          "beta2", "theta2"))) {
                          ang.control <- (angles["theta1"] -
                            angles["beta1"])/2
                          if (abs(ang.control) < 1e-05) {
                            case <- 2
                          }
                          else {
                            case <- 1
                          }
                        }
                        if (case == 2) {
                          ang.middle2 <- (angles["theta2"] -
                            angles["beta2"])/2
                          v.new2 <- rotation(c(1, 0), -arcs[watch,
                            "theta"] + ang.middle2 - ang.OX)
                          cutp <- rbind(cutp, arcs[watch, 1:2] +
                            arcs[watch, 3] * rotation(v.new2,
                              ang.middle2))
                          inn <- dim(cutp)[1]
                          arcs[watch, 4:6] <- c(v.new2, ang.middle2)
                          indp[watch, 1] <- inn
                          pmaux <- (cutp[inn, ] + cutp[indp[j,
                            1], ]) * 0.5
                          dmaux <- pmaux - cutp[indp[j, 1], ]
                          ndmaux <- sqrt(sum(dmaux^2))
                          vaux <- pmaux - arcs[j, 1:2]
                          nvaux <- sqrt(sum(vaux^2))
                          th <- atan(ndmaux/nvaux)
                          arcs[j, 4:6] <- c(vaux/nvaux, th)
                          indp[j, 2] <- inn
                        }
                      }
                      else if (indp[watch, 2] == indp[j, 1]) {
                        if (all(order == c("theta1", "theta2",
                          "beta1", "beta2"))) {
                          case <- 1
                        }
                        else if (all(order == c("theta1", "beta1",
                          "beta2", "theta2"))) {
                          case <- 2
                        }
                        else if (all(order == c("theta1", "beta1",
                          "theta2", "beta2"))) {
                          ang.control <- (angles["theta2"] -
                            angles["beta2"])/2
                          if (abs(ang.control) < 1e-05) {
                            case <- 2
                          }
                          else {
                            case <- 1
                          }
                        }
                        if (case == 2) {
                          ang.middle <- (angles["beta1"] - angles["theta1"])/2
                          v.new <- rotation(c(1, 0), arcs[watch,
                            "theta"] - ang.middle - ang.OX)
                          cutp <- rbind(cutp, arcs[watch, 1:2] +
                            arcs[watch, 3] * rotation(v.new,
                              -ang.middle))
                          inn <- dim(cutp)[1]
                          arcs[watch, 4:6] <- c(v.new, ang.middle)
                          indp[watch, 2] <- inn
                          pmaux <- (cutp[inn, ] + cutp[indp[j,
                            2], ]) * 0.5
                          dmaux <- pmaux - cutp[indp[j, 2], ]
                          ndmaux <- sqrt(sum(dmaux^2))
                          vaux <- pmaux - arcs[j, 1:2]
                          nvaux <- sqrt(sum(vaux^2))
                          th <- atan(ndmaux/nvaux)
                          arcs[j, 4:6] <- c(vaux/nvaux, th)
                          indp[j, 1] <- inn
                        }
                      }
                    }
                    else if (all(order == c("theta1", "beta1",
                      "beta2", "theta2"))) {
                      v.arcj <- c(arcs[j, "v.x"], arcs[j, "v.y"])
                      if (v.arcj[2] >= 0) {
                        ang.OXj <- acos(v.arcj[1])
                      }
                      else {
                        ang.OXj <- 2 * pi - acos(v.arcj[1])
                      }
                      v.arc.rotj <- rotation(v.arcj, ang.OXj)
                      v.intj <- intersection$v2
                      v.int.rotj <- rotation(v.intj, ang.OXj)
                      if (v.int.rotj[2] >= 0) {
                        ang.v.int.rot.OXj <- acos(pmax(pmin(v.int.rotj[1], 1), -1))
                        anglesj <- c(-arcs[j, "theta"], arcs[j,
                          "theta"], ang.v.int.rot.OXj - intersection$theta2,
                          ang.v.int.rot.OXj + intersection$theta2)
                        names(anglesj) <- c("theta1", "theta2",
                          "beta1", "beta2")
                        orderj <- names(sort.int(anglesj))
                        theta1j <- ang.OXj - arcs[j, "theta"]
                        theta2j <- ang.OXj + arcs[j, "theta"]
                        beta1j <- ang.v.int.rot.OXj + ang.OXj -
                          intersection$theta2
                        beta2j <- ang.v.int.rot.OXj + ang.OXj +
                          intersection$theta2
                      }
                      else {
                        ang.v.int.rot.OXj <- acos(pmax(pmin(v.int.rotj[1], 1), -1))
                        anglesj <- c(-arcs[j, "theta"], arcs[j,
                          "theta"], -ang.v.int.rot.OXj - intersection$theta2,
                          -ang.v.int.rot.OXj + intersection$theta2)
                        names(anglesj) <- c("theta1", "theta2",
                          "beta1", "beta2")
                        orderj <- names(sort.int(anglesj))
                        theta1j <- ang.OXj - arcs[j, "theta"]
                        theta2j <- ang.OXj + arcs[j, "theta"]
                        beta1j <- -ang.v.int.rot.OXj + ang.OXj -
                          intersection$theta2
                        beta2j <- -ang.v.int.rot.OXj + ang.OXj +
                          intersection$theta2
                      }
                      if (all(orderj == c("theta1", "beta1",
                        "beta2", "theta2"))) {
                        ang.middle <- (angles["beta1"] - angles["theta1"])/2
                        v.new <- rotation(c(1, 0), arcs[watch,
                          "theta"] - ang.middle - ang.OX)
                        ang.middle2 <- (angles["theta2"] - angles["beta2"])/2
                        v.new2 <- rotation(c(1, 0), -arcs[watch,
                          "theta"] + ang.middle2 - ang.OX)
                        ensureCapacity()
                        n.arc <- n.arc + 1
                        arcs[n.arc, ] <- c(arcs[watch, 1], arcs[watch,
                          2], arcs[watch, 3], v.new2[1], v.new2[2],
                          ang.middle2)
                        arcs[watch, ] <- c(arcs[watch, 1], arcs[watch,
                          2], arcs[watch, 3], v.new[1], v.new[2],
                          ang.middle)
                        np1 <- arcs[watch, 1:2] + arcs[watch,
                          3] * rotation(v.new, -ang.middle)
                        np2 <- arcs[watch, 1:2] + arcs[watch,
                          3] * rotation(v.new2, ang.middle2)
                        indold <- indp[watch, 2]
                        inn1 <- dim(cutp)[1] + 1
                        inn2 <- dim(cutp)[1] + 2
                        indp[watch, 2] <- inn1
                        indp[n.arc, ] <- c(inn2, indold)
                        cutp <- rbind(cutp, np1)
                        cutp <- rbind(cutp, np2)
                        indold <- indp[j, 1]
                        indp[j, 1] <- inn1
                        ensureCapacity()
                        indp[n.arc + 1, ] <- c(indold, inn2)
                        pmaux <- (cutp[inn1, ] + cutp[indp[j,
                          2], ]) * 0.5
                        dmaux <- pmaux - cutp[indp[j, 2], ]
                        ndmaux <- sqrt(sum(dmaux^2))
                        vaux <- pmaux - arcs[j, 1:2]
                        nvaux <- sqrt(sum(vaux^2))
                        th <- atan(ndmaux/nvaux)
                        arcs[j, 4:6] <- c(vaux/nvaux, th)
                        pmaux <- (cutp[inn2, ] + cutp[indold,
                          ]) * 0.5
                        dmaux <- pmaux - cutp[indold, ]
                        ndmaux <- sqrt(sum(dmaux^2))
                        vaux <- pmaux - arcs[j, 1:2]
                        nvaux <- sqrt(sum(vaux^2))
                        th <- atan(ndmaux/nvaux)
                        ensureCapacity()
                        n.arc <- n.arc + 1
                        arcs[n.arc, ] <- c(arcs[j, 1:3], vaux/nvaux,
                          th)
                      }
                    }
                  }
                  }
                }
                case <- 0
                j <- j + 1
            }
            watch <- watch + 1
            j <- 1
        }
        arcs <- arcs[seq_len(n.arc), , drop = FALSE]
        indp <- indp[seq_len(n.arc), , drop = FALSE]
        ord.old <- 1:dim(indp)[1]
        ord.new <- numeric()
        while (length(ord.new) < length(ord.old)) {
            if (length(ord.new) == 0) {
                ord.new <- 1
            }
            else {
                ord.new <- c(ord.new, ord.old[-ord.new][1])
            }
            coinc <- match(indp[ord.new[length(ord.new)], 2],
                indp[-ord.new, 1])
            while (!is.na(coinc)) {
                ord.new <- c(ord.new, ord.old[-ord.new][coinc])
                coinc <- match(indp[ord.new[length(ord.new)],
                  2], indp[-ord.new, 1])
            }
        }
        indp <- indp[ord.new, ]
        ahull.arcs <- cbind(arcs[ord.new, ], indp)
        colnames(ahull.arcs) <- c("c1", "c2", "r", "v.x", "v.y",
            "theta", "end1", "end2")
        lengthah <- lengthahull(arcs)
        addp <- noforget[is.na(match(noforget, indp))]
        num <- length(addp)
        if (num > 0) {
            mat.noforget <- cbind(matrix(alphaShape$x[addp, 1:2],
                ncol = 2, byrow = FALSE), rep(0, num), rep(0,
                num), rep(0, num), rep(0, num), addp, addp)
            ahull.arcs <- rbind(ahull.arcs, mat.noforget)
        }
    }
    else {
        num <- length(noforget)
        ahull.arcs <- cbind(alphaShape$x[noforget, ], rep(0,
            num), rep(0, num), rep(0, num), rep(0, num), noforget,
            noforget)
        colnames(ahull.arcs) <- c("c1", "c2", "r", "v.x", "v.y",
            "theta", "end1", "end2")
        lengthah <- 0
    }
    alphaHull <- list(arcs = ahull.arcs, xahull = cutp, length = lengthah,
        complement = compl, alpha = alpha, alphaShape = alphaShape)
    class(alphaHull) <- "ahull"
    invisible(alphaHull)
}
