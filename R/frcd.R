frcd <-
function (net, seed = seed, maxiter, drp, scl, mov, ...) 
{
    if (isTRUE(is.data.frame(net) == TRUE) == TRUE) {
        net <- as.matrix(net)
    }
    net <- replace(net, net == Inf, 0L)
    n <- dim(net)[1]
    ifelse(missing(maxiter) == TRUE, maxiter <- (60L + n), NA)
    if (missing(drp) == FALSE && is.numeric(drp) == TRUE) {
        netdrp <- replace(net, net <= drp, 0)
        netd <- multiplex::dichot(netdrp, c = 1L)
    }
    else {
        netd <- multiplex::dichot(net, c = 1L)
        netdrp <- net
    }
    mat <- multiplex::mnplx(netd, directed = FALSE)
    nds <- matrix(0, nrow = n, ncol = 2)
    if (is.null(rownames(mat)) == TRUE) {
        lbs <- seq_len(n)
        rownames(nds) <- lbs
    }
    else {
        lbs <- rownames(mat)
        rownames(nds) <- make.unique(lbs)
    }
    cmps <- multiplex::comps(netd)
    set.seed(seed)
    for (k in 1:length(cmps$com)) {
        if (is.null(cmps$com[[k]]) == FALSE) {
            com <- which(lbs %in% cmps$com[[k]])
            if (isTRUE(length(cmps$com[[k]]) == 2L) == TRUE) {
                tmp <- rbind(stats::rcauchy(2L, location = -20, 
                  scale = 2), stats::rcauchy(2L, location = 20, 
                  scale = 2))/(n * (length(cmps$com[[k]])^4))
                locx <- tmp[, 1] + (stats::rcauchy(1L, location = 1, 
                  scale = 2) * 1L)
                locy <- tmp[, 2]
                K <- 1
            }
            else if (isTRUE(length(cmps$com[[k]]) > 2L) == TRUE) {
                mt <- mat[com, com]
                ndc <- data.frame(X = round(stats::runif(dim(mt)[1]) * 
                  2L, 5), Y = round(stats::runif(dim(mt)[2]) * 
                  2L, 5))
                locx <- ndc[, 1]
                locy <- ndc[, 2]
                K <- 0.75 * sqrt(((max(locx) - min(locx)) * (max(locy) - 
                  min(locy)))/dim(mt)[1])
                forcex <- rep(0, dim(mt)[1])
                forcey <- rep(0, dim(mt)[1])
                for (niter in seq_len(maxiter)) {
                  for (i in seq_len(dim(mt)[1])) {
                    forcevx <- 0
                    forcevy <- 0
                    for (j in which(!(seq_len(dim(mt)[1]) %in% 
                      i))) {
                      dx <- locx[j] - locx[i]
                      dy <- locy[j] - locy[i]
                      d <- sqrt(dx^2 + dy^2)
                      ifelse(isTRUE(mt[i, j] == 0) == FALSE, 
                        Fd <- (d/K) - K^2/d^2, Fd <- (-K^2/d^2))
                      forcevx <- forcevx + Fd * dx
                      forcevy <- forcevy + Fd * dy
                    }
                    rm(j)
                    forcex[i] <- forcevx
                    forcey[i] <- forcevy
                  }
                  rm(i)
                  TEMP <- 2L/niter
                  for (l in seq_len(dim(mt)[1])) {
                    forcemag <- sqrt(forcex[l]^2 + forcey[l]^2)
                    scala <- min(forcemag, TEMP)/forcemag
                    locx[l] <- locx[l] + (forcex[l] * scala)
                    locy[l] <- locy[l] + (forcey[l] * scala)
                  }
                  rm(l)
                }
                rm(niter)
            }
            else {
                NA
            }
            if (all(nds == 0) == FALSE && isTRUE(length(cmps$com) > 
                1) == TRUE) {
                ndsc <- nds[which(nds[, 1] != 0 | nds[, 2] != 
                  0), ]
                if (isTRUE(length(which(ndsc[, 1] < mean(ndsc[, 
                  1]))) > length(which(ndsc[, 1] > mean(ndsc[, 
                  1])))) == TRUE) {
                  dirx <- +1L
                  locx <- max(ndsc[, 1]) + (max(locx) - locx)
                }
                else {
                  dirx <- -1L
                  locx <- min(ndsc[, 1]) - (max(locx) - locx)
                }
                if (isTRUE(length(which(ndsc[, 2] < mean(ndsc[, 
                  2]))) > length(which(ndsc[, 2] > mean(ndsc[, 
                  2])))) == TRUE) {
                  diry <- +1L
                  locy <- max(ndsc[, 2]) + (max(locy) - locy)
                }
                else {
                  diry <- -1L
                  locy <- min(ndsc[, 2]) - (min(locy) - locy)
                }
                if (isTRUE(any(duplicated(rbind(ndsc, cbind(locx, 
                  locy))))) == TRUE) {
                  locx <- locx + (dirx * 0.25)
                  locy <- locy + (diry * 0.25)
                }
                else {
                  NA
                }
            }
            nds[com, ] <- (cbind(locx, locy))
        }
        else {
            K <- 2
        }
    }
    rm(k)
    ifelse(isTRUE(sum(nds) == 0) == TRUE, rat <- 1L, rat <- (max(nds[, 
        1]) - min(nds[, 1]))/(max(nds[, 2]) - min(nds[, 2])))
    ni <- length(cmps$isol)
    if (isTRUE(ni > 1) == TRUE) {
        nds[, 1] <- (nds[, 1] - min(nds[, 1]))/(max(nds[, 1]) - 
            min(nds[, 1]))
        ifelse(isTRUE(rat > 0) == TRUE, nds[, 2] <- ((nds[, 2] - 
            min(nds[, 2]))/(max(nds[, 2]) - min(nds[, 2]))) * 
            (1L/rat), nds[, 2] <- ((nds[, 2] - min(nds[, 2]))/(max(nds[, 
            2]) - min(nds[, 2]))) * (rat))
        nds <- as.matrix((nds))
        ndst <- nds[which(nds[, 1] != 0), ]
        tmpi <- popl(ni, seed = seed) * (ni/10L)
        if (is.null(cmps$com) == FALSE) {
            lcx <- ((tmpi[, 1]) - (min(ndst[, 1])))
            ifelse(isTRUE(rat > 0) == TRUE, lcy <- ((min(ndst[, 
                2])) - abs(tmpi[, 2]) - 0), lcy <- ((max(ndst[, 
                2])) + abs(tmpi[, 2]) + 0))
            ndst.chull <- grDevices::chull(ndst)
            ndst.chull <- ndst[ndst.chull, ]
            ifelse(isTRUE(length(which(ndst.chull[, 1] < mean(ndst.chull[, 
                1]))) > length(which(ndst.chull[, 1] > mean(ndst.chull[, 
                1])))) == TRUE, lcx <- lcx + abs(K * 0.5), lcx <- lcx * 
                -1)
            ifelse(isTRUE(length(which(ndst.chull[, 2] < mean(ndst.chull[, 
                2]))) > length(which(ndst.chull[, 2] > mean(ndst.chull[, 
                2])))) == TRUE, lcy <- lcy - abs(K * 0.5), lcy <- lcy * 
                -1)
        }
        else {
            lcx <- (tmpi[, 1])
            lcy <- (tmpi[, 2])
        }
        nds[which(lbs %in% cmps$isol), ] <- cbind(lcx, lcy)
    }
    else if (isTRUE(ni == 1L) == TRUE) {
        ndst <- nds[which(nds[, 1] != 0), ]
        lcx <- max(ndst[, 1]) + (K/(n - 1L))
        ifelse(isTRUE(rat < 0) == TRUE, lcy <- max(ndst[, 2]) + 
            (K/(n - 1L)), lcy <- min(ndst[, 2]) - (K/(n - 1L)))
        nds[which(lbs %in% cmps$isol), ] <- (cbind(lcx, lcy)) * 
            K
    }
    nds[, 2] <- nds[, 2] * -1
    if (missing(scl) == FALSE) {
        nds[, 1] <- nds[, 1] * scl[1]
        nds[, 2] <- nds[, 2] * scl[2]
    }
    if (missing(mov) == FALSE) {
        nds[, 1] <- nds[, 1] + mov[1]
        nds[, 2] <- nds[, 2] + mov[2]
    }
    as.data.frame(nds)
}
