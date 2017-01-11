frcd <-
function (net, seed = seed, maxiter, drp, ...) 
{
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
    ifelse(is.null(rownames(mat)) == TRUE, lbs <- 1:n, lbs <- rownames(mat))
    cmps <- multiplex::comps(netd)
    nds <- matrix(0, nrow = n, ncol = 2)
    rownames(nds) <- lbs
    set.seed(seed)
    for (k in 1:length(cmps$com)) {
        if (is.null(cmps$com[[k]]) == FALSE) {
            com <- which(lbs %in% cmps$com[[k]])
            if (isTRUE(length(cmps$com[[k]]) == 2L) == TRUE) {
                tmp <- popl(2L, seed = seed)/3L
                locx <- tmp[, 1]
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
                  min(locy)))/nrow(mt))
                forcex <- rep(0, nrow(mt))
                forcey <- rep(0, nrow(mt))
                for (niter in 1:maxiter) {
                  for (i in 1:nrow(mt)) {
                    forcevx <- 0
                    forcevy <- 0
                    for (j in 1:nrow(mt)) {
                      if (isTRUE(i != j) == TRUE) {
                        dx <- locx[j] - locx[i]
                        dy <- locy[j] - locy[i]
                        d <- sqrt(dx^2 + dy^2)
                        if (isTRUE(mt[i, j] == 0) == FALSE) {
                          Fd <- (d/K) - K^2/d^2
                        }
                        else if (isTRUE(mt[i, j] == 0) == TRUE) {
                          Fd <- -K^2/d^2
                        }
                        forcevx <- forcevx + Fd * dx
                        forcevy <- forcevy + Fd * dy
                      }
                    }
                    rm(j)
                    forcex[i] <- forcevx
                    forcey[i] <- forcevy
                  }
                  rm(i)
                  TEMP <- 2L/niter
                  for (i in 1:nrow(mt)) {
                    forcemag <- sqrt(forcex[i]^2 + forcey[i]^2)
                    scala <- min(forcemag, TEMP)/forcemag
                    locx[i] <- locx[i] + (forcex[i] * scala)
                    locy[i] <- locy[i] + (forcey[i] * scala)
                  }
                  rm(i)
                }
                rm(niter)
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
    ifelse(isTRUE(sum(nds) == 0) == TRUE, rat <- 1, rat <- (max(nds[, 
        1]) - min(nds[, 1]))/(max(nds[, 2]) - min(nds[, 2])))
    if (isTRUE(length(cmps$isol) > 1) == TRUE) {
        nds[, 1] <- (nds[, 1] - min(nds[, 1]))/(max(nds[, 1]) - 
            min(nds[, 1]))
        ifelse(isTRUE(rat > 0) == TRUE, nds[, 2] <- ((nds[, 2] - 
            min(nds[, 2]))/(max(nds[, 2]) - min(nds[, 2]))) * 
            (1L/rat), nds[, 2] <- ((nds[, 2] - min(nds[, 2]))/(max(nds[, 
            2]) - min(nds[, 2]))) * (rat))
        nds <- as.matrix((nds))
        ndst <- nds[which(nds[, 1] != 0), ]
        tmpi <- popl(length(cmps$isol), seed = seed)/(length(cmps$isol) * 
            2) * K * length(cmps$isol)
        if (is.null(cmps$com) == FALSE) {
            locx <- ((tmpi[, 1]/3L) - (min(ndst[, 1])) - 0)
            ifelse(isTRUE(rat > 0) == TRUE, locy <- ((min(ndst[, 
                2])) - (tmpi[, 2]/3L) - 0), locy <- ((max(ndst[, 
                2])) + (tmpi[, 2]/3L) + 0))
            ndst.chull <- grDevices::chull(ndst)
            ndst.chull <- ndst[ndst.chull, ]
            ifelse(isTRUE(length(which(ndst.chull[, 1] < mean(ndst.chull[, 
                1]))) > length(which(ndst.chull[, 1] > mean(ndst.chull[, 
                1])))) == TRUE, locx <- locx + (K/10L), locx <- locx + 
                ((K/10L) * -1))
            ifelse(isTRUE(length(which(ndst.chull[, 2] < mean(ndst.chull[, 
                2]))) > length(which(ndst.chull[, 2] > mean(ndst.chull[, 
                2])))) == TRUE, locy <- locy - (K/10L), locy <- locy - 
                ((K/10L) * -1))
        }
        else {
            locx <- (tmpi[, 1])
            locy <- (tmpi[, 2])
        }
        nds[which(lbs %in% cmps$isol), ] <- (cbind(locx, locy))
    }
    else if (isTRUE(length(cmps$isol) == 1L) == TRUE) {
        ndst <- nds[which(nds[, 1] != 0), ]
        locx <- max(ndst[, 1]) + (K/10L)
        ifelse(isTRUE(rat < 0) == TRUE, locy <- max(ndst[, 2]) + 
            (K/10L), locy <- min(ndst[, 2]) - (K/10L))
        nds[which(lbs %in% cmps$isol), ] <- (cbind(locx, locy))
    }
    nds[, 1] <- (nds[, 1] - min(nds[, 1]))/(max(nds[, 1]) - min(nds[, 
        1]))
    ifelse(isTRUE(rat > 0) == TRUE, nds[, 2] <- ((nds[, 2] - 
        min(nds[, 2]))/(max(nds[, 2]) - min(nds[, 2]))) * (1L/rat), 
        nds[, 2] <- ((nds[, 2] - min(nds[, 2]))/(max(nds[, 2]) - 
            min(nds[, 2]))) * (rat))
    nds[, 2] <- nds[, 2] * -1
    as.data.frame(nds)
}
