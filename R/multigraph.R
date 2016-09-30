multigraph <-
function (net, layout = c("circ", "stress", "rand"), directed = TRUE, 
    collRecip = TRUE, coord = NULL, tcol = 1, bwd = 1, clu = NULL, 
    cex = NULL, tcex = NULL, showLbs = TRUE, showAtts = TRUE, 
    att = NULL, lbat = "1", seed = NULL, maxiter = 100, alpha = c(1, 
        1, 1), main = NULL, cex.main, bg, mar, pos, asp, lwd, 
    pch, lty, ecol, vcol, vcol0, hds, vedist, rot, mirrorX, mirrorY, 
    col, ...) 
{
    if (isTRUE(is.array(net)) == FALSE) 
        stop("\"net\" should be an array.")
    ifelse(isTRUE(dim(net)[3] == 1) == TRUE, net <- net[, , 1], 
        NA)
    n <- dim(net)[1]
    if (is.array(att) == TRUE) {
        if (isTRUE(n != dim(att)[1]) == TRUE) {
            warning("Dimensions in \"net\" and \"att\" differ. No attributes are shown.")
            showAtts <- FALSE
        }
    }
    ifelse(isTRUE(directed == TRUE) == TRUE, NA, directed <- FALSE)
    ifelse(is.null(dimnames(net)[[1]]) == TRUE, lbs <- as.character(1:n), 
        lbs <- dimnames(net)[[1]])
    if (isTRUE(length(alpha) < 2) == TRUE) {
        alfa <- 1
        alpha <- rep(alpha, 3)
    }
    else {
        alfa <- alpha[2]
    }
    if (isTRUE(length(alpha) < 3) == TRUE) 
        alpha <- append(alpha, 0.1)
    ifelse(missing(asp) == TRUE, asp <- 1, NA)
    ifelse(missing(lwd) == TRUE, lwd <- 1, NA)
    ifelse(missing(pch) == TRUE, pch <- 21, NA)
    if (missing(pos) == TRUE) {
        if (match.arg(layout) == "circ") {
            cdp <- data.frame(X = sin(2L * pi * ((0:(n - 1L))/n)), 
                Y = cos(2L * pi * ((0:(n - 1L))/n)))
            cdp[which(cdp[, 2] == -1), 1] <- 0
            cdp <- replace(cdp, cdp < 0, -1)
            cdp <- replace(cdp, cdp > 0, 1)
            pos <- vector()
            for (i in 1:nrow(cdp)) {
                if (all(as.vector(cdp[i, ]) == c(-1, -1)) == 
                  TRUE | all(as.vector(cdp[i, ]) == c(-1, 0)) == 
                  TRUE | all(as.vector(cdp[i, ]) == c(-1, 1)) == 
                  TRUE) {
                  pos <- append(pos, 2)
                }
                else if (all(as.vector(cdp[i, ]) == c(1, -1)) == 
                  TRUE | all(as.vector(cdp[i, ]) == c(1, 0)) == 
                  TRUE | all(as.vector(cdp[i, ]) == c(1, 1)) == 
                  TRUE) {
                  pos <- append(pos, 4)
                }
                else if (all(as.vector(cdp[i, ]) == c(0, 1)) == 
                  TRUE) {
                  pos <- append(pos, 3)
                }
                else {
                  pos <- append(pos, 1)
                }
            }
            rm(i)
            if (missing(mirrorX) == FALSE && isTRUE(mirrorX == 
                TRUE) == TRUE) {
                pos[which(pos == 2)] <- 0
                pos[which(pos == 4)] <- 2
                pos[which(pos == 0)] <- 4
            }
            else if (missing(mirrorY) == FALSE && isTRUE(mirrorY == 
                TRUE) == TRUE) {
                pos[which(pos == 1)] <- 0
                pos[which(pos == 3)] <- 1
                pos[which(pos == 0)] <- 3
            }
            else {
                NA
            }
        }
        else {
            pos <- 4
        }
    }
    else {
        if (isTRUE(pos < 0L) == TRUE | isTRUE(pos > 4L) == TRUE) 
            stop("Invalid \"pos\" value.")
    }
    ifelse(missing(bg) == TRUE, bg <- graphics::par()$bg, NA)
    ifelse(missing(mar) == TRUE, mar <- graphics::par()$mar, 
        NA)
    ifelse(missing(cex.main) == TRUE, cex.main <- graphics::par()$cex.main, 
        NA)
    ifelse(isTRUE(bwd > 1L) == TRUE, bwd <- 1L, NA)
    ifelse(isTRUE(bwd <= 0L) == TRUE, bwd <- 0L, NA)
    if (!(missing(hds))) {
        ifelse(isTRUE(hds == 0L) == TRUE, hds <- 0.01, NA)
    }
    else {
        hds <- 1L
    }
    ifelse(missing(vedist) == TRUE, vedist <- 0, NA)
    ifelse(missing(rot) == TRUE, NA, rot <- rot * -1)
    ifelse(isTRUE(vedist > 1L) == TRUE, vedist <- 1L, NA)
    net <- multiplex::dichot(net, c = 1L)
    if (isTRUE(directed == FALSE) == TRUE) {
        if ((is.na(dim(net)[3]) == TRUE | isTRUE(dim(net)[3] == 
            1) == TRUE)) {
            nt <- net + t(net)
            rcp <- multiplex::dichot(nt, c = 2L)
            if (isTRUE(collRecip == TRUE) == TRUE) {
                rcp[lower.tri(rcp)] <- 0L
            }
        }
        else {
            nt <- array(0L, dim = c(n, n, dim(net)[3]))
            dimnames(nt)[[1]] <- dimnames(nt)[[2]] <- dimnames(net)[[1]]
            dimnames(nt)[[3]] <- dimnames(net)[[3]]
            for (i in 1:dim(net)[3]) nt[, , i] <- multiplex::dichot(net[, 
                , i], c = 1L) + t(multiplex::dichot(net[, , i], 
                c = 1L))
            rcp <- multiplex::dichot(nt, c = 2L)
            if (isTRUE(collRecip == TRUE) == TRUE) {
                for (i in 1:dim(net)[3]) rcp[, , i][lower.tri(rcp[, 
                  , i])] <- 0L
            }
        }
        if (isTRUE(collRecip == TRUE) == TRUE) {
            ucnet <- net - rcp
        }
        else {
            ucnet <- net
        }
    }
    if (isTRUE(directed == TRUE) == FALSE) {
        bd <- multiplex::bundles(ucnet, loops = FALSE, lb2lb = FALSE, 
            collapse = FALSE)
    }
    else {
        bd <- multiplex::bundles(net, loops = FALSE, lb2lb = FALSE, 
            collapse = FALSE)
    }
    ifelse((is.na(dim(net)[3]) == TRUE | isTRUE(dim(net)[3] == 
        1L) == TRUE), r <- 1L, r <- length(bd[[1]]))
    bds <- multiplex::summaryBundles(bd, byties = TRUE)
    ifelse(isTRUE(length(bds) == 0) == TRUE, showAtts <- FALSE, 
        NA)
    ifelse(missing(ecol) == TRUE, ecol <- grDevices::gray.colors(r), 
        NA)
    ifelse(missing(lty) == TRUE, lty <- 1:r, NA)
    if ((is.na(dim(net)[3]) == TRUE | isTRUE(dim(net)[3] == 1) == 
        TRUE)) {
        Lt <- lty[1]
        vecol <- ecol[1]
    }
    else {
        ifelse(isTRUE(length(lty) == 1L) == TRUE, Lt <- 1:r, 
            Lt <- rep(lty, r)[1:r])
        ifelse(isTRUE(length(ecol) == 1L) == TRUE, vecol <- rep(ecol, 
            dim(net)[3]), vecol <- rep(ecol, dim(net)[3])[1:dim(net)[3]])
        if (isTRUE(length(lty) == length(Lt)) == FALSE) {
            Ltc <- seq_along(vecol)
        }
        else {
            ifelse(isTRUE(seq(lty) == lty) == TRUE, Ltc <- Lt, 
                Ltc <- 1:r)
        }
    }
    if (is.null(clu) == FALSE) {
        if (is.vector(clu) == FALSE) 
            stop("'clu' must be a vector")
        if (is.character(clu) == TRUE) {
            tmpclu <- clu
            for (i in 1:nlevels(factor(clu))) {
                clu[which(levels(factor(tmpclu))[i] == clu)] <- i
            }
            rm(i)
            clu <- methods::as(clu, "numeric")
            rm(tmpclu)
        }
        clu[which(is.na(clu))] <- 0
        nclu <- nlevels(factor(clu))
    }
    else {
        nclu <- 1L
    }
    flgcx <- FALSE
    if (is.null(cex) == TRUE) {
        if (isTRUE(length(bds) == 0) == TRUE) {
            cex <- 1L
        }
        else {
            cex <- length(bds[[1]])/2L
            if (isTRUE(length(bds) > 1L) == TRUE) {
                for (i in 2:length(bds)) ifelse(isTRUE(cex < 
                  (length(bds[[i]])/2L)) == TRUE, cex <- (length(bds[[i]])/2L), 
                  NA)
            }
            cex <- ceiling(cex)
        }
    }
    else {
        if (is.vector(cex) == FALSE) 
            stop("'cex' must be a vector")
        cex[which(is.na(cex))] <- 0
    }
    if (isTRUE(length(cex) == 1L) == TRUE) {
        cex <- rep(cex, n)
    }
    else if (isTRUE(length(cex) == nclu) == TRUE | isTRUE(length(cex) == 
        n) == TRUE) {
        flgcx <- TRUE
    }
    else {
        cex <- rep(cex[1], n)
    }
    if (isTRUE(flgcx == TRUE) == TRUE) {
        ifelse(isTRUE(max(cex) >= 10L) == TRUE, mxc <- 9L, mxc <- max(cex))
        cex <- round(((cex - min(cex))/(max(cex) - min(cex))) * 
            mxc) + 1L
        ifelse(isTRUE(min(cex) == 0) == TRUE, cex <- cex + 1L + 
            (2L/n), NA)
        rm(mxc)
    }
    else {
        ifelse(isTRUE(max(cex) >= 21L) == TRUE, cex <- 20L, NA)
    }
    if (is.null(pch) == TRUE) {
        pch <- rep(1, n)
    }
    if (isTRUE(length(pch) == 1L) == TRUE) {
        pch <- rep(pch, n)
    }
    else if (isTRUE(length(pch) == nclu) == TRUE) {
        if (identical(pch, clu) == FALSE) {
            tmppch <- rep(0, n)
            for (i in 1:nclu) {
                tmppch[which(clu == as.numeric(levels(factor(clu))[i]))] <- pch[i]
            }
            rm(i)
            pch <- tmppch
            rm(tmppch)
        }
    }
    else if (isTRUE(length(pch) != n) == TRUE) {
        pch <- rep(pch[1], n)
    }
    if (missing(vcol) == TRUE) {
        vcol <- grDevices::gray.colors(nclu)
        ifelse(missing(col) == TRUE, NA, vcol <- col)
    }
    else {
        if (isTRUE(length(vcol) == 1L) == TRUE) {
            vcol <- rep(vcol, n)
        }
        else if (isTRUE(length(vcol) == nclu) == TRUE) {
            if (identical(vcol, clu) == FALSE) {
                tmpvcol <- rep(0, n)
                for (i in 1:nclu) {
                  tmpvcol[which(clu == as.numeric(levels(factor(clu))[i]))] <- vcol[i]
                }
                rm(i)
                vcol <- tmpvcol
                rm(tmpvcol)
            }
        }
        else if (isTRUE(length(vcol) != n) == TRUE & isTRUE(nclu == 
            1) == TRUE) {
            vcol <- rep(vcol[1], n)
        }
        vcol[which(is.na(vcol))] <- graphics::par()$bg
        vcol[which(vcol == 0)] <- graphics::par()$bg
    }
    if (isTRUE(any(pch %in% 21:25)) == TRUE) {
        if (missing(vcol0) == TRUE) {
            vcol0 <- vcol
        }
        else {
            ifelse(missing(vcol0) == TRUE, NA, vcol0[which(is.na(vcol0))] <- 1)
        }
        if (isTRUE(length(vcol0) == 1L) == TRUE) {
            vcol0 <- rep(vcol0, n)
        }
        else if (isTRUE(length(vcol0) == nclu) == TRUE) {
            if (identical(vcol0, clu) == FALSE) {
                tmpvcol0 <- rep(0, n)
                for (i in 1:nclu) {
                  tmpvcol0[which(clu == as.numeric(levels(factor(clu))[i]))] <- vcol0[i]
                }
                rm(i)
                vcol0 <- tmpvcol0
                rm(tmpvcol0)
            }
        }
        else if (isTRUE(length(vcol0) != n) == TRUE | isTRUE(nclu == 
            1) == TRUE) {
            vcol0 <- rep(vcol0[1], n)
        }
    }
    else {
        vcol0 <- vcol
    }
    if (isTRUE(flgcx == TRUE) == FALSE) {
        ifelse(isTRUE(directed == TRUE) == TRUE, fds <- 130L, 
            fds <- 140L)
    }
    else if (isTRUE(flgcx == TRUE) == TRUE) {
        fds <- 85L + (n * 2L)
    }
    if (is.null(coord) == FALSE) {
        if (isTRUE(nrow(coord) == n) == FALSE) 
            stop("Length of 'coord' does not match network order.")
        flgcrd <- TRUE
    }
    else if (is.null(coord) == TRUE) {
        flgcrd <- FALSE
        switch(match.arg(layout), stress = {
            coord <- stss(net, seed = seed, maxiter = maxiter)
            ifelse(isTRUE(flgcx == TRUE) == TRUE, fds <- fds - 
                15L, NA)
        }, circ = {
            coord <- data.frame(X = sin(2L * pi * ((0:(n - 1L))/n)), 
                Y = cos(2L * pi * ((0:(n - 1L))/n)))
            ifelse(isTRUE(flgcx == TRUE) == TRUE, fds <- fds - 
                10L, NA)
        }, rand = {
            set.seed(seed)
            coord <- data.frame(X = round(stats::runif(n) * 1L, 
                5), Y = round(stats::runif(n) * 1L, 5))
        })
    }
    if (missing(rot) == FALSE) {
        coord[, 1:2] <- xyrt(coord[, 1:2], as.numeric(rot))
        coord[, 1:2] <- coord[, 1:2] - min(coord[, 1:2])
    }
    rat <- (max(coord[, 1]) - min(coord[, 1]))/(max(coord[, 2]) - 
        min(coord[, 2]))
    coord[, 1] <- (coord[, 1] - min(coord[, 1]))/(max(coord[, 
        1]) - min(coord[, 1]))
    ifelse(isTRUE(rat > 0) == TRUE, coord[, 2] <- ((coord[, 2] - 
        min(coord[, 2]))/(max(coord[, 2]) - min(coord[, 2]))) * 
        (1L/rat), coord[, 2] <- ((coord[, 2] - min(coord[, 2]))/(max(coord[, 
        2]) - min(coord[, 2]))) * (rat))
    if (isTRUE(flgcrd == TRUE) == TRUE) {
        if (isTRUE(ncol(coord) > 2) == TRUE) {
            lbgml <- tolower(as.vector(coord[, 3]))
            lbnet <- tolower(as.vector(dimnames(net)[[1]]))
            lbp <- vector()
            for (i in 1:nrow(coord)) {
                lbp <- append(lbp, which(lbnet[i] == lbgml))
            }
            rm(i)
            if (isTRUE(ncol(coord) > 3) == TRUE) {
                atgml <- as.vector(coord[, 4])
                atgml[which(is.na(atgml))] <- ""
                atts <- atgml[lbp]
            }
            nds <- data.frame(X = as.numeric(as.vector(coord[lbp, 
                1])), Y = as.numeric(as.vector(coord[lbp, 2])))
        }
        else {
            nds <- data.frame(X = as.numeric(as.vector(coord[, 
                1])), Y = as.numeric(as.vector(coord[, 2])))
        }
    }
    else {
        nds <- data.frame(X = as.numeric(as.vector(coord[, 1])), 
            Y = as.numeric(as.vector(coord[, 2])))
    }
    nds <- (2L/max(nds)) * nds
    if (isTRUE(flgcx == TRUE) == TRUE && isTRUE(sqrt(((max(nds[, 
        1]) - min(nds[, 1])) * (max(nds[, 2]) - min(nds[, 2])))/nrow(nds)) < 
        (1/3)) == TRUE) {
        nds <- nds * (2.223 - (4.45 * (sqrt(((max(nds[, 1]) - 
            min(nds[, 1])) * (max(nds[, 2]) - min(nds[, 2])))/n))))
    }
    else {
        nds <- nds * (0.5)
    }
    ifelse(isTRUE(max(cex) < 2) == TRUE, NA, bwd <- bwd * 0.75)
    ifelse(missing(mirrorX) == FALSE && isTRUE(mirrorX == TRUE) == 
        TRUE, nds[, 1] <- nds[, 1] * cos(pi) - nds[, 2] * sin(pi), 
        NA)
    ifelse(missing(mirrorY) == FALSE && isTRUE(mirrorY == TRUE) == 
        TRUE, nds[, 2] <- nds[, 2] * cos(pi) - nds[, 1] * sin(pi), 
        NA)
    opm <- graphics::par()$mar
    ifelse(all(mar == c(5.1, 4.1, 4.1, 2.1)) == TRUE, mar <- rep(0, 
        4), NA)
    ifelse(is.null(main) == TRUE, graphics::par(mar = mar), graphics::par(mar = mar + 
        c(0, 0, 2, 0)))
    if (match.arg(layout) == "circ") {
        ylim <- c(min(nds[, 2]) - (max(cex)/50L), max(nds[, 2]) + 
            (max(cex)/50L))
        xlim <- c(min(nds[, 1]) - (max(cex)/50L), max(nds[, 1]) + 
            (max(cex)/50L))
    }
    else {
        ylim <- c(min(nds[, 2]) - (max(cex)/100L), max(nds[, 
            2]) + (max(cex)/100L))
        xlim <- c(min(nds[, 1]) - (max(cex)/100L), max(nds[, 
            1]) + (max(cex)/100L))
    }
    obg <- graphics::par()$bg
    graphics::par(bg = grDevices::adjustcolor(bg, alpha = alpha[3]))
    graphics::plot(nds, type = "n", axes = FALSE, xlab = "", 
        ylab = "", ylim = ylim, xlim = xlim, asp = asp, main = main, 
        cex.main = cex.main)
    tlbs <- vector()
    for (i in 1:length(attr(bds, "names"))) {
        ifelse(isTRUE(length(multiplex::dhc(attr(bds, "names")[i], 
            prsep = "")) > 4L) == TRUE, tlbs <- append(tlbs, 
            tolower(paste(multiplex::dhc(attr(bds, "names")[i], 
                prsep = "")[1:4], collapse = ""))), tlbs <- append(tlbs, 
            tolower(attr(bds, "names"))[i]))
    }
    rm(i)
    if (isTRUE(length(bds) > 0) == TRUE) {
        for (k in 1:length(bds)) {
            prs <- as.numeric(multiplex::dhc(bds[[k]]))
            if (isTRUE(directed == TRUE) == TRUE | isTRUE(collRecip == 
                TRUE) == FALSE) {
                pars <- as.matrix(nds[as.numeric(levels(factor(multiplex::dhc(bds[[k]])))), 
                  ])
            }
            else {
                pars <- as.matrix(nds[prs, ])
                ifelse(isTRUE(tlbs[k] == "txch") == TRUE, pars <- pars[c(1:2, 
                  4:3), ], NA)
                if (isTRUE(tlbs[k] %in% c("mixd", "full")) == 
                  TRUE && isTRUE(collRecip == TRUE) == TRUE) {
                  ifelse(isTRUE(prs[1] == prs[length(prs)]) == 
                    TRUE, pars[c(1:2, 3:nrow(pars)), ] <- pars[c(2:1, 
                    3:nrow(pars)), ], NA)
                }
                else {
                  NA
                }
            }
            rr <- length(bds[[k]])
            if (isTRUE(rr > 0L) == TRUE) {
                q <- which(tlbs[k] == attr(bd, "names"))
                if ((is.na(dim(net)[3]) == TRUE | isTRUE(dim(net)[3] == 
                  1L) == TRUE)) {
                  vlt <- rep(Lt, rr)
                  vecol <- rep(ecol[1], rr)
                  tbnd <- as.vector(unlist(bd[q]))
                  if (isTRUE(length(tbnd) > 0L) == TRUE) {
                    ifelse(isTRUE(any(tbnd %in% bds[[k]])) == 
                      TRUE, vlt <- append(vlt, rep(Lt, q)), NA)
                  }
                }
                else {
                  vlt <- vector()
                  for (i in seq_along(Lt)) {
                    tbnd <- as.vector(unlist(bd[[q]][i]))
                    if (isTRUE(length(tbnd) > 0L) == TRUE) {
                      ifelse(isTRUE(any(tbnd %in% bds[[k]])) == 
                        TRUE, vlt <- append(vlt, rep(Lt[i], length(which(tbnd %in% 
                        bds[[k]])))), NA)
                    }
                  }
                  rm(i)
                  if (isTRUE(length(lty) == 1L) == TRUE) {
                    vlt1 <- rep(lty, length(vlt))
                    vltc <- vlt
                  }
                  else {
                    vltc <- vector()
                    if (isTRUE(Lt == Ltc) == FALSE) {
                      for (i in seq_along(Ltc)) {
                        tbnd <- as.vector(unlist(bd[[q]][i]))
                        if (isTRUE(length(tbnd) > 0L) == TRUE) {
                          ifelse(isTRUE(any(tbnd %in% bds[[k]])) == 
                            TRUE, vltc <- append(vltc, rep(Ltc[i], 
                            length(which(tbnd %in% bds[[k]])))), 
                            NA)
                        }
                      }
                      rm(i)
                    }
                    else {
                      if (isTRUE(seq(lty) == lty) == TRUE) {
                        vltc <- vlt
                      }
                      else {
                        for (i in seq_along(lty)) {
                          vltc <- append(vltc, replace(vlt[which(vlt == 
                            lty[i])], vlt[which(vlt == lty[i])] != 
                            i, i))
                        }
                        rm(i)
                      }
                    }
                  }
                }
                if (isTRUE(flgcx == TRUE) == TRUE) {
                  if (isTRUE(directed == TRUE) == TRUE && isTRUE(tlbs[k] %in% 
                    c("asym", "tent", "txch", "mixd")) == TRUE) {
                    if (isTRUE(multiplex::dhc(bds[[k]][1])[1] > 
                      multiplex::dhc(bds[[k]][1])[2]) == TRUE) {
                      Prs <- vector()
                      length(Prs) <- length(prs)
                      Prs[which(1:length(prs)%%2L == 1L)] <- prs[which(1:length(prs)%%2L == 
                        0L)]
                      Prs[which(1:length(prs)%%2L == 0L)] <- prs[which(1:length(prs)%%2L == 
                        1L)]
                      cx <- cex[Prs]
                    }
                    else {
                      cx <- cex[prs]
                    }
                  }
                  else {
                    cx <- cex[prs]
                  }
                  ifelse(isTRUE(match.arg(layout) == "rand") == 
                    TRUE, NA, cx[which(cx == max(cx))] <- max(cx) * 
                    0.85)
                }
                else {
                  cx <- rep(cex[1], 2)
                }
                fds <- fds + (vedist * -10)
                if ((isTRUE(is.na(dim(net)[3]) == TRUE | dim(net)[3] == 
                  1) == TRUE)) {
                  mbnd(pars, rr, bds[[k]], vlt, cx, lwd, vecol, 
                    directed, asp, bwd, alfa, fds, flgcx)
                }
                else {
                  ifelse(isTRUE(length(lty) == 1L) == TRUE, mbnd(pars, 
                    rr, bds[[k]], vlt1, cx, lwd, vecol[vltc], 
                    directed, asp, bwd, alfa, fds, flgcx), mbnd(pars, 
                    rr, bds[[k]], vlt, cx, lwd, vecol[vltc], 
                    directed, asp, bwd, alfa, fds, flgcx))
                }
            }
            else {
                NA
            }
        }
        rm(k)
    }
    if (all(pch %in% 21:25) == TRUE) {
        suppressWarnings(graphics::points(nds[, 1], nds[, 2], 
            pch = pch, cex = cex, col = grDevices::adjustcolor(vcol0, 
                alpha = alpha[1]), bg = grDevices::adjustcolor(vcol, 
                alpha = alpha[1])))
    }
    else {
        suppressWarnings(graphics::points(nds[, 1], nds[, 2], 
            pch = pch, cex = cex, col = grDevices::adjustcolor(vcol, 
                alpha = alpha[1]), bg = grDevices::adjustcolor(vcol, 
                alpha = alpha[1])))
    }
    if (isTRUE(showLbs == TRUE) == TRUE) {
        if (is.null(tcex) == TRUE) {
            ifelse(isTRUE(max(cex) < 2) == TRUE, tcex <- cex * 
                0.66, tcex <- cex * 0.33)
        }
        if (isTRUE(pos == 0) == TRUE) {
            graphics::text(nds, labels = lbs, cex = tcex, adj = 0.5, 
                col = tcol)
        }
        else {
            graphics::text(nds, lbs, cex = tcex, pos = pos, col = tcol, 
                offset = (cex/4L), adj = c(0.5, 1))
        }
    }
    if (isTRUE(showAtts == TRUE) == TRUE) {
        if (is.null(coord) == FALSE && isTRUE(ncol(coord) > 3L) == 
            TRUE) {
            NA
        }
        else if (is.null(att) == FALSE) {
            if (is.array(att) == TRUE) {
                atts <- rep("", length(lbs))
                if (is.na(dim(att)[3]) == TRUE | isTRUE(dim(att)[3] == 
                  1) == TRUE) {
                  atts[which(diag(att) != 0)] <- lbat
                }
                else {
                  neta <- multiplex::zbind(net, att)
                  clss <- multiplex::expos(multiplex::rel.sys(neta, 
                    att = (dim(net)[3] + 1L):dim(neta)[3]), classes = TRUE)$Classes
                  attr(clss, "names")[which(attr(clss, "names") == 
                    "ALL")] <- multiplex::jnt(dimnames(att)[[3]], 
                    prsep = "")
                  for (i in 2:length(clss)) {
                    atts[which(lbs %in% clss[[i]])] <- attr(clss, 
                      "names")[i]
                  }
                  rm(i)
                }
            }
            else if (is.vector(att) == TRUE) {
                ifelse(isTRUE(length(att) == n) == TRUE, atts <- att, 
                  atts <- rep("", length(lbs)))
            }
            else {
                atts <- rep("", length(lbs))
                atts[which((att) != 0)] <- lbat
            }
        }
        else {
            atts <- rep("", nrow(nds))
        }
        if (isTRUE(pos == 4) == TRUE) {
            graphics::text(nds, labels = atts, cex = tcex, pos = 3, 
                col = tcol, offset = (cex/4L), adj = c(1.5, 2))
        }
        else if (isTRUE(length(pos) == 1) == TRUE) {
            graphics::text(nds, labels = atts, cex = tcex, pos = pos + 
                1L, col = tcol, offset = (cex/4L), adj = c(1.5, 
                2))
        }
    }
    graphics::par(mar = opm)
    graphics::par(bg = obg)
}
