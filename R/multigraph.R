multigraph <-
function (net, layout = c("circ", "force", "stress", "conc", 
    "rand"), outline, directed = TRUE, main = NULL, seed = NULL, 
    maxiter = 100, alpha = c(1, 1, 1), collRecip, showLbs, showAtts, 
    weighted, cex.main, coord, clu, cex, lwd, pch, lty, bwd, 
    tcol, tcex, att, bg, mar, pos, asp, ecol, vcol, vcol0, hds, 
    vedist, rot, mirrorX, mirrorY, col, lbat, drp, ...) 
{
    if (isTRUE(is.array(net) == TRUE) == FALSE) 
        stop("\"net\" should be an array.")
    ifelse(isTRUE(dim(net)[3] == 1) == TRUE, net <- net[, , 1], 
        NA)
    ifelse(missing(weighted) == FALSE && isTRUE(weighted == TRUE) == 
        TRUE, weighted <- TRUE, weighted <- FALSE)
    ifelse(missing(collRecip) == FALSE && isTRUE(collRecip == 
        FALSE) == TRUE, collRecip <- FALSE, collRecip <- TRUE)
    ifelse(missing(showLbs) == FALSE && isTRUE(showLbs == FALSE) == 
        TRUE, showLbs <- FALSE, showLbs <- TRUE)
    ifelse(missing(showAtts) == FALSE && isTRUE(showAtts == FALSE) == 
        TRUE, showAtts <- FALSE, showAtts <- TRUE)
    ifelse(isTRUE(directed == FALSE) == TRUE, directed <- FALSE, 
        NA)
    if (missing(outline) == FALSE) {
        if (isTRUE(is.list(outline) == TRUE) == FALSE) 
            stop("\"outline\" should be a list or a vector of lists.")
        outline <- list(outline)
        ifelse(is.null(outline[[1]]) == TRUE, outline <- outline[2:length(outline)], 
            NA)
        if (isTRUE(length(outline) > 1) == TRUE && isTRUE(names(outline[1]) == 
            "coord") == TRUE) {
            outline <- outline[length(outline):1]
            flgrev <- TRUE
        }
        else {
            flgrev <- FALSE
        }
        tmp <- outline[[1]]
        if (isTRUE(length(outline) > 1) == TRUE && isTRUE(length(outline[[1]]) > 
            1) == TRUE) {
            for (k in 2:length(outline)) {
                tmp[length(tmp) + 1L] <- as.list(outline[k])
                names(tmp)[length(tmp)] <- attr(outline[k], "names")
            }
            rm(k)
        }
        else if (isTRUE(length(outline) > 1) == TRUE) {
            names(tmp) <- attr(outline[1], "names")
            for (k in 2:length(outline)) {
                if (is.list(outline[[k]]) == TRUE && is.data.frame(outline[[k]]) == 
                  FALSE) {
                  for (j in 1:length(outline[[k]])) {
                    tmp[length(tmp) + 1L] <- as.list(outline[[k]][j])
                    names(tmp)[length(tmp)] <- attr(outline[[k]][j], 
                      "names")
                  }
                  rm(j)
                }
                else if (is.data.frame(outline[[k]]) == FALSE) {
                  tmp[length(tmp) + 1L] <- as.list(outline[k])
                  names(tmp)[length(tmp)] <- attr(outline[k], 
                    "names")
                }
                else if (is.data.frame(outline[[k]]) == TRUE) {
                  tmp[length(tmp) + 1L] <- as.vector(outline[k])
                  names(tmp)[length(tmp)] <- attr(outline[k], 
                    "names")
                }
                else {
                  NA
                }
            }
            rm(k)
        }
        else {
            tmp <- outline[[1]]
        }
        ifelse(isTRUE(flgrev == TRUE) == TRUE, outline <- tmp[length(tmp):1], 
            outline <- tmp)
        for (i in 1:length(outline)) {
            if (isTRUE(names(outline)[i] %in% c("seed", "main")) == 
                TRUE) {
                tmpi <- as.vector(outline[[i]])
                assign(names(outline)[i], get("tmpi"))
            }
            else {
                if (is.null((outline[[i]])) == FALSE) {
                  tmpi <- as.vector(outline[[i]])
                  ifelse(isTRUE(names(outline)[i] != "") == TRUE, 
                    assign(names(outline)[i], get("tmpi")), NA)
                }
            }
        }
        rm(i)
    }
    else {
        NA
    }
    n <- dim(net)[1]
    ifelse(isTRUE(is.na(dim(net)[3]) == TRUE) == TRUE, z <- 1, 
        z <- dim(net)[3])
    ifelse(is.null(dimnames(net)[[1]]) == TRUE, lbs <- as.character(1:n), 
        lbs <- dimnames(net)[[1]])
    if (missing(att) == FALSE && is.array(att) == TRUE) {
        if (isTRUE(n != dim(att)[1]) == TRUE) {
            warning("Dimensions in \"net\" and \"att\" differ. No attributes are shown.")
            showAtts <- FALSE
        }
    }
    ifelse(missing(asp) == TRUE, asp <- 1, NA)
    ifelse(missing(lwd) == TRUE, lwd <- 1, NA)
    ifelse(missing(pch) == TRUE, pch <- 21, NA)
    ifelse(missing(tcol) == TRUE, tcol <- 1, NA)
    ifelse(missing(bwd) == TRUE, bwd <- 1, NA)
    ifelse(isTRUE(bwd < 0L) == TRUE, bwd <- 0L, NA)
    ifelse(missing(bg) == TRUE, bg <- graphics::par()$bg, NA)
    ifelse(missing(mar) == TRUE, mar <- graphics::par()$mar, 
        NA)
    ifelse(missing(cex.main) == TRUE, cex.main <- graphics::par()$cex.main, 
        NA)
    if (isTRUE(length(alpha) < 2) == TRUE) {
        alfa <- 1
        alpha <- rep(alpha, 3)
    }
    else {
        alfa <- alpha[2]
    }
    if (isTRUE(length(alpha) < 3) == TRUE) 
        alpha <- append(alpha, 0.1)
    if (!(missing(hds))) {
        ifelse(isTRUE(hds == 0L) == TRUE, hds <- 0.01, NA)
    }
    else {
        hds <- 1L
    }
    ifelse(missing(vedist) == TRUE, vedist <- 0, NA)
    ifelse(missing(rot) == TRUE, NA, rot <- rot[1] * -1)
    ifelse(isTRUE(vedist > 1L) == TRUE, vedist <- 1L, NA)
    if (missing(drp) == FALSE && is.numeric(drp) == TRUE) {
        netdrp <- replace(net, net <= drp, 0)
        netd <- multiplex::dichot(netdrp, c = 1L)
    }
    else {
        netd <- multiplex::dichot(net, c = 1L)
        netdrp <- net
    }
    if (isTRUE(directed == FALSE) == TRUE && isTRUE(collRecip == 
        TRUE) == TRUE && isTRUE(weighted == TRUE) == TRUE) {
        if (isTRUE(z == 1L) == TRUE) {
            netdrp <- netdrp + t(netdrp)
        }
        else {
            for (i in 1:z) {
                netdrp[, , i] <- netdrp[, , i] + t(netdrp[, , 
                  i])
            }
            rm(i)
        }
    }
    else {
        NA
    }
    if (isTRUE(directed == FALSE) == TRUE && isTRUE(collRecip == 
        TRUE) == TRUE) {
        if (isTRUE(z == 1L) == TRUE) {
            nt <- netd + t(netd)
            rcp <- multiplex::dichot(nt, c = 2L)
            rcp[lower.tri(rcp)] <- 0L
        }
        else {
            nt <- array(0L, dim = c(n, n, z))
            dimnames(nt)[[1]] <- dimnames(nt)[[2]] <- lbs
            dimnames(nt)[[3]] <- dimnames(net)[[3]]
            for (i in 1:z) {
                nt[, , i] <- netd[, , i] + t(netd[, , i])
            }
            rm(i)
            rcp <- multiplex::dichot(nt, c = 2L)
            for (i in 1:z) rcp[, , i][lower.tri(rcp[, , i])] <- 0L
        }
        ucnet <- netd - rcp
    }
    else {
        ucnet <- netd
    }
    if (isTRUE(directed == TRUE) == FALSE) {
        bd <- multiplex::bundles(ucnet, loops = FALSE, lb2lb = FALSE, 
            collapse = FALSE)
    }
    else {
        bd <- multiplex::bundles(netd, loops = FALSE, lb2lb = FALSE, 
            collapse = FALSE)
    }
    ifelse(isTRUE(z == 1L) == TRUE, r <- 1L, r <- length(bd[[1]]))
    bds <- multiplex::summaryBundles(bd, byties = TRUE)
    ifelse(missing(ecol) == TRUE, ecol <- grDevices::gray.colors(r), 
        NA)
    if (isTRUE(weighted == TRUE) == TRUE) {
        ifelse(missing(lty) == TRUE, lty <- rep(1, r), lty <- rep(lty, 
            r))
    }
    else {
        ifelse(missing(lty) == TRUE, lty <- 1:r, NA)
    }
    if (isTRUE(z == 1L) == TRUE) {
        Lt <- lty[1]
        vecol <- ecol[1]
    }
    else {
        ifelse(isTRUE(length(lty) == 1L) == TRUE, Lt <- 1:r, 
            Lt <- rep(lty, r)[1:r])
        ifelse(isTRUE(length(ecol) == 1L) == TRUE, vecol <- rep(ecol, 
            z), vecol <- rep(ecol, z)[1:z])
        if (isTRUE(length(lty) == length(Lt)) == FALSE) {
            Ltc <- seq_along(vecol)
        }
        else {
            ifelse(isTRUE(seq(lty) == lty) == TRUE, Ltc <- Lt, 
                Ltc <- 1:r)
        }
    }
    vltz <- Lt
    if (missing(clu) == FALSE) {
        if (is.vector(as.vector(clu)) == FALSE) 
            stop("'clu' must be a vector")
        if (is.factor(clu) == TRUE) {
            tmpclu <- clu
            for (i in 1:nlevels(factor(clu))) {
                levels(clu) <- c(levels(clu), i)
                clu[which(levels(factor(tmpclu))[i] == clu)] <- i
            }
            rm(i)
            clu <- methods::as(as.vector(clu), "numeric")
            rm(tmpclu)
        }
        else if (is.character(clu) == TRUE) {
            tmpclu <- clu
            clu[which(clu == clu[1])] <- 1
            for (i in 1:(nlevels(factor(tmpclu)) - 1L)) {
                clu[which(clu == clu[which((clu %in% tmpclu) == 
                  TRUE)[(i - 0)]])] <- (i + 1L)
            }
            rm(i)
            clu[which((clu %in% tmpclu) == TRUE)] <- nlevels(factor(tmpclu))
            clu <- methods::as(as.vector(clu), "numeric")
            rm(tmpclu)
        }
        else {
            NA
        }
        clu[which(is.na(clu))] <- 0
        nclu <- nlevels(factor(clu))
    }
    else {
        nclu <- 1L
    }
    flgcx <- FALSE
    if (missing(cex) == TRUE) {
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
        cex <- (((cex - min(cex))/(max(cex) - min(cex))) * mxc) + 
            1L
        ifelse(isTRUE(min(cex) == 0) == TRUE, cex <- cex + 1L + 
            (2L/n), NA)
        rm(mxc)
    }
    else {
        ifelse(isTRUE(max(cex) >= 21L) == TRUE, cex <- 20L, NA)
    }
    if (missing(tcex) == TRUE) {
        ifelse(isTRUE(max(cex) < 2) == TRUE, tcex <- cex * 0.66, 
            tcex <- cex * 0.33)
    }
    else {
        NA
    }
    ifelse(isTRUE(weighted == FALSE) == TRUE && isTRUE(bwd > 
        1L) == TRUE, bwd <- 1L, NA)
    ifelse(isTRUE(max(cex) < 2) == TRUE, NA, bwd <- bwd * 0.75)
    if (isTRUE(length(pch) == 1L) == TRUE) {
        pch <- rep(pch, n)
    }
    else if (isTRUE(length(pch) == nclu) == TRUE) {
        if (identical(pch, clu) == FALSE) {
            tmppch <- rep(0, n)
            for (i in 1:nclu) {
                tmppch[which(clu == (levels(factor(clu))[i]))] <- pch[i]
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
                  tmpvcol[which(clu == (levels(factor(clu))[i]))] <- vcol[i]
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
                  tmpvcol0[which(clu == (levels(factor(clu))[i]))] <- vcol0[i]
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
    if (isTRUE(flgcx == FALSE) == TRUE) {
        ifelse(isTRUE(directed == TRUE) == TRUE, fds <- 130L, 
            fds <- 140L)
    }
    else if (isTRUE(flgcx == TRUE) == TRUE) {
        fds <- 120L
    }
    if (missing(coord) == FALSE) {
        if (isTRUE(nrow(coord) == n) == FALSE) 
            stop("Length of 'coord' does not match network order.")
        flgcrd <- TRUE
        crd <- coord
    }
    else if (missing(coord) == TRUE) {
        flgcrd <- FALSE
        switch(match.arg(layout), force = {
            crd <- frcd(netd, seed = seed, maxiter = maxiter)
            ifelse(isTRUE(flgcx == TRUE) == TRUE, fds <- fds - 
                15L, NA)
        }, stress = {
            crd <- stsm(netd, seed = seed, maxiter = maxiter, 
                ...)
            ifelse(isTRUE(flgcx == TRUE) == TRUE, fds <- fds - 
                15L, NA)
        }, circ = {
            crd <- data.frame(X = sin(2L * pi * ((0:(n - 1L))/n)), 
                Y = cos(2L * pi * ((0:(n - 1L))/n)))
            ifelse(isTRUE(flgcx == TRUE) == TRUE, fds <- fds - 
                10L, NA)
        }, conc = {
            crd <- conc(netd, ...)
            ifelse(isTRUE(flgcx == TRUE) == TRUE, fds <- fds - 
                10L, NA)
        }, rand = {
            set.seed(seed)
            crd <- data.frame(X = round(stats::runif(n) * 1L, 
                5), Y = round(stats::runif(n) * 1L, 5))
        })
    }
    if (missing(rot) == FALSE) {
        crd[, 1:2] <- xyrt(crd[, 1:2], as.numeric(rot))
        crd[, 1:2] <- crd[, 1:2] - min(crd[, 1:2])
        cnt <- 1L
    }
    else {
        cnt <- 0
    }
    if (isTRUE(flgcrd == FALSE) == TRUE) {
        if (match.arg(layout) == "circ" && missing(pos) == TRUE) {
            angl <- vector()
            length(angl) <- n
            for (i in 1:n) {
                ifelse((atan2((crd[i, 2] - cnt), (crd[i, 1] - 
                  cnt)) * (180L/pi)) >= 0, angl[i] <- atan2((crd[i, 
                  2] - cnt), (crd[i, 1] - cnt)) * (180L/pi), 
                  angl[i] <- ((atan2((crd[i, 2] - cnt), (crd[i, 
                    1] - cnt)) * (180L/pi))%%180L) + 180L)
            }
            rm(i)
            pos <- vector()
            for (i in 1:length(angl)) {
                if (isTRUE(65 < angl[i]) == TRUE && isTRUE(115 > 
                  angl[i]) == TRUE) {
                  pos <- append(pos, 3)
                }
                else if (isTRUE(115 <= angl[i]) == TRUE && isTRUE(245 >= 
                  angl[i]) == TRUE) {
                  pos <- append(pos, 2)
                }
                else if (isTRUE(245 < angl[i]) == TRUE && isTRUE(295 > 
                  angl[i]) == TRUE) {
                  pos <- append(pos, 1)
                }
                else {
                  pos <- append(pos, 4)
                }
            }
            rm(i)
        }
    }
    if (missing(pos) == TRUE) {
        pos <- 4
    }
    else {
        if (isTRUE(pos < 0L) == TRUE | isTRUE(pos > 4L) == TRUE) 
            stop("Invalid \"pos\" value.")
    }
    ifelse(missing(mirrorX) == FALSE && isTRUE(mirrorX == TRUE) == 
        TRUE, crd[, 1] <- crd[, 1] * cos(pi) - crd[, 2] * sin(pi), 
        mirrorX <- FALSE)
    ifelse(missing(mirrorY) == FALSE && isTRUE(mirrorY == TRUE) == 
        TRUE, crd[, 2] <- crd[, 2] * cos(pi) - crd[, 1] * sin(pi), 
        mirrorY <- FALSE)
    if (match.arg(layout) == "circ" && isTRUE(flgcrd == FALSE) == 
        TRUE) {
        if (isTRUE(mirrorX == TRUE) == TRUE && isTRUE(length(pos) == 
            n) == TRUE) {
            pos[which(pos == 2)] <- 0
            pos[which(pos == 4)] <- 2
            pos[which(pos == 0)] <- 4
        }
        else if (isTRUE(mirrorY == TRUE) == TRUE && isTRUE(length(pos) == 
            n) == TRUE) {
            pos[which(pos == 1)] <- 0
            pos[which(pos == 3)] <- 1
            pos[which(pos == 0)] <- 3
        }
        else {
            NA
        }
    }
    rat <- (max(crd[, 1]) - min(crd[, 1]))/(max(crd[, 2]) - min(crd[, 
        2]))
    crd[, 1] <- (crd[, 1] - min(crd[, 1]))/(max(crd[, 1]) - min(crd[, 
        1]))
    ifelse(isTRUE(rat > 0) == TRUE, crd[, 2] <- ((crd[, 2] - 
        min(crd[, 2]))/(max(crd[, 2]) - min(crd[, 2]))) * (1L/rat), 
        crd[, 2] <- ((crd[, 2] - min(crd[, 2]))/(max(crd[, 2]) - 
            min(crd[, 2]))) * (rat))
    if (isTRUE(flgcrd == TRUE) == TRUE && isTRUE(ncol(crd) > 
        2) == TRUE) {
        lbgml <- tolower(as.vector(crd[, 3]))
        lbnet <- tolower(as.vector(lbs))
        lbp <- vector()
        for (i in 1:nrow(crd)) {
            lbp <- append(lbp, which(lbnet[i] == lbgml))
        }
        rm(i)
        if (isTRUE(ncol(crd) > 3) == TRUE) {
            atgml <- as.vector(crd[, 4])
            atgml[which(is.na(atgml))] <- ""
            atts <- atgml[lbp]
        }
        nds <- data.frame(X = as.numeric(as.vector(crd[lbp, 1])), 
            Y = as.numeric(as.vector(crd[lbp, 2])))
    }
    else {
        nds <- data.frame(X = as.numeric(as.vector(crd[, 1])), 
            Y = as.numeric(as.vector(crd[, 2])))
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
    opm <- graphics::par()$mar
    ifelse(all(mar == c(5.1, 4.1, 4.1, 2.1)) == TRUE, mar <- rep(0, 
        4), NA)
    ifelse(is.null(main) == TRUE, graphics::par(mar = mar), graphics::par(mar = mar + 
        c(0, 0, cex.main, 0)))
    if (match.arg(layout) == "circ" && isTRUE(flgcrd == FALSE) == 
        TRUE) {
        ylim <- c(min(nds[, 2]) - (max(cex)/50L), max(nds[, 2]) + 
            (max(cex)/50L))
        xlim <- c(min(nds[, 1]) - (max(cex)/50L), max(nds[, 1]) + 
            (max(cex)/50L))
    }
    else if (isTRUE(flgcx == TRUE) == TRUE) {
        ylim <- c(min(nds[, 2]) - (max(cex)/500L), max(nds[, 
            2]) + (max(cex)/500L))
        xlim <- c(min(nds[, 1]) - (max(cex)/500L), max(nds[, 
            1]) + (max(cex)/500L))
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
            sep = "")) > 4L) == TRUE, tlbs <- append(tlbs, tolower(paste(multiplex::dhc(attr(bds, 
            "names")[i], sep = "")[1:4], collapse = ""))), tlbs <- append(tlbs, 
            tolower(attr(bds, "names"))[i]))
    }
    rm(i)
    ifelse(isTRUE(weighted == TRUE) == TRUE && isTRUE(max(netdrp) > 
        10L) == TRUE, fnnetdrp <- (norm(as.matrix(netdrp), type = "F")), 
        NA)
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
                if (isTRUE(z == 1L) == TRUE) {
                  vlt <- rep(Lt, rr)
                  vecol <- rep(ecol[1], rr)
                  tbnd <- as.vector(unlist(bd[q]))
                  if (isTRUE(length(tbnd) > 0L) == TRUE) {
                    ifelse(isTRUE(any(tbnd %in% bds[[k]])) == 
                      TRUE, vlt <- append(vlt, rep(Lt, q)), NA)
                    ifelse(isTRUE(any(tbnd %in% bds[[k]])) == 
                      TRUE, vltz <- append(vltz, rep(Lt, q)), 
                      NA)
                  }
                  vltc <- vlt[1]
                }
                else {
                  vlt <- vector()
                  for (i in seq_along(Lt)) {
                    tbnd <- as.vector(unlist(bd[[q]][i]))
                    if (isTRUE(length(tbnd) > 0L) == TRUE) {
                      ifelse(isTRUE(any(tbnd %in% bds[[k]])) == 
                        TRUE, vlt <- append(vlt, rep(Lt[i], length(which(tbnd %in% 
                        bds[[k]])))), NA)
                      ifelse(isTRUE(any(tbnd %in% bds[[k]])) == 
                        TRUE, vltz <- append(vltz, rep(Lt[i], 
                        length(which(tbnd %in% bds[[k]])))), 
                        NA)
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
                  if ((isTRUE(directed == TRUE) == TRUE && isTRUE(tlbs[k] %in% 
                    c("asym", "tent", "txch", "mixd")) == TRUE) | 
                    (isTRUE(collRecip == FALSE) == TRUE && isTRUE(tlbs[k] %in% 
                      c("asym", "tent", "txch", "mixd")) == TRUE)) {
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
                else if (isTRUE(flgcx == FALSE) == TRUE) {
                  cx <- rep(cex[1], 2)
                }
                fds <- fds + (vedist * -10)
                if (isTRUE(weighted == TRUE) == TRUE) {
                  Lw <- vector()
                  i <- 1
                  for (j in 1:length(bds[[k]])) {
                    qn <- c(prs[i], prs[(i + 1)])
                    if (isTRUE(directed == FALSE) == TRUE && 
                      isTRUE(collRecip == TRUE) == TRUE) {
                      ifelse(isTRUE(z == 1L) == TRUE, Lw <- append(Lw, 
                        netdrp[qn[1], qn[2]]), Lw <- append(Lw, 
                        netdrp[qn[1], qn[2], vltc[j]] + t(netdrp[qn[1], 
                          qn[2], vltc[j]])))
                    }
                    else if (isTRUE(directed == TRUE) == TRUE | 
                      isTRUE(collRecip == FALSE) == TRUE) {
                      ifelse(isTRUE(z == 1L) == TRUE, Lw <- append(Lw, 
                        netdrp[qn[1], qn[2]]), Lw <- append(Lw, 
                        netdrp[qn[1], qn[2], vltc[j]]))
                    }
                    i <- i + 2L
                  }
                  rm(j)
                  rm(i)
                  if (isTRUE(max(netdrp) > 10L) == TRUE) {
                    lwd <- (Lw/fnnetdrp) * (10L * 5L)
                  }
                  else {
                    lwd <- Lw
                  }
                }
                else {
                  lwd <- rep(lwd[1], rr)
                }
                if (isTRUE(z == 1L) == TRUE) {
                  mbnd(pars, rr, bds[[k]], vlt, cx, lwd, vecol, 
                    directed, asp, bwd, alfa, fds, flgcx, weighted)
                }
                else {
                  ifelse(isTRUE(length(lty) == 1L) == TRUE, mbnd(pars, 
                    rr, bds[[k]], vlt1, cx, lwd, vecol[vltc], 
                    directed, asp, bwd, alfa, fds, flgcx, weighted), 
                    mbnd(pars, rr, bds[[k]], vlt, cx, lwd, vecol[vltc], 
                      directed, asp, bwd, alfa, fds, flgcx, weighted))
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
    if (isTRUE(showLbs == TRUE) == TRUE && isTRUE(is.null(dimnames(net)[[1]]) == 
        FALSE)) {
        if (isTRUE(length(pos) == 1) == TRUE) {
            if (isTRUE(pos == 0) == TRUE) {
                graphics::text(nds, labels = lbs, cex = tcex, 
                  adj = 0.5, col = tcol)
            }
            else {
                graphics::text(nds, lbs, cex = tcex, pos = pos, 
                  col = tcol, offset = (cex/4L), adj = c(0.5, 
                    1))
            }
        }
        else if (isTRUE(length(pos) == n) == TRUE) {
            graphics::text(nds, lbs, cex = tcex, pos = pos, col = tcol[1], 
                offset = (cex/4L), adj = c(0.5, 1))
        }
        else {
            if (isTRUE(pos[1] == 0) == TRUE) {
                graphics::text(nds, labels = lbs, cex = tcex, 
                  adj = 0.5, col = tcol)
            }
            else {
                graphics::text(nds, lbs, cex = tcex, pos = pos[1], 
                  col = tcol, offset = (cex/4L), adj = c(0.5, 
                    1))
            }
        }
    }
    if (isTRUE(showAtts == TRUE) == TRUE) {
        if (isTRUE(flgcrd == TRUE) == TRUE && isTRUE(ncol(coord) > 
            3L) == TRUE) {
            NA
        }
        else {
            atts <- rep("", nrow(nds))
            if (missing(att) == FALSE) {
                if (is.array(att) == TRUE) {
                  if (is.na(dim(att)[3]) == TRUE | isTRUE(dim(att)[3] == 
                    1) == TRUE) {
                    ifelse(missing(lbat) == FALSE, atts[which((att) != 
                      0)] <- lbat, atts[which((att) != 0)] <- "1")
                  }
                  else {
                    if (missing(lbat) == FALSE) {
                      atts[which(diag(multiplex::mnplx(netd, 
                        diag.incl = TRUE)) != 0)] <- lbat
                    }
                    else {
                      dimnames(netd)[[3]] <- NULL
                      neta <- multiplex::zbind(netd, att)
                      clss <- multiplex::expos(multiplex::rel.sys(neta, 
                        att = (z + 1L):dim(neta)[3]), classes = TRUE)$Classes
                      attr(clss, "names")[which(attr(clss, "names") == 
                        "ALL")] <- multiplex::jnt(dimnames(att)[[3]], 
                        sep = "")
                      for (i in 2:length(clss)) {
                        atts[which(lbs %in% clss[[i]])] <- attr(clss, 
                          "names")[i]
                      }
                      rm(i)
                    }
                  }
                }
                else if (is.vector(att) == TRUE | is.factor(att) == 
                  TRUE) {
                  ifelse(isTRUE(length(att) == n) == TRUE, atts <- as.vector(att), 
                    atts <- rep("", length(lbs)))
                }
                else {
                  atts <- rep("", length(lbs))
                }
            }
            else {
                NA
            }
        }
        if (isTRUE(flgcx == FALSE) == TRUE) {
            graphics::text(nds, labels = atts, cex = tcex, pos = pos%%4 + 
                1L, col = tcol, offset = (cex/4L), adj = c(0.5, 
                1))
        }
        else if (isTRUE(flgcx == TRUE) == TRUE) {
            graphics::text(nds, labels = atts, cex = tcex, pos = pos%%4 + 
                1L, col = tcol, offset = (min(cex)/4L), adj = c(0.5, 
                1))
        }
    }
    graphics::par(mar = opm)
    graphics::par(bg = obg)
    graphics::par(lend = 0)
}
