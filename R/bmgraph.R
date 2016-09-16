bmgraph <-
function (net, layout = c("bip", "bip3", "bip3e", "bip4", "stress", 
    "rand", "circ"), coord = NULL, tcex = NULL, alpha = c(1, 
    1, 1), showLbs = TRUE, showAtts = TRUE, att = NULL, lbat = "1", 
    main = NULL, cex.main, bg, mar, cex, pos, lwd, lty, ecol, 
    vcol, vcol0, asp, directed, collRecip, seed = NULL, maxiter, 
    bwd, clu, pch, tcol, hds, vedist, rot, mirrorX, mirrorY, 
    col, ...) 
{
    net <- multiplex::dichot(net, c = 1L)
    ifelse(missing(directed) == TRUE, directed <- FALSE, NA)
    ifelse(missing(collRecip) == TRUE, collRecip <- TRUE, NA)
    ifelse(missing(tcol) == TRUE, tcol <- c(1, 1), NA)
    ifelse(missing(pch) == TRUE, pch <- 1:0, NA)
    if (missing(clu)) {
        ifelse(isTRUE(length(pch) > 1L) == TRUE, pch <- pch[1:2], 
            pch <- rep(pch, 2))
    }
    ifelse(isTRUE(dim(net)[1] == 1L) == TRUE && ((match.arg(layout) == 
        "bip3") | (match.arg(layout) == "bip4")), layout <- "bip", 
        NA)
    ifelse(isTRUE(dim(net)[2] == 1L) == TRUE && ((match.arg(layout) == 
        "bip3e") | (match.arg(layout) == "bip4")), layout <- "bip", 
        NA)
    nn <- dim(net)[1]
    mm <- dim(net)[2]
    if (is.na(dim(net)[3]) == TRUE) {
        if (isTRUE(directed == TRUE) == FALSE) {
            bmnet <- (rbind(cbind(matrix(0, ncol = nn, nrow = nn, 
                dimnames = list(rownames(net), rownames(net))), 
                net), cbind(t(net), matrix(0, ncol = mm, nrow = mm, 
                dimnames = list(colnames(net), colnames(net))))))
        }
        else if (isTRUE(directed == TRUE) == TRUE) {
            bmnet <- (rbind(cbind(matrix(0, ncol = nn, nrow = nn, 
                dimnames = list(rownames(net), rownames(net))), 
                net), cbind(matrix(0, ncol = nn, nrow = mm, dimnames = list(colnames(net), 
                rownames(net))), matrix(0, ncol = mm, nrow = mm, 
                dimnames = list(colnames(net), colnames(net))))))
        }
    }
    else {
        bmnetdf <- data.frame(matrix(ncol = (nn + mm)^2, nrow = 0))
        for (k in 1:dim(net)[3]) {
            if (isTRUE(directed == TRUE) == FALSE) {
                temp <- (rbind(cbind(matrix(0, ncol = nn, nrow = nn, 
                  dimnames = list(rownames(net[, , k]), rownames(net[, 
                    , k]))), net[, , k]), cbind(t(net[, , k]), 
                  matrix(0, ncol = mm, nrow = mm, dimnames = list(colnames(net[, 
                    , k]), colnames(net[, , k]))))))
            }
            else if (isTRUE(directed == TRUE) == TRUE) {
                temp <- (rbind(cbind(matrix(0, ncol = nn, nrow = nn, 
                  dimnames = list(rownames(net[, , k]), rownames(net[, 
                    , k]))), net[, , k]), cbind(matrix(0, ncol = nn, 
                  nrow = mm, dimnames = list(colnames(net[, , 
                    k]), rownames(net[, , k]))), matrix(0, ncol = mm, 
                  nrow = mm, dimnames = list(colnames(net[, , 
                    k]), colnames(net[, , k]))))))
            }
            bmnetdf[(nrow(bmnetdf) + 1), ] <- as.vector(temp)
        }
        rm(k)
        rm(temp)
        bmnet <- array(dim = c((nn + mm), (nn + mm), nrow(bmnetdf)))
        for (i in 1:nrow(bmnetdf)) {
            bmnet[, , i][1:(nn + mm)^2] <- as.numeric(bmnetdf[i, 
                ])
        }
        rm(i)
    }
    n <- dim(bmnet)[1]
    ifelse(is.null(dimnames(net)[[1]]) == TRUE, nlbs <- as.character(1:nn), 
        nlbs <- dimnames(net)[[1]])
    ifelse(is.null(dimnames(net)[[2]]) == TRUE, mlbs <- as.character(nn + 
        1:mm), mlbs <- dimnames(net)[[2]])
    lbs <- dimnames(bmnet)[[1]] <- dimnames(bmnet)[[2]] <- c(nlbs, 
        mlbs)
    if (is.na(dim(net)[3]) == FALSE && is.null(dimnames(net)[[3]]) == 
        FALSE) {
        ifelse((is.na(dim(net)[3]) == TRUE | isTRUE(dim(net)[3] == 
            1) == TRUE), NA, dimnames(bmnet)[[3]] <- dimnames(net)[[3]])
    }
    if (isTRUE(length(alpha) < 2) == TRUE) {
        alfa <- 1
        alpha <- rep(alpha, 3)
    }
    else {
        alfa <- alpha[2]
    }
    if (isTRUE(length(alpha) < 3) == TRUE) 
        alpha <- append(alpha, 0.1)
    ifelse(missing(lwd) == TRUE, lwd <- 1, NA)
    ifelse(missing(bg) == TRUE, bg <- graphics::par()$bg, NA)
    ifelse(missing(mar) == TRUE, mar <- graphics::par()$mar, 
        NA)
    ifelse(missing(cex.main) == TRUE, cex.main <- graphics::par()$cex.main, 
        NA)
    if (missing(bwd) == TRUE) {
        bwd <- 1L
    }
    else {
        ifelse(isTRUE(bwd > 1L) == TRUE, bwd <- 1L, NA)
        ifelse(isTRUE(bwd <= 0L) == TRUE, bwd <- 0L, NA)
    }
    if (!(missing(hds))) {
        ifelse(isTRUE(hds == 0L) == TRUE, hds <- 1L, NA)
    }
    else {
        hds <- 1L
    }
    ifelse(missing(rot) == TRUE, NA, rot <- rot * -1)
    if (isTRUE(directed == FALSE) == TRUE) {
        if ((is.na(dim(bmnet)[3]) == TRUE | isTRUE(dim(bmnet)[3] == 
            1) == TRUE)) {
            bmnet <- multiplex::dichot(bmnet, c = 1L)
            nt <- bmnet + t(bmnet)
            rcp <- multiplex::dichot(nt, c = 2L)
            if (isTRUE(collRecip == TRUE) == TRUE) {
                rcp[lower.tri(rcp)] <- 0L
                ucbmnet <- bmnet - rcp
            }
            else {
                ucbmnet <- bmnet
            }
        }
        else {
            ucbmnet <- multiplex::dichot(bmnet, c = 1L)
            for (k in 1:dim(bmnet)[3]) {
                nt <- bmnet[, , k] + t(bmnet[, , k])
                rcp <- multiplex::dichot(nt, c = 2L)
                if (isTRUE(collRecip == TRUE) == TRUE) {
                  rcp[upper.tri(rcp)] <- 0L
                  ucbmnet[, , k] <- bmnet[, , k] - rcp
                }
                else {
                  ucbmnet[, , k] <- bmnet[, , k]
                }
            }
            rm(k)
        }
    }
    else {
        NA
    }
    if (isTRUE(directed == FALSE) == TRUE) {
        ifelse((is.na(dim(bmnet)[3]) == TRUE | isTRUE(dim(bmnet)[3] == 
            1) == TRUE), bd <- multiplex::bundles(as.array(as.matrix(ucbmnet)), 
            loops = FALSE, lb2lb = FALSE, collapse = FALSE), 
            bd <- multiplex::bundles((ucbmnet), loops = FALSE, 
                lb2lb = FALSE, collapse = FALSE))
    }
    else {
        ifelse((is.na(dim(bmnet)[3]) == TRUE | isTRUE(dim(bmnet)[3] == 
            1) == TRUE), bd <- multiplex::bundles(as.array(as.matrix(bmnet)), 
            loops = FALSE, lb2lb = FALSE, collapse = FALSE), 
            bd <- multiplex::bundles((bmnet), loops = FALSE, 
                lb2lb = FALSE, collapse = FALSE))
    }
    ifelse((is.na(dim(bmnet)[3]) == TRUE | isTRUE(dim(bmnet)[3] == 
        1L) == TRUE), r <- 1L, r <- length(bd[[1]]))
    bds <- multiplex::summaryBundles(bd, byties = TRUE)
    ifelse(isTRUE(length(bds) == 0) == TRUE, showAtts <- FALSE, 
        NA)
    ifelse(missing(ecol) == TRUE, ecol <- grDevices::gray.colors(r), 
        NA)
    ifelse(missing(lty) == TRUE, lty <- 1:r, NA)
    if ((is.na(dim(bmnet)[3]) == TRUE | isTRUE(dim(bmnet)[3] == 
        1) == TRUE)) {
        Lt <- lty[1]
        vecol <- ecol[1]
    }
    else {
        ifelse(isTRUE(length(lty) == 1L) == TRUE, Lt <- 1:r, 
            Lt <- rep(lty, r)[1:r])
        ifelse(isTRUE(length(ecol) == 1L) == TRUE, vecol <- rep(ecol, 
            dim(bmnet)[3]), vecol <- rep(ecol, dim(bmnet)[3])[1:dim(bmnet)[3]])
        if (isTRUE(length(lty) == length(Lt)) == FALSE) {
            Ltc <- seq_along(vecol)
        }
        else {
            ifelse(isTRUE(seq(lty) == lty) == TRUE, Ltc <- Lt, 
                Ltc <- 1:r)
        }
    }
    if (!(missing(clu))) {
        flgclu <- TRUE
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
        flgclu <- FALSE
        clu <- c(rep(1, nn), rep(2, mm))
        nclu <- 2L
    }
    flgcx <- FALSE
    if (missing(cex) == TRUE) {
        if (isTRUE(match.arg(layout) == "stress") == TRUE | isTRUE(match.arg(layout) == 
            "rand") == TRUE) {
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
            cex <- ceiling(24L/max(dim(net)))
        }
    }
    else if (missing(cex) == FALSE) {
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
        ifelse(isTRUE(max(cex) >= 10L) == TRUE, mxc <- 9L, mxc <- max(cex))
        cex <- round(((cex - min(cex))/(max(cex) - min(cex))) * 
            mxc) + 1L
        ifelse(isTRUE(min(cex) == 0) == TRUE, cex <- cex + 1L + 
            (2L/n), NA)
        rm(mxc)
    }
    else {
        cex <- rep(cex[1], n)
        ifelse(isTRUE(max(cex) >= 21L) == TRUE, cex <- 20L, NA)
    }
    if (isTRUE(length(pch) == nclu) == TRUE) {
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
    if (missing(vcol) == TRUE) {
        ifelse(isTRUE(flgclu == TRUE) == TRUE, vcol <- grDevices::gray.colors(nclu), 
            vcol <- rep(1, 2))
        ifelse(missing(col) == TRUE, NA, vcol <- col)
    }
    else if (isTRUE(length(vcol) == 1L) == TRUE) {
        ifelse(isTRUE(vcol == 0) == TRUE, vcol <- "transparent", 
            NA)
        vcol <- rep(vcol, n)
    }
    if (isTRUE(length(vcol) == nclu) == TRUE) {
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
    m <- n
    if (is.null(coord) == FALSE) {
        if (isTRUE(nrow(coord) == n) == FALSE) 
            stop("Length of 'coord' does not match network order.")
        ifelse(isTRUE(flgcx == TRUE) == TRUE, fds <- 100L + (m * 
            2L), fds <- 110L)
        if (missing(rot) == FALSE) {
            coord[, 1:2] <- xyrt(coord[, 1:2], as.numeric(rot))
            coord[, 1:2] <- coord[, 1:2] - min(coord[, 1:2])
        }
        rat <- (max(coord[, 1]) - min(coord[, 1]))/(max(coord[, 
            2]) - min(coord[, 2]))
        coord[, 1] <- (coord[, 1] - min(coord[, 1]))/(max(coord[, 
            1]) - min(coord[, 1]))
        ifelse(isTRUE(rat > 0) == TRUE, coord[, 2] <- ((coord[, 
            2] - min(coord[, 2]))/(max(coord[, 2]) - min(coord[, 
            2]))) * (1L/rat), coord[, 2] <- ((coord[, 2] - min(coord[, 
            2]))/(max(coord[, 2]) - min(coord[, 2]))) * (rat))
        if (isTRUE(ncol(coord) > 2) == TRUE) {
            lbgml <- tolower(as.vector(coord[, 3]))
            lbbmnet <- tolower(as.vector(dimnames(bmnet)[[1]]))
            lbp <- vector()
            for (i in 1:nrow(coord)) {
                lbp <- append(lbp, which(lbbmnet[i] == lbgml))
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
        nds <- (2L/max(nds)) * nds
        if (isTRUE(flgcx == TRUE) == TRUE && isTRUE(area(nds) < 
            (1/3)) == TRUE) {
            nds <- nds * (2.223 - (4.45 * (sqrt(((max(nds[, 1]) - 
                min(nds[, 1])) * (max(nds[, 2]) - min(nds[, 2])))/n))))
        }
        else {
            nds <- nds * (0.5)
        }
        are <- 50L + (1/area(nds))
        m <- n
        ifelse(isTRUE(max(cex) < 2) == TRUE, fds <- fds + (stats::median(cex) * 
            2), NA)
    }
    else if (is.null(coord) == TRUE) {
        switch(match.arg(layout), stress = {
            ifelse(isTRUE(flgcx == TRUE) == TRUE, fds <- 110L - 
                (m * 2L), fds <- 90L)
            ifelse(missing(maxiter) == TRUE, maxiter <- 99L, 
                NA)
            if (missing(seed) == TRUE) {
                cds <- stss(as.matrix(multiplex::mnplx(bmnet)), 
                  seed = set.seed(NULL), maxiter = maxiter)
            } else {
                cds <- stss(as.matrix(multiplex::mnplx(bmnet)), 
                  seed = seed, maxiter = maxiter)
            }
            if (missing(rot) == FALSE) {
                cds <- xyrt(cds, as.numeric(rot))
            }
            if (isTRUE(m > 3) == TRUE) {
                rat <- (max(cds[, 1]) - min(cds[, 1]))/(max(cds[, 
                  2]) - min(cds[, 2]))
                cds[, 1] <- (cds[, 1] - min(cds[, 1]))/(max(cds[, 
                  1]) - min(cds[, 1]))
                ifelse(isTRUE(rat > 0) == TRUE, cds[, 2] <- ((cds[, 
                  2] - min(cds[, 2]))/(max(cds[, 2]) - min(cds[, 
                  2]))) * (1L/rat), cds[, 2] <- ((cds[, 2] - 
                  min(cds[, 2]))/(max(cds[, 2]) - min(cds[, 2]))) * 
                  (rat))
            }
            nds <- data.frame(X = as.numeric(as.vector(cds[, 
                1])), Y = as.numeric(as.vector(cds[, 2])))
            if (isTRUE(n > 6) == TRUE && isTRUE(area(nds) < (1/2.8)) == 
                TRUE) {
                nds <- nds * (2.223 - (4.45 * area(nds)))
            } else {
                nds <- nds * (1 - area(nds))
            }
            are <- 50L + (1/area(nds))
            ifelse(isTRUE(max(cex) < 2) == TRUE, NA, fds <- fds + 
                (mean(cex) * 3))
        }, bip = {
            act <- nrm(rng(nn))
            evt <- nrm(rng(mm))
            Act <- cbind(rep(0, nrow(net)), act)
            Evt <- cbind(rep(1, ncol(net)), evt)
            nds <- rbind(Act, Evt)
            nds[which(is.nan(nds))] <- 0.5
            nds[, 2] <- nds[, 2] * cos(pi) - nds[, 1] * sin(pi)
            rownames(nds) <- lbs
            if (missing(rot) == FALSE) {
                nds <- as.data.frame(xyrt(nds, as.numeric(rot)))
            }
            if (isTRUE(flgcx == TRUE) == TRUE) {
                nds <- nrm(nds)
                rat <- (max(nds[, 1]) - min(nds[, 1]))/(max(nds[, 
                  2]) - min(nds[, 2]))
                nds[, 1] <- (nds[, 1] - min(nds[, 1]))/(max(nds[, 
                  1]) - min(nds[, 1]))
                ifelse(isTRUE(rat > 0) == TRUE, nds[, 2] <- ((nds[, 
                  2] - min(nds[, 2]))/(max(nds[, 2]) - min(nds[, 
                  2]))) * (1L/rat), nds[, 2] <- ((nds[, 2] - 
                  min(nds[, 2]))/(max(nds[, 2]) - min(nds[, 2]))) * 
                  (rat))
                fds <- 120L
            } else {
                fds <- 60L
            }
            are <- 30L + (lwd * 10L)
        }, bip3 = {
            act1 <- nrm(rng(ceiling(nn/2)))
            act2 <- nrm(rng(floor(nn/2)))
            evt <- nrm(rng(mm))
            Act1 <- cbind(rep(0, ceiling(nrow(net)/2)), act1)
            Act2 <- cbind(rep(2, floor(nrow(net)/2)), act2)
            Evt <- cbind(rep(1, ncol(net)), evt)
            nds <- rbind(Act1, Act2, Evt)
            nds[which(is.nan(nds))] <- 0.5
            nds[, 2] <- nds[, 2] * cos(pi) - nds[, 1] * sin(pi)
            rownames(nds) <- lbs
            if (missing(rot) == FALSE) {
                nds <- as.data.frame(xyrt(nds, as.numeric(rot)))
            }
            if (isTRUE(flgcx == TRUE) == TRUE) {
                nds <- nrm(nds)
                rat <- (max(nds[, 1]) - min(nds[, 1]))/(max(nds[, 
                  2]) - min(nds[, 2]))
                nds[, 1] <- (nds[, 1] - min(nds[, 1]))/(max(nds[, 
                  1]) - min(nds[, 1]))
                ifelse(isTRUE(rat > 0) == TRUE, nds[, 2] <- ((nds[, 
                  2] - min(nds[, 2]))/(max(nds[, 2]) - min(nds[, 
                  2]))) * (1L/rat), nds[, 2] <- ((nds[, 2] - 
                  min(nds[, 2]))/(max(nds[, 2]) - min(nds[, 2]))) * 
                  (rat))
                fds <- 120L
            } else {
                fds <- 60L
            }
            are <- 30L + (lwd * 10L)
        }, bip3e = {
            act <- nrm(rng(nn))
            Act <- cbind(rep(1, nrow(net)), act)
            evt1 <- nrm(rng(ceiling(mm/2)))
            evt2 <- nrm(rng(floor(mm/2)))
            Evt1 <- cbind(rep(0, ceiling(ncol(net)/2)), evt1)
            Evt2 <- cbind(rep(2, floor(ncol(net)/2)), evt2)
            nds <- rbind(Act, Evt1, Evt2)
            nds[which(is.nan(nds))] <- 0.5
            nds[, 2] <- nds[, 2] * cos(pi) - nds[, 1] * sin(pi)
            rownames(nds) <- lbs
            if (missing(rot) == FALSE) {
                nds <- as.data.frame(xyrt(nds, as.numeric(rot)))
            }
            if (isTRUE(flgcx == TRUE) == TRUE) {
                nds <- nrm(nds)
                rat <- (max(nds[, 1]) - min(nds[, 1]))/(max(nds[, 
                  2]) - min(nds[, 2]))
                nds[, 1] <- (nds[, 1] - min(nds[, 1]))/(max(nds[, 
                  1]) - min(nds[, 1]))
                ifelse(isTRUE(rat > 0) == TRUE, nds[, 2] <- ((nds[, 
                  2] - min(nds[, 2]))/(max(nds[, 2]) - min(nds[, 
                  2]))) * (1L/rat), nds[, 2] <- ((nds[, 2] - 
                  min(nds[, 2]))/(max(nds[, 2]) - min(nds[, 2]))) * 
                  (rat))
                fds <- 120L
            } else {
                fds <- 60L
            }
            are <- 30L + (lwd * 10L)
        }, bip4 = {
            qd <- n/(n + 6)
            act1 <- nrm(rng(ceiling(nn/2)) + 1L)
            act2 <- nrm(rng(floor(nn/2)))
            evt1 <- nrm(rng(ceiling(mm/2)) - 1L)
            evt2 <- nrm(rng(floor(mm/2)))
            Act2 <- cbind((act2 * qd) + ((0.05/qd)), rep(0, floor(nn/2)))
            Evt2 <- cbind((evt2 * qd) + (0.95 - qd), rep(1, floor(mm/2)))
            Act1 <- cbind(rep(0, ceiling(nn/2)), (act1 * 0.85) + 
                0.075)
            Evt1 <- cbind(rep(1, ceiling(mm/2)), (evt1 * 0.85) + 
                0.075)
            nds <- as.data.frame(rbind(Act1, Act2, Evt1, Evt2))
            nds[which(is.na(nds))] <- 0.5
            nds[, 2] <- nds[, 2] * cos(pi) - nds[, 1] * sin(pi)
            if (missing(rot) == FALSE) {
                nds <- as.data.frame(xyrt(nds, as.numeric(rot)))
            }
            if (isTRUE(flgcx == TRUE) == TRUE) {
                nds <- nrm(nds)
                rat <- (max(nds[, 1]) - min(nds[, 1]))/(max(nds[, 
                  2]) - min(nds[, 2]))
                nds[, 1] <- (nds[, 1] - min(nds[, 1]))/(max(nds[, 
                  1]) - min(nds[, 1]))
                ifelse(isTRUE(rat > 0) == TRUE, nds[, 2] <- ((nds[, 
                  2] - min(nds[, 2]))/(max(nds[, 2]) - min(nds[, 
                  2]))) * (1L/rat), nds[, 2] <- ((nds[, 2] - 
                  min(nds[, 2]))/(max(nds[, 2]) - min(nds[, 2]))) * 
                  (rat))
                fds <- 120L
            } else {
                fds <- 120L
            }
            are <- 40L + (1/area(nds))
        }, rand = {
            ifelse(isTRUE(flgcx == TRUE) == TRUE, fds <- 100L + 
                (m * 2L), fds <- 100L + (n * 2L))
            if (missing(seed) == FALSE) {
                set.seed(seed)
            } else {
                NA
            }
            nds <- data.frame(X = round(stats::runif(n) * 1L, 
                5), Y = round(stats::runif(n) * 1L, 5))
            if (isTRUE(n == 3) == TRUE) {
                if (isTRUE(area(nds) < 1/6) == TRUE) {
                  nds <- nds * (n)
                } else if (isTRUE(area(nds) < 1/5) == TRUE) {
                  nds <- nds * (2)
                } else if (isTRUE(area(nds) < 1/4) == TRUE) {
                  nds <- nds * (1.5)
                } else {
                  NA
                }
            }
            are <- 50L + (1/area(nds))
            ifelse(isTRUE(max(cex) < 2) == TRUE, NA, fds <- fds + 
                (mean(cex) * 3))
            if (missing(rot) == FALSE) {
                nds <- as.data.frame(xyrt(nds, as.numeric(rot)))
            }
        }, circ = {
            nds <- data.frame(X = sin(2L * pi * ((0:(n - 1L))/n)), 
                Y = cos(2L * pi * ((0:(n - 1L))/n)))
            if (missing(rot) == FALSE) {
                nds <- as.data.frame(xyrt(nds, as.numeric(rot)))
            }
            if (isTRUE(flgcx == TRUE) == TRUE) {
                rat <- (max(nds[, 1]) - min(nds[, 1]))/(max(nds[, 
                  2]) - min(nds[, 2]))
                nds[, 1] <- (nds[, 1] - min(nds[, 1]))/(max(nds[, 
                  1]) - min(nds[, 1]))
                ifelse(isTRUE(rat > 0) == TRUE, nds[, 2] <- ((nds[, 
                  2] - min(nds[, 2]))/(max(nds[, 2]) - min(nds[, 
                  2]))) * (1L/rat), nds[, 2] <- ((nds[, 2] - 
                  min(nds[, 2]))/(max(nds[, 2]) - min(nds[, 2]))) * 
                  (rat))
                fds <- 120L
            } else {
                fds <- 70L
            }
        })
    }
    if (missing(pos) == TRUE) {
        flgpos <- TRUE
        ifelse(match.arg(layout) == "bip3" | match.arg(layout) == 
            "bip3e", pos <- c(2, 1, 4), pos <- c(2, 4))
        ifelse(match.arg(layout) == "bip4", pos <- c(2, 4, 3, 
            1), NA)
        ifelse(match.arg(layout) == "circ", pos <- c(4, 2), NA)
    }
    else {
        flgpos <- FALSE
        if (isTRUE(pos < 0L) == TRUE | isTRUE(pos > 4L) == TRUE) 
            stop("Invalid \"pos\" value.")
        if (isTRUE(length(pos) == 1) == TRUE) {
            ifelse(match.arg(layout) == "bip3" | match.arg(layout) == 
                "bip3e", pos <- rep(pos, 3), pos <- rep(pos, 
                2))
        }
        else if (isTRUE(length(pos) > 1) == TRUE) {
            ifelse(match.arg(layout) == "bip3" | match.arg(layout) == 
                "bip3e", pos <- pos[c(1, 2, 1)], pos <- pos[1:2])
        }
        else {
            NA
        }
    }
    if (missing(mirrorX) == FALSE && isTRUE(mirrorX == TRUE) == 
        TRUE) {
        nds[, 1] <- nds[, 1] * cos(pi) - nds[, 2] * sin(pi)
        ifelse(isTRUE(flgpos == TRUE) == TRUE && match.arg(layout) != 
            "bip4", pos <- pos[length(pos):1], NA)
        ifelse(isTRUE(flgpos == TRUE) == TRUE && match.arg(layout) == 
            "bip4", pos <- pos[c(2:1, 3:4)], NA)
    }
    else {
        NA
    }
    if (missing(mirrorY) == FALSE && isTRUE(mirrorY == TRUE) == 
        TRUE) {
        nds[, 2] <- nds[, 2] * cos(pi) - nds[, 1] * sin(pi)
        ifelse(isTRUE(flgpos == TRUE) == TRUE && match.arg(layout) == 
            "bip4", pos <- pos[c(1:2, 4:3)], NA)
    }
    else {
        NA
    }
    if (match.arg(layout) == "stress" | match.arg(layout) == 
        "rand" | match.arg(layout) == "circ") {
        xlim <- c(min(nds[, 1]) - (max(cex)/100L) - (0), max(nds[, 
            1]) + (max(cex)/100L) + (0))
        ylim <- c(min(nds[, 2]) - (max(cex)/100L), max(nds[, 
            2]) + (max(cex)/100L))
        ifelse(missing(asp) == TRUE, asp <- 1, NA)
    }
    else if (match.arg(layout) == "bip4") {
        xlim <- c(min(nds[, 1]) - (max(cex)/100L) - (0), max(nds[, 
            1]) + (max(cex)/100L) + (0))
        ylim <- c(min(nds[, 2]) - (max(cex)/100L), max(nds[, 
            2]) + (max(cex)/100L))
        ifelse(missing(asp) == TRUE, asp <- 1.2, NA)
    }
    else {
        mx <- ceiling(length(dhc(lbs[1:nn], ""))/nn) * 0.02
        xlim <- c(min(nds[, 1]) - (max(cex)/100L) - (mx), max(nds[, 
            1]) + (max(cex)/100L) + (mx))
        ylim <- c(min(nds[, 2]) - (max(cex)/100L), max(nds[, 
            2]) + (max(cex)/100L))
        if (match.arg(layout) == "bip") {
            ifelse(missing(asp) == TRUE, asp <- 2L, asp <- asp[1] * 
                2L)
        }
        else {
            ifelse(missing(asp) == TRUE, asp <- 2.5, asp <- asp[1] * 
                2.5)
        }
        if (missing(rot) == FALSE) {
            ifelse(isTRUE(abs(rot) >= 360L) == TRUE, agx <- abs(rot)%%360L, 
                agx <- abs(rot))
            if (isTRUE(agx >= 270L) == TRUE) {
                gx <- abs(90L + (agx%%270L))
            }
            else if (isTRUE(agx > 180L) == TRUE) {
                ifelse(isTRUE(agx > 180L) == TRUE, lx <- abs(180L - 
                  (agx%%180L)), lx <- agx%%180L)
                gx <- abs(90L - (lx%%90L))
            }
            else {
                ifelse(isTRUE(agx > 180L) == TRUE, gx <- abs(180L - 
                  (x%%180L)), gx <- agx%%180L)
            }
            ifelse(isTRUE(gx >= 90L) == TRUE, gxx <- abs(90L - 
                (gx%%90L)), gxx <- gx%%90L)
            xx <- seq(0L, 90L, 1L)
            yy <- c(seq(2L, 1L, length.out = 46), seq(1L, 0.5, 
                length.out = 46)[2:46])
            ifelse(isTRUE(gx > 0) == TRUE, asp <- stats::predict(stats::nls(yy ~ 
                (a * exp(-b * xx)), start = list(a = 2L, b = (2L * 
                log(2L)/2L))))[gxx], asp <- stats::predict(stats::nls(yy ~ 
                (a * exp(-b * xx)), start = list(a = 2L, b = (2L * 
                log(2L)/2L))))[1])
            if (isTRUE(gxx > 45) == TRUE && isTRUE(flgpos == 
                TRUE) == TRUE) {
                pos <- (pos + 1L)%%4L
                ifelse(missing(mirrorX) == FALSE && isTRUE(mirrorX == 
                  TRUE) == TRUE, pos <- (pos - 2L)%%4L, NA)
                ifelse(missing(mirrorY) == FALSE && isTRUE(mirrorY == 
                  TRUE) == TRUE, pos <- (pos - 2L)%%4L, NA)
            }
            else {
                NA
            }
        }
        else {
            NA
        }
        fds <- fds + (r * 2)
    }
    opm <- graphics::par()$mar
    ifelse(all(mar == c(5.1, 4.1, 4.1, 2.1)) == TRUE, mar <- rep(0, 
        4L), NA)
    ifelse(is.null(main) == TRUE, graphics::par(mar = mar), graphics::par(mar = mar + 
        c(0, 0, 2L, 0)))
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
                if ((is.na(dim(bmnet)[3]) == TRUE | isTRUE(dim(bmnet)[3] == 
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
                if ((match.arg(layout) == "bip" | match.arg(layout) == 
                  "bip3" | match.arg(layout) == "bip3e") && missing(vedist) == 
                  FALSE) {
                  ifelse(isTRUE(vedist > 1L) == TRUE, vedist <- 1L, 
                    NA)
                  ifelse(isTRUE(vedist < 0) == TRUE, fds <- fds + 
                    (vedist * -1), fds <- fds + (vedist * 1))
                }
                else if (match.arg(layout) == "bip4" && missing(vedist) == 
                  FALSE) {
                  ifelse(isTRUE(vedist > 1L) == TRUE, vedist <- 1L, 
                    NA)
                  ifelse(isTRUE(vedist < 0) == TRUE, fds <- fds - 
                    (vedist * -1), fds <- fds - (vedist * 1))
                }
                else {
                  NA
                }
                ifelse(isTRUE(hds != 0) == TRUE, are <- (1/hds) * 
                  50, NA)
                if ((isTRUE(is.na(dim(bmnet)[3]) == TRUE | dim(bmnet)[3] == 
                  1) == TRUE)) {
                  mbnd(pares = pars, r = rr, b = bds[[k]], vlt, 
                    cx, lwd, ecol = vecol, directed, asp, bwd, 
                    alfa, fds, flgcx)
                }
                else {
                  ifelse(isTRUE(length(lty) == 1L) == TRUE, mbnd(pares = pars, 
                    r = rr, b = bds[[k]], vlt1, cx, lwd, ecol = vecol[vltc], 
                    directed, asp, bwd, alfa, fds, flgcx), mbnd(pares = pars, 
                    r = rr, b = bds[[k]], vlt, cx, lwd, ecol = vecol[vltc], 
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
    if (isTRUE(length(tcol) > 2) == TRUE | isTRUE(length(tcol) == 
        1) == TRUE) {
        tcol <- rep(tcol[1], 2)
    }
    else {
        NA
    }
    if (isTRUE(showLbs == TRUE) == TRUE) {
        if (is.null(tcex) == TRUE) {
            ifelse(isTRUE(max(cex) < 2) == TRUE, tcex <- cex * 
                0.66, tcex <- cex * 0.33)
        }
        if (match.arg(layout) == "bip3") {
            if (isTRUE(pos[1] == 0) == TRUE) {
                graphics::text(nds[1:length(act1), ], labels = lbs[1:length(act1)], 
                  cex = tcex, adj = 0.5, col = tcol[1])
            }
            else {
                graphics::text(nds[1:length(act1), ], lbs[1:length(act1)], 
                  cex = tcex, pos = pos[1], col = tcol[1], offset = (cex/4L), 
                  adj = c(0.5, 1))
            }
            if (isTRUE(pos[2] == 0) == TRUE) {
                graphics::text(nds[(nn + 1L):dim(bmnet)[1], ], 
                  labels = lbs[(nn + 1L):dim(bmnet)[1]], cex = tcex, 
                  adj = 0.5, col = tcol[2])
            }
            else {
                graphics::text(nds[(nn + 1L):dim(bmnet)[1], ], 
                  lbs[(nn + 1L):dim(bmnet)[1]], cex = tcex, pos = pos[2], 
                  col = tcol[2], offset = (cex/4L), adj = c(0.5, 
                    1))
            }
            if (isTRUE(pos[3] == 0) == TRUE) {
                graphics::text(nds[(length(act1) + 1):nn, ], 
                  labels = lbs[(length(act1) + 1):nn], cex = tcex, 
                  adj = 0.5, col = tcol[1])
            }
            else {
                graphics::text(nds[(length(act1) + 1):nn, ], 
                  lbs[(length(act1) + 1):nn], cex = tcex, pos = pos[3], 
                  col = tcol[1], offset = (cex/4L), adj = c(0.5, 
                    1))
            }
        }
        else if (match.arg(layout) == "bip3e") {
            if (isTRUE(pos[2] == 0) == TRUE) {
                graphics::text(nds[1:length(act), ], labels = lbs[1:length(act)], 
                  cex = tcex, adj = 0.5, col = tcol[2])
            }
            else {
                graphics::text(nds[1:length(act), ], lbs[1:length(act)], 
                  cex = tcex, pos = pos[2], col = tcol[2], offset = (cex/4L), 
                  adj = c(0.5, 1))
            }
            if (isTRUE(pos[1] == 0) == TRUE) {
                graphics::text(nds[(nn + 1L):(nn + length(evt1)), 
                  ], labels = lbs[(nn + 1L):(nn + length(evt1))], 
                  cex = tcex, adj = 0.5, col = tcol[1])
            }
            else {
                graphics::text(nds[(nn + 1L):(nn + length(evt1)), 
                  ], lbs[(nn + 1L):(nn + length(evt1))], cex = tcex, 
                  pos = pos[1], col = tcol[1], offset = (cex/4L), 
                  adj = c(0.5, 1))
            }
            if (isTRUE(pos[3] == 0) == TRUE) {
                graphics::text(nds[(nn + 1 + length(evt1)):dim(bmnet)[1], 
                  ], labels = lbs[(nn + 1 + length(evt1)):dim(bmnet)[1]], 
                  cex = tcex, adj = 0.5, col = tcol[1])
            }
            else {
                graphics::text(nds[(nn + 1 + length(evt1)):dim(bmnet)[1], 
                  ], lbs[(nn + 1 + length(evt1)):dim(bmnet)[1]], 
                  cex = tcex, pos = pos[3], col = tcol[1], offset = (cex/4L), 
                  adj = c(0.5, 1))
            }
        }
        else if (match.arg(layout) == "bip4") {
            if (isTRUE(pos[1] == 0) == TRUE) {
                graphics::text(nds[1:length(act1), ], labels = lbs[1:length(act1)], 
                  cex = tcex, adj = 0.5, col = tcol[1])
                graphics::text(nds[(length(act1) + 1):nn, ], 
                  labels = lbs[(length(act1) + 1):nn], cex = tcex, 
                  adj = 0.5, col = tcol[1])
            }
            else {
                graphics::text(nds[1:length(act1), ], lbs[1:length(act1)], 
                  cex = tcex, pos = pos[1], col = tcol[1], offset = (cex/4L), 
                  adj = c(0.5, 1))
                graphics::text(nds[(length(act1) + 1):nn, ], 
                  labels = lbs[(length(act1) + 1):nn], cex = tcex, 
                  pos = pos[3], col = tcol[1], offset = (cex/4L), 
                  adj = c(0.5, 1))
            }
            if (isTRUE(pos[2] == 0) == TRUE) {
                graphics::text(nds[(nn + 1L):(nn + length(evt1)), 
                  ], labels = lbs[(nn + 1L):(nn + length(evt1))], 
                  cex = tcex, adj = 0.5, col = tcol[2])
                graphics::text(nds[(nn + 1 + length(evt1)):dim(bmnet)[1], 
                  ], labels = lbs[(nn + 1 + length(evt1)):dim(bmnet)[1]], 
                  cex = tcex, adj = 0.5, col = tcol[2])
            }
            else {
                graphics::text(nds[(nn + 1L):(nn + length(evt1)), 
                  ], labels = lbs[(nn + 1L):(nn + length(evt1))], 
                  cex = tcex, pos = pos[2], col = tcol[2], offset = (cex/4L), 
                  adj = c(0.5, 1))
                graphics::text(nds[(nn + 1 + length(evt1)):dim(bmnet)[1], 
                  ], lbs[(nn + 1 + length(evt1)):dim(bmnet)[1]], 
                  cex = tcex, pos = pos[4], col = tcol[2], offset = (cex/4L), 
                  adj = c(0.5, 1))
            }
        }
        else {
            if (isTRUE(length(pos) == 2) == TRUE) {
                if (isTRUE(pos[1] == 0) == TRUE) {
                  graphics::text(nds[1:nn, ], labels = lbs[1:nn], 
                    cex = tcex, adj = 0.5, col = tcol[1])
                }
                else {
                  graphics::text(nds[1:nn, ], lbs[1:nn], cex = tcex, 
                    pos = pos[1], col = tcol[1], offset = (cex/4L), 
                    adj = c(0.5, 1))
                }
                if (isTRUE(pos[2] == 0) == TRUE) {
                  graphics::text(nds[(nn + 1L):dim(bmnet)[1], 
                    ], labels = lbs[(nn + 1L):dim(bmnet)[1]], 
                    cex = tcex, adj = 0.5, col = tcol[2])
                }
                else {
                  graphics::text(nds[(nn + 1L):dim(bmnet)[1], 
                    ], lbs[(nn + 1L):dim(bmnet)[1]], cex = tcex, 
                    pos = pos[2], col = tcol[2], offset = (cex/4L), 
                    adj = c(0.5, 1))
                }
            }
            else if (isTRUE(length(pos) == 1) == TRUE) {
                if (isTRUE(pos == 0) == TRUE) {
                  graphics::text(nds, labels = lbs, cex = tcex, 
                    adj = 0.5, col = tcol[1])
                }
                else {
                  graphics::text(nds, lbs, cex = tcex, pos = pos, 
                    col = tcol[1], offset = (cex/4L), adj = c(0.5, 
                      1))
                }
            }
            else {
                NA
            }
        }
    }
    graphics::par(mar = opm)
    graphics::par(bg = obg)
    x <- NULL
    rm(x)
}
