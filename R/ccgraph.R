ccgraph <-
function (x, main = NULL, seed = 0, maxiter = 100, alpha = c(1, 
    1, 1), scope, loops, collRecip, undRecip, showLbs, cex.main, 
    conc, coord, clu, cex, lwd, pch, lty, bwd, bwd2, att, bg, 
    mar, pos, asp, ecol, vcol, vcol0, lbs, col, lbat, swp, swp2, 
    scl, mirrorX, mirrorY, mirrorD, mirrorL, mirrorV, mirrorH, 
    rot, hds, vedist, ffamily, fstyle, fsize, fcol, ...) 
{
    pclu <- NULL
    if (isTRUE("Semigroup" %in% attr(x, "class")) == TRUE) {
        ifelse(is.null(x$ord) == FALSE, n <- x$ord, n <- dim(x$S)[1])
        if (all(x$st %in% dimnames(x$S)[[1]]) == TRUE) {
            Lbs <- x$st
        }
        else {
            if (isTRUE(unique(unlist(x$S)) %in% dimnames(x$S)[[1]]) == 
                FALSE) 
                stop("Semigroup labels do not match table elements.")
            Lbs <- dimnames(x$S)[[1]]
        }
        if (any(is.na(x$gens)) == TRUE || is.null(x$gens) == 
            TRUE) {
            warning("Generators are not provided, and first element of 'x' is taken.")
            x$gens <- 1
        }
        else {
            NA
        }
        if (is.array(x$gens) == TRUE) {
            cgm <- array(0, dim = c(n, n, dim(x$gens)[3]), dimnames = c(list(Lbs, 
                Lbs), list(attr(x$gens, "dimnames")[[3]])))
        }
        else if (is.vector(x$gens) == TRUE) {
            if (identical(which(Lbs %in% x$gens), seq_len(length(x$gens))) == 
                FALSE) {
                pclu <- rep(1, n)
                pclu[which(Lbs %in% x$gens)] <- 0
                Sp <- as.data.frame(multiplex::perm(as.matrix(x$S), 
                  clu = pclu))
                cgm <- array(0, dim = c(n, n, length(x$gens)), 
                  dimnames = c(dimnames(Sp), list(x$gens)))
            }
            else {
                cgm <- array(0, dim = c(n, n, length(x$gens)), 
                  dimnames = c(list(Lbs, Lbs), list(x$gens)))
            }
        }
        else if (is.numeric(x$gens) == TRUE) {
            cgm <- array(0, dim = c(n, n, x$gens), dimnames = c(list(Lbs, 
                Lbs), list(seq_len(x$gens))))
        }
        if (is.null(pclu) == TRUE) {
            for (k in seq_len(dim(cgm)[3])) {
                for (i in seq_len(n)) {
                  cgm[i, which(x$S[i, k] == Lbs), k] <- 1
                }
                rm(i)
            }
            rm(k)
        }
        else if (is.null(pclu) == FALSE) {
            for (k in seq_len(dim(cgm)[3])) {
                for (i in seq_len(n)) {
                  cgm[i, which(Sp[i, k] == colnames(Sp)), k] <- 1
                }
                rm(i)
            }
            rm(k)
        }
    }
    else if (isTRUE(attr(x, "class") == "EdgeTable") == TRUE) {
        if (is.null(x$ET) == FALSE) {
            xet <- x$ET
        }
        else {
            ifelse(is.list(x) == TRUE, xet <- x[[1]], stop("'x' must be a data frame in a list."))
        }
        ifelse(is.null(x$ord) == FALSE, n <- x$ord, n <- nrow(xet))
        lb <- rownames(xet)
        ifelse(is.null(x$gens) == FALSE, gens <- x$gens, gens <- colnames(xet))
        ifelse(identical(as.numeric(gens), as.numeric(colnames(xet))) == 
            TRUE, cgm <- array(0, dim = c(n, n, length(gens)), 
            dimnames = list(lb, lb, gens)), cgm <- array(0, dim = c(n, 
            n, dim(gens)[3]), dimnames = list(lb, lb, attr(gens, 
            "dimnames")[[3]])))
        for (k in seq_len(dim(cgm)[3])) {
            for (i in seq_len(n)) {
                cgm[i, which(xet[i, k] == lb), k] <- 1
            }
            rm(i)
        }
        rm(k)
    }
    else {
        if (is.array(x) == FALSE) 
            stop("Data must be at least a stacked array of square matrices.")
        x <- multiplex::semigroup(x, type = "symbolic")
        cgm <- array(0, dim = c(x$ord, x$ord, dim(x$gens)[3]), 
            dimnames = c(list(x$st, x$st), list(attr(x$gens, 
                "dimnames")[[3]])))
        for (k in seq_len(dim(cgm)[3])) {
            for (i in seq_len(x$ord)) {
                cgm[i, which(x$S[i, k] == x$st), k] <- 1
            }
            rm(i)
        }
        rm(k)
    }
    net <- cgm
    ifelse(is.null(pclu) == TRUE, NA, net <- multiplex::perm(net, 
        clu = c(which(pclu == 0), which(pclu == 1))))
    ifelse(isTRUE(dim(net)[3] == 1) == TRUE, net <- net[, , 1], 
        NA)
    if (isTRUE(n == 1L) == TRUE) 
        stop("1-element semigroup detected, and not yet supported")
    ifelse(missing(undRecip) == FALSE && isTRUE(undRecip == TRUE) == 
        TRUE, undRecip <- TRUE, undRecip <- FALSE)
    ifelse(missing(loops) == FALSE && isTRUE(loops == FALSE) == 
        TRUE, loops <- FALSE, loops <- TRUE)
    ifelse(missing(collRecip) == FALSE && isTRUE(collRecip == 
        FALSE) == TRUE, collRecip <- FALSE, collRecip <- TRUE)
    ifelse(missing(conc) == FALSE && isTRUE(conc == TRUE) == 
        TRUE, conc <- TRUE, conc <- FALSE)
    ifelse(missing(mirrorH) == FALSE && isTRUE(mirrorH == TRUE) == 
        TRUE, mirrorY <- TRUE, NA)
    ifelse(missing(mirrorV) == FALSE && isTRUE(mirrorV == TRUE) == 
        TRUE, mirrorX <- TRUE, NA)
    if (missing(showLbs) == FALSE && isTRUE(showLbs == TRUE) == 
        TRUE) {
        showLbs <- TRUE
    }
    else if (missing(showLbs) == FALSE && isTRUE(showLbs == FALSE) == 
        TRUE) {
        showLbs <- FALSE
    }
    else {
        ifelse(is.null(dimnames(net)[[1]]) == FALSE, showLbs <- TRUE, 
            showLbs <- FALSE)
    }
    ifelse(missing(swp) == FALSE && isTRUE(swp == TRUE) == TRUE, 
        swp <- TRUE, swp <- FALSE)
    ifelse(missing(swp2) == FALSE && isTRUE(swp2 == TRUE) == 
        TRUE, swp2 <- TRUE, swp2 <- FALSE)
    if (missing(scope) == FALSE) {
        if (isTRUE(is.list(scope) == TRUE) == FALSE) 
            stop("\"scope\" should be a list or a vector of lists.")
        scope <- list(scope)
        ifelse(is.null(scope[[1]]) == TRUE, scope <- scope[2:length(scope)], 
            NA)
        if (isTRUE(length(scope) > 1L) == TRUE && isTRUE(names(scope[1]) == 
            "coord") == TRUE) {
            scope <- scope[rev(seq_len(length(scope)))]
            flgrev <- TRUE
        }
        else {
            flgrev <- FALSE
        }
        tmp <- scope[[1]]
        if (isTRUE(length(scope) > 1L) == TRUE && isTRUE(length(scope[[1]]) > 
            1L) == TRUE) {
            for (k in 2:length(scope)) {
                tmp[length(tmp) + 1L] <- as.list(scope[k])
                names(tmp)[length(tmp)] <- attr(scope[k], "names")
            }
            rm(k)
        }
        else if (isTRUE(length(scope) > 1L) == TRUE) {
            names(tmp) <- attr(scope[1], "names")
            for (k in 2:length(scope)) {
                if (is.list(scope[[k]]) == TRUE && is.data.frame(scope[[k]]) == 
                  FALSE) {
                  for (j in seq_len(length(scope[[k]]))) {
                    tmp[length(tmp) + 1L] <- as.list(scope[[k]][j])
                    names(tmp)[length(tmp)] <- attr(scope[[k]][j], 
                      "names")
                  }
                  rm(j)
                }
                else if (is.data.frame(scope[[k]]) == FALSE) {
                  tmp[length(tmp) + 1L] <- as.list(scope[k])
                  names(tmp)[length(tmp)] <- attr(scope[k], "names")
                }
                else if (is.data.frame(scope[[k]]) == TRUE) {
                  tmp[length(tmp) + 1L] <- as.vector(scope[k])
                  names(tmp)[length(tmp)] <- attr(scope[k], "names")
                }
                else {
                  NA
                }
            }
            rm(k)
        }
        else {
            tmp <- scope[[1]]
        }
        ifelse(isTRUE(flgrev == TRUE) == TRUE, scope <- tmp[rev(seq_len(length(tmp)))], 
            scope <- tmp)
        for (i in seq_len(length(scope))) {
            if (isTRUE(names(scope)[i] %in% c("seed", "main")) == 
                TRUE) {
                tmpi <- as.vector(scope[[i]])
                assign(names(scope)[i], get("tmpi"))
            }
            else {
                if (is.null((scope[[i]])) == FALSE) {
                  tmpi <- as.vector(scope[[i]])
                  ifelse(isTRUE(names(scope)[i] != "") == TRUE, 
                    assign(names(scope)[i], get("tmpi")), NA)
                }
            }
        }
        rm(i)
    }
    else {
        NA
    }
    ifelse(missing(asp) == TRUE, asp <- 1, NA)
    ifelse(missing(lwd) == TRUE, lwd <- 1, NA)
    ifelse(missing(pch) == TRUE, pch <- 21, NA)
    ifelse(missing(fcol) == TRUE, fcol <- 1, NA)
    ifelse(missing(bwd) == TRUE, bwd <- 1, NA)
    ifelse(isTRUE(bwd < 0L) == TRUE, bwd <- 0L, NA)
    ifelse(missing(bg) == TRUE, bg <- graphics::par()$bg, NA)
    ifelse(missing(cex.main) == TRUE, cex.main <- graphics::par()$cex.main, 
        NA)
    ifelse(missing(rot) == TRUE, NA, rot <- rot[1] * -1)
    if (isTRUE(length(alpha) < 2) == TRUE) {
        alfa <- 1
        alpha <- rep(alpha, 3)
    }
    else {
        alfa <- alpha[2]
    }
    if (isTRUE(length(alpha) < 3) == TRUE) 
        alpha <- append(alpha, 0.1)
    if (!(missing(hds)) && missing(scl) == TRUE) {
        if (isTRUE(hds > 1L) == TRUE) {
            hds <- (hds/1.5)
        }
        else if (isTRUE(hds < 1L) == TRUE) {
            hds <- (hds/(hds + 0.15))
        }
        else if (isTRUE(hds == 0L) == TRUE) {
            hds <- 0.01
        }
        else {
            NA
        }
    }
    else {
        ifelse(missing(scl) == TRUE, hds <- 1L, hds <- 1L * scl)
    }
    ifelse(isTRUE(dim(net)[1] > 8) == TRUE || isTRUE(lwd >= 3) == 
        TRUE, hds <- hds * 0.75, NA)
    ifelse(missing(scl) == TRUE, scl <- rep(1, 2), NA)
    ifelse(isTRUE(length(scl) == 1) == TRUE, scl <- rep(scl, 
        2), scl <- scl[1:2])
    ifelse(missing(vedist) == TRUE, vedist <- 0, NA)
    if (isTRUE(vedist > 10L) == TRUE) {
        vedist <- 10L
    }
    else if (isTRUE(vedist < (-10L)) == TRUE) {
        vedist <- -10L
    }
    if (missing(lbs) == TRUE) {
        ifelse(is.null(dimnames(net)[[1]]) == TRUE, lbs <- as.character(seq_len(dim(net)[1])), 
            lbs <- dimnames(net)[[1]])
    }
    else {
        NA
    }
    n <- dim(net)[1]
    ifelse(isTRUE(is.na(dim(net)[3]) == TRUE) == TRUE, z <- 1L, 
        z <- dim(net)[3])
    ifelse(isTRUE(swp == TRUE) == TRUE && isTRUE(z > 1L) == TRUE, 
        net <- net[, , rev(seq_len(z))], NA)
    netd <- multiplex::dichot(net, c = 1L)
    if (isTRUE(collRecip == TRUE) == TRUE) {
        if (isTRUE(z == 1L) == TRUE) {
            nt <- netd + t(netd)
            rcp <- multiplex::dichot(nt, c = 2L)
            rcp[lower.tri(rcp, diag = TRUE)] <- 0L
        }
        else {
            nt <- array(0L, dim = c(n, n, z))
            dimnames(nt)[[1]] <- dimnames(nt)[[2]] <- lbs
            dimnames(nt)[[3]] <- dimnames(net)[[3]]
            for (i in seq_len(z)) {
                nt[, , i] <- netd[, , i] + t(netd[, , i])
            }
            rm(i)
            rcp <- multiplex::dichot(nt, c = 2L)
            for (i in seq_len(z)) {
                rcp[, , i][lower.tri(rcp[, , i], diag = TRUE)] <- 0L
            }
            rm(i)
        }
        ucnet <- netd - rcp
    }
    else {
        ucnet <- netd
    }
    bd <- multiplex::bundles(netd, loops = loops, lb2lb = FALSE, 
        collapse = FALSE)
    ifelse(isTRUE(z == 1L) == TRUE, r <- 1L, r <- length(bd[[1]]))
    ifelse(isTRUE(sum(net) == 0) == TRUE && isTRUE(loops == TRUE) == 
        TRUE, bd$loop <- character(0), NA)
    bds <- multiplex::summaryBundles(bd, byties = TRUE)
    ifelse(missing(ecol) == TRUE, ecol <- grDevices::gray.colors(r, 
        start = 0.1, end = 0.5), NA)
    ifelse(missing(lty) == TRUE, lty <- seq_len(r), NA)
    if (isTRUE(z == 1L) == TRUE) {
        Lt <- lty[1]
        vecol <- ecol[1]
    }
    else {
        ifelse(isTRUE(length(ecol) == 1L) == TRUE, vecol <- rep(ecol, 
            z), vecol <- rep(ecol, z)[seq_len(z)])
        ifelse(isTRUE(length(lty) == 1L) == TRUE, Lt <- seq_len(r), 
            Lt <- rep(lty, r)[seq_len(r)])
        if (isTRUE(length(lty) == length(Lt)) == FALSE) {
            Ltc <- seq_along(vecol)
        }
        else {
            if (isTRUE(seq(lty) == lty) == TRUE) {
                Ltc <- Lt
            }
            else {
                ifelse(isTRUE(swp == TRUE) == TRUE, Ltc <- rev(seq_len(r)), 
                  Ltc <- seq_len(r))
            }
        }
    }
    vltz <- Lt
    if (missing(clu) == FALSE) {
        if (is.vector(as.vector(clu)) == FALSE) 
            stop("'clu' must be a vector")
        if (is.factor(clu) == TRUE) {
            tmpclu <- clu
            for (i in seq_len(nlevels(factor(clu)))) {
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
            for (i in seq_len(nlevels(factor(tmpclu)) - 1L)) {
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
    if (missing(cex) == TRUE && isTRUE(loops == FALSE) == TRUE) {
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
    else if (missing(cex) == TRUE) {
        cex <- 1L
    }
    if (isTRUE(length(cex) == 1L) == TRUE) {
        cex <- rep(cex, n)
    }
    else {
        if (is.vector(cex) == FALSE) 
            stop("'cex' must be a vector")
        cex[which(is.na(cex))] <- 0
        cex <- cex[seq_len(n)]
        flgcx <- TRUE
    }
    if (isTRUE(flgcx == TRUE) == TRUE && isTRUE(max(cex) > 10L) == 
        TRUE) {
        if (isTRUE(mean(cex) > 20L) == TRUE) {
            cex <- (((cex - min(cex))/(max(cex) - min(cex))) * 
                10L)
        }
        else {
            cex <- (cex/(norm(as.matrix(cex), type = "M"))) * 
                10L
        }
        ifelse(isTRUE(min(cex) == 0) == TRUE, cex <- cex + 1L + 
            (2L/n), NA)
    }
    else if (isTRUE(flgcx == FALSE) == TRUE) {
        ifelse(isTRUE(max(cex) >= 21L) == TRUE, cex <- 20L, NA)
    }
    else {
        NA
    }
    if (missing(fsize) == TRUE) {
        ifelse(isTRUE(max(cex) < 2) == TRUE, fsize <- cex * 0.66, 
            fsize <- cex * 0.33)
    }
    else {
        fsize <- fsize/10
    }
    ifelse(isTRUE(bwd > 1L) == TRUE, bwd <- 1L, NA)
    ifelse(isTRUE(max(cex) < 2) == TRUE, NA, bwd <- bwd * 0.75)
    if (isTRUE(length(pch) == 1L) == TRUE) {
        pch <- rep(pch, n)
    }
    else if (isTRUE(length(pch) == nclu) == TRUE) {
        if (identical(pch, clu) == FALSE) {
            tmppch <- rep(0, n)
            for (i in seq_len(nclu)) {
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
                for (i in seq_len(nclu)) {
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
                for (i in seq_len(nclu)) {
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
    ifelse(isTRUE(n > 20) == TRUE, ffds <- 0.2, ffds <- 0)
    fds <- 130L - (n * ffds)
    if (isTRUE(flgcx == TRUE) == TRUE) {
        fds <- fds - 10L
    }
    else if (isTRUE(flgcx == FALSE) == TRUE) {
        NA
    }
    if (isTRUE(max(scl) < 1) == TRUE) {
        fds <- fds - (1/(mean(scl)/30L))
    }
    else if (isTRUE(max(scl) > 1) == TRUE) {
        fds <- fds + (mean(scl) * 20L)
    }
    else {
        NA
    }
    if (missing(coord) == FALSE) {
        if (isTRUE(nrow(coord) == n) == FALSE) 
            stop("Length of 'coord' does not match network order.")
        flgcrd <- TRUE
        crd <- coord
    }
    else if (missing(coord) == TRUE) {
        flgcrd <- FALSE
        ifelse(isTRUE(conc == TRUE) == FALSE, crd <- frcd(netd, 
            seed = seed, maxiter = maxiter), crd <- conc(netd, 
            ...))
    }
    if (missing(rot) == FALSE) {
        crd[, 1:2] <- xyrt(crd[, 1:2], as.numeric(rot))
        crd[, 1:2] <- crd[, 1:2] - min(crd[, 1:2])
        cnt <- 1L
        ifelse(isTRUE(n == 2) == TRUE && isTRUE(rot == -90) == 
            TRUE, rot <- -89.9, NA)
    }
    else {
        cnt <- 0
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
    if (missing(mirrorL) == FALSE && isTRUE(mirrorL == TRUE) == 
        TRUE) {
        crd[, 1:2] <- xyrt(crd[, 1:2], as.numeric(45))
        crd[, 1:2] <- crd[, 1:2] - min(crd[, 1:2])
        crd[, 2] <- crd[, 2] * cos(pi) - crd[, 1] * sin(pi)
        crd[, 1:2] <- xyrt(crd[, 1:2], as.numeric(-45))
        crd[, 1:2] <- crd[, 1:2] - min(crd[, 1:2])
    }
    else if (missing(mirrorD) == FALSE && isTRUE(mirrorD == TRUE) == 
        TRUE) {
        crd[, 1:2] <- xyrt(crd[, 1:2], as.numeric(-45))
        crd[, 1:2] <- crd[, 1:2] - min(crd[, 1:2])
        crd[, 2] <- crd[, 2] * cos(pi) - crd[, 1] * sin(pi)
        crd[, 1:2] <- xyrt(crd[, 1:2], as.numeric(45))
        crd[, 1:2] <- crd[, 1:2] - min(crd[, 1:2])
    }
    else {
        NA
    }
    if (isTRUE(n > 1) == TRUE) {
        rat <- (max(crd[, 1]) - min(crd[, 1]))/(max(crd[, 2]) - 
            min(crd[, 2]))
        crd[, 1] <- (crd[, 1] - min(crd[, 1]))/(max(crd[, 1]) - 
            min(crd[, 1]))
        ifelse(isTRUE(rat > 0) == TRUE, crd[, 2] <- ((crd[, 2] - 
            min(crd[, 2]))/(max(crd[, 2]) - min(crd[, 2]))) * 
            (1L/rat), crd[, 2] <- ((crd[, 2] - min(crd[, 2]))/(max(crd[, 
            2]) - min(crd[, 2]))) * (rat))
    }
    else {
        NA
    }
    fds <- fds + (vedist * -10)
    if (isTRUE(flgcrd == TRUE) == TRUE && isTRUE(ncol(crd) > 
        2) == TRUE) {
        lbgml <- tolower(as.vector(crd[, 3]))
        lbnet <- tolower(as.vector(lbs))
        lbp <- vector()
        for (i in seq_len(nrow(crd))) {
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
    nds <- ((2L/max(nds * (0.75))) * (nds * 0.75)) * (0.5)
    mscl <- mean(scl)
    cex <- cex * mscl
    fsize <- fsize * mscl
    omr <- graphics::par()$mar
    omi <- graphics::par()$mai
    if (missing(mar) == TRUE) {
        mar <- c(0, 0, 0, 0)
    }
    else {
        mar <- omr
    }
    ifelse(is.null(main) == TRUE, graphics::par(mar = mar), graphics::par(mar = mar + 
        c(0, 0, cex.main, 0)))
    obg <- graphics::par()$bg
    graphics::par(bg = grDevices::adjustcolor(bg, alpha = alpha[3]))
    if (isTRUE(loops == TRUE) == TRUE) {
        ylim <- c(min(nds[, 2]) - (max(cex)/100L), max(nds[, 
            2]) + (max(cex)/70L))
        xlim <- c(min(nds[, 1]) - (max(cex)/100L), max(nds[, 
            1]) + (max(cex)/100L))
    }
    else if (isTRUE(flgcx == TRUE) == TRUE) {
        ylim <- c(min(nds[, 2]) - (max(cex)/500L), max(nds[, 
            2]) + (max(cex)/500L))
        xlim <- c(min(nds[, 1]) - (max(cex)/500L), max(nds[, 
            1]) + (max(cex)/500L))
    }
    else {
        ylim <- c(min(nds[, 2]) - ((cex[1])/200L), max(nds[, 
            2]) + ((cex[1])/200L))
        xlim <- c(min(nds[, 1]) - ((cex[1])/200L), max(nds[, 
            1]) + ((cex[1])/200L))
    }
    suppressWarnings(graphics::plot(nds, type = "n", axes = FALSE, 
        xlab = "", ylab = "", ylim = ylim, xlim = xlim, asp = asp, 
        main = main, cex.main = cex.main, ...))
    tlbs <- vector()
    if (isTRUE(length(bds) > 0) == TRUE) {
        for (i in seq_len(length(attr(bds, "names")))) {
            ifelse(isTRUE(length(multiplex::dhc(attr(bds, "names")[i], 
                sep = "")) > 4L) == TRUE, tlbs <- append(tlbs, 
                tolower(paste(multiplex::dhc(attr(bds, "names")[i], 
                  sep = "")[1:4], collapse = ""))), tlbs <- append(tlbs, 
                tolower(attr(bds, "names"))[i]))
        }
        rm(i)
    }
    if (isTRUE(loops == TRUE) == TRUE) {
        tlbslp <- tlbs
        tlbs <- tlbs[which(tlbs != "loop")]
    }
    else {
        NA
    }
    if (isTRUE(swp == TRUE) == TRUE) {
        Lt <- Lt[rev(seq_len(length(Lt)))]
        lwd <- lwd[length(lwd):1]
        alfa <- alfa[rev(seq_len(length(alfa)))]
    }
    if (isTRUE(collRecip == TRUE) == TRUE) {
        trcp <- multiplex::transf(rcp, type = "tolist")
    }
    else {
        NA
    }
    if (isTRUE(length(tlbs) > 0) == TRUE) {
        for (k in seq_len(length(tlbs))) {
            prs <- as.numeric(multiplex::dhc(bds[[k]]))
            pars <- as.matrix(nds[as.numeric(levels(factor(multiplex::dhc(bds[[k]])))), 
                ])
            rbds <- length(bds[[k]])
            if (isTRUE(rbds > 0L) == TRUE) {
                q <- which(tlbs[k] == attr(bd, "names"))
                if (isTRUE(z == 1L) == TRUE) {
                  vlt <- rep(Lt, rbds)
                  vecol <- rep(ecol[1], rbds)
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
                cx <- cex
                lw <- rep(lwd[1], rbds)
                ifelse(isTRUE(swp2 == TRUE) == TRUE && isTRUE(tlbs[k] %in% 
                  c("recp")) == TRUE, bds[[k]] <- multiplex::swp(bds[[k]]), 
                  NA)
                if (isTRUE(collRecip == TRUE) == TRUE && isTRUE(tlbs[k] %in% 
                  c("recp")) == TRUE) {
                  bw <- 0L
                  ifelse(isTRUE(undRecip == TRUE) == TRUE, hd <- 0L, 
                    hd <- hds)
                  lw <- (lw + (3L/lw)) * mscl
                }
                else if (isTRUE(collRecip == TRUE) == FALSE && 
                  isTRUE(tlbs[k] %in% c("recp")) == TRUE) {
                  ifelse(isTRUE(undRecip == TRUE) == TRUE, hd <- 0L, 
                    hd <- hds)
                }
                else {
                  bw <- bwd
                  hd <- hds
                  lw <- lw * mscl
                }
                if (isTRUE(collRecip == TRUE) == TRUE && isTRUE(tlbs[k] %in% 
                  c("recp")) == FALSE) {
                  flgcr <- numeric()
                  sbds <- multiplex::swp(bds[[k]])
                  if (any(sbds %in% unlist(trcp)) == TRUE) {
                    for (i in seq_len(z)) {
                      ifelse(any(sbds %in% trcp[[i]]) == TRUE, 
                        flgcr <- append(flgcr, as.numeric(i)), 
                        NA)
                    }
                    rm(i)
                  }
                }
                else {
                  flgcr <- rep(0L, z)
                }
                pars[, 1] <- pars[, 1] * scl[1]
                pars[, 2] <- pars[, 2] * scl[2]
                if (isTRUE(z == 1L) == TRUE) {
                  ccbnd(pars, rbds, bds[[k]], vlt, cx * mscl, 
                    lw, vecol, bw, alfa, fds, flgcx, flgcr, hd, 
                    n)
                }
                else {
                  ifelse(isTRUE(length(lty) == 1L) == TRUE, ccbnd(pars, 
                    rbds, bds[[k]], vlt1, cx * mscl, lw, vecol[vltc], 
                    bw, alfa, fds, flgcx, flgcr, hd, n), ccbnd(pars, 
                    rbds, bds[[k]], vlt, cx * mscl, lw, vecol[vltc], 
                    bw, alfa, fds, flgcx, flgcr, hd, n))
                }
            }
            else {
                NA
            }
        }
        rm(k)
    }
    else {
        NA
    }
    if (isTRUE(loops == TRUE) == TRUE) {
        ifelse(isTRUE(swp == TRUE) == TRUE && isTRUE(z == 2L) == 
            TRUE, Lt <- rev(Lt), NA)
        bdlp <- bd$loop
        ndss <- nds
        ndss[, 1] <- ndss[, 1] * scl[1]
        ndss[, 2] <- ndss[, 2] * scl[2]
        dz <- (rng(z) + abs(min(rng(z))))/(10L)
        if (isTRUE(z == 1L) == TRUE) {
            lp <- as.vector(which(diag(net) > 0))
            if (isTRUE(length(lp) > 0) == TRUE) {
                for (i in seq_len(length(lp))) {
                  if (isTRUE(n < 3) == TRUE) {
                    dcx <- (cex[lp[i]] * 0.0075)
                    lpsz <- (cex[lp[i]] * 0.005) - (dz)
                  }
                  else {
                    dcx <- (cex[lp[i]] * 0.01)
                    lpsz <- (cex[lp[i]] * 0.0075) - (dz)
                  }
                  hc(ndss[lp[i], 1], ndss[lp[i], 2] + (dcx), 
                    lpsz, col = vecol, lty = Lt, lwd = lwd)
                }
                rm(i)
            }
            else {
                NA
            }
        }
        else if (isTRUE(z > 1) == TRUE) {
            if (missing(bwd2) == TRUE) {
                NA
            }
            else {
                ifelse(isTRUE(bwd2 < 1L) == TRUE && isTRUE(bwd2 == 
                  0) == FALSE, bwd2 <- 1L, NA)
                ifelse(isTRUE(bwd2 > 2L) == TRUE, bwd2 <- 2L, 
                  NA)
                if (isTRUE(bwd2 == 0) == TRUE || (any(duplicated(unlist(bdlp))) == 
                  FALSE && isTRUE(bwd2 == 1L) == TRUE)) {
                  dz <- rep(0, z)
                }
                else {
                  ifelse(missing(bwd2) == TRUE, dz <- (rng(z) + 
                    abs(min(rng(z))))/(10L), dz <- (bwd2) * (rng(z) + 
                    abs(min(rng(z))))/(1L))
                }
            }
            ifelse(isTRUE(cx > 3L) == TRUE, fcex <- 3L, fcex <- floor(cx))
            for (k in seq_len(length(bdlp))) {
                lp <- as.numeric(unique(multiplex::dhc(bdlp)[k][[1]]))
                if (isTRUE(length(lp) > 0) == TRUE) {
                  for (i in seq_len(length(lp))) {
                    ifelse(isTRUE(fcex[lp[i]] <= 3L) == TRUE | 
                      isTRUE(n < 3L) == TRUE, dz <- dz * 0.75, 
                      NA)
                    if (isTRUE(n < 3L) == TRUE) {
                      dcx <- fcex[lp[i]]/110L
                      lpsz <- abs((fcex[lp[i]] * 0.007) - dz[k])
                    }
                    else {
                      dcx <- fcex[lp[i]]/100L
                      lpsz <- abs((fcex[lp[i]] * 0.0075) - dz[k])
                    }
                    ifelse(isTRUE(length(lty) == 1) == TRUE, 
                      hc(ndss[lp[i], 1], ndss[lp[i], 2] + (dcx), 
                        lpsz, col = grDevices::adjustcolor(vecol[k], 
                          alpha = alfa), lty = lty, lwd = lwd[k]), 
                      hc(ndss[lp[i], 1], ndss[lp[i], 2] + (dcx), 
                        lpsz, col = grDevices::adjustcolor(vecol[k], 
                          alpha = alfa), lty = Lt[k], lwd = lwd[k]))
                  }
                  rm(i)
                }
                else {
                  dz <- append(0, dz)
                }
            }
            rm(k)
        }
    }
    else {
        NA
    }
    if (all(pch %in% 21:25) == TRUE) {
        graphics::points(nds[, 1] * scl[1], nds[, 2] * scl[2], 
            pch = pch, cex = cex, col = grDevices::adjustcolor(vcol0, 
                alpha = alpha[1]), bg = grDevices::adjustcolor(vcol, 
                alpha = alpha[1]))
    }
    else {
        graphics::points(nds[, 1] * scl[1], nds[, 2] * scl[2], 
            pch = pch, cex = cex, col = grDevices::adjustcolor(vcol, 
                alpha = alpha[1]), bg = grDevices::adjustcolor(vcol, 
                alpha = alpha[1]))
    }
    if (isTRUE(showLbs == TRUE) == TRUE) {
        ndss <- nds
        ndss[, 1] <- ndss[, 1] * scl[1]
        ndss[, 2] <- ndss[, 2] * scl[2]
        ifelse(missing(ffamily) == FALSE && isTRUE(ffamily %in% 
            names(grDevices::postscriptFonts())) == TRUE, graphics::par(family = ffamily), 
            NA)
        if (isTRUE(length(pos) == 1) == TRUE) {
            if (isTRUE(pos == 0) == TRUE) {
                if (missing(fstyle) == TRUE || (missing(fstyle) == 
                  FALSE && isTRUE(fstyle %in% c("italic", "bold", 
                  "bolditalic") == FALSE))) {
                  graphics::text(ndss, labels = lbs, cex = fsize, 
                    adj = 0.5, col = fcol)
                }
                else if (missing(fstyle) == FALSE) {
                  if (isTRUE(fstyle == "italic") == TRUE) {
                    graphics::text(ndss, labels = as.expression(lapply(lbs, 
                      function(x) bquote(italic(.(x))))), cex = fsize, 
                      adj = 0.5, col = fcol)
                  }
                  else if (isTRUE(fstyle == "bold") == TRUE) {
                    graphics::text(ndss, labels = as.expression(lapply(lbs, 
                      function(x) bquote(bold(.(x))))), cex = fsize, 
                      adj = 0.5, col = fcol)
                  }
                  else if (isTRUE(fstyle == "bolditalic") == 
                    TRUE) {
                    graphics::text(ndss, labels = as.expression(lapply(lbs, 
                      function(x) bquote(bolditalic(.(x))))), 
                      cex = fsize, adj = 0.5, col = fcol)
                  }
                }
            }
            else {
                if (missing(fstyle) == TRUE || (missing(fstyle) == 
                  FALSE && isTRUE(fstyle %in% c("italic", "bold", 
                  "bolditalic") == FALSE))) {
                  graphics::text(ndss, lbs, cex = fsize, pos = pos, 
                    col = fcol, offset = (cex/4L), adj = c(0.5, 
                      1))
                }
                else if (missing(fstyle) == FALSE) {
                  if (isTRUE(fstyle == "italic") == TRUE) {
                    graphics::text(ndss, as.expression(lapply(lbs, 
                      function(x) bquote(italic(.(x))))), cex = fsize, 
                      pos = pos, col = fcol, offset = (cex/4L), 
                      adj = c(0.5, 1))
                  }
                  else if (isTRUE(fstyle == "bold") == TRUE) {
                    graphics::text(ndss, as.expression(lapply(lbs, 
                      function(x) bquote(bold(.(x))))), cex = fsize, 
                      pos = pos, col = fcol, offset = (cex/4L), 
                      adj = c(0.5, 1))
                  }
                  else if (isTRUE(fstyle == "bolditalic") == 
                    TRUE) {
                    graphics::text(ndss, as.expression(lapply(lbs, 
                      function(x) bquote(bolditalic(.(x))))), 
                      cex = fsize, pos = pos, col = fcol, offset = (cex/4L), 
                      adj = c(0.5, 1))
                  }
                }
            }
        }
        else if (isTRUE(length(pos) == n) == TRUE) {
            if (missing(fstyle) == TRUE || (missing(fstyle) == 
                FALSE && isTRUE(fstyle %in% c("italic", "bold", 
                "bolditalic") == FALSE))) {
                graphics::text(ndss, lbs, cex = fsize, pos = pos, 
                  col = fcol[1], offset = (cex/4L), adj = c(0.5, 
                    1))
            }
            else if (missing(fstyle) == FALSE) {
                if (isTRUE(fstyle == "italic") == TRUE) {
                  graphics::text(ndss, as.expression(lapply(lbs, 
                    function(x) bquote(italic(.(x))))), cex = fsize, 
                    pos = pos, col = fcol[1], offset = (cex/4L), 
                    adj = c(0.5, 1))
                }
                else if (isTRUE(fstyle == "bold") == TRUE) {
                  graphics::text(ndss, as.expression(lapply(lbs, 
                    function(x) bquote(bold(.(x))))), cex = fsize, 
                    pos = pos, col = fcol[1], offset = (cex/4L), 
                    adj = c(0.5, 1))
                }
                else if (isTRUE(fstyle == "bolditalic") == TRUE) {
                  graphics::text(ndss, as.expression(lapply(lbs, 
                    function(x) bquote(bolditalic(.(x))))), cex = fsize, 
                    pos = pos, col = fcol[1], offset = (cex/4L), 
                    adj = c(0.5, 1))
                }
            }
        }
        else {
            if (isTRUE(pos[1] == 0) == TRUE) {
                if (missing(fstyle) == TRUE || (missing(fstyle) == 
                  FALSE && isTRUE(fstyle %in% c("italic", "bold", 
                  "bolditalic") == FALSE))) {
                  graphics::text(ndss, labels = lbs, cex = fsize, 
                    adj = 0.5, col = fcol)
                }
                else if (missing(fstyle) == FALSE) {
                  if (isTRUE(fstyle == "italic") == TRUE) {
                    graphics::text(ndss, labels = as.expression(lapply(lbs, 
                      function(x) bquote(italic(.(x))))), cex = fsize, 
                      adj = 0.5, col = fcol)
                  }
                  else if (isTRUE(fstyle == "bold") == TRUE) {
                    graphics::text(ndss, labels = as.expression(lapply(lbs, 
                      function(x) bquote(bold(.(x))))), cex = fsize, 
                      adj = 0.5, col = fcol)
                  }
                  else if (isTRUE(fstyle == "bolditalic") == 
                    TRUE) {
                    graphics::text(ndss, labels = as.expression(lapply(lbs, 
                      function(x) bquote(bolditalic(.(x))))), 
                      cex = fsize, adj = 0.5, col = fcol)
                  }
                }
            }
            else {
                if (missing(fstyle) == TRUE || (missing(fstyle) == 
                  FALSE && isTRUE(fstyle %in% c("italic", "bold", 
                  "bolditalic") == FALSE))) {
                  graphics::text(ndss, lbs, cex = fsize, pos = pos[1], 
                    col = fcol, offset = (cex/4L), adj = c(0.5, 
                      1))
                }
                else if (missing(fstyle) == FALSE) {
                  if (isTRUE(fstyle == "italic") == TRUE) {
                    graphics::text(ndss, as.expression(lapply(lbs, 
                      function(x) bquote(italic(.(x))))), cex = fsize, 
                      pos = pos[1], col = fcol, offset = (cex/4L), 
                      adj = c(0.5, 1))
                  }
                  else if (isTRUE(fstyle == "bold") == TRUE) {
                    graphics::text(ndss, as.expression(lapply(lbs, 
                      function(x) bquote(bold(.(x))))), cex = fsize, 
                      pos = pos[1], col = fcol, offset = (cex/4L), 
                      adj = c(0.5, 1))
                  }
                  else if (isTRUE(fstyle == "bolditalic") == 
                    TRUE) {
                    graphics::text(ndss, as.expression(lapply(lbs, 
                      function(x) bquote(bolditalic(.(x))))), 
                      cex = fsize, pos = pos[1], col = fcol, 
                      offset = (cex/4L), adj = c(0.5, 1))
                  }
                }
            }
        }
    }
    graphics::par(mar = omr)
    graphics::par(bg = obg)
    graphics::par(lend = 0)
    graphics::par(mai = omi)
}
