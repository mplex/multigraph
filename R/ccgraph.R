## (GPL-3)  Antonio Rivero Ostoic, multiplex@post.com
ccgraph <-
function (x, main = NULL, seed = 0, maxiter = 100, alpha = c(1, 
    1, 1), scope, loops, collRecip, undRecip, showLbs, cex.main, 
    col.main, font.main, conc, coord, clu, cex, lwd, pch, lty, 
    bwd, bwd2, att, bg, mar, pos, asp, ecol, vcol, vcol0, lbs, 
    col, lbat, swp, swp2, rot, vedist, vldist, mirrorX, mirrorY, 
    mirrorD, mirrorL, mirrorV, mirrorH, ffamily, fstyle, fsize, 
    fcol, scl, hds, nr, gens, flip, ...) 
{
    pclu <- NULL
    if (is.matrix(x) == TRUE && missing(gens) == TRUE) 
        stop(paste("Semigroup generators in matrix ", sQuote(deparse(substitute(x))), 
            " are missing.", sep = ""))
    ifelse(is.matrix(x) == TRUE && missing(gens) == FALSE, x <- multiplex::as.semigroup(x, 
        gens = gens), NA)
    if (isTRUE(tolower(class(x)[1]) == "semigroup") == TRUE) {
        ifelse(is.null(x$ord) == FALSE, n <- x$ord, n <- dim(x$S)[1])
        if (any(dimnames(x$S)[[1]] != dimnames(x$S)[[2]]) == 
            TRUE) 
            stop("Dimnames in semigroup multiplication table differ.")
        if (all(x$st %in% dimnames(x$S)[[1]]) == TRUE) {
            Lbs <- x$st
        }
        else {
            Lbs <- dimnames(x$S)[[1]]
        }
        if (missing(gens) == FALSE) {
            x$gens <- gens
        }
        else if (any(is.na(x$gens)) == TRUE || is.null(x$gens) == 
            TRUE) {
            warning(paste("Generators of semigroup ", sQuote(deparse(substitute(x))), 
                " are missing and assume first element.", sep = ""), 
                call. = FALSE)
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
                message("Apply automorphism on semigroup multiplication table.")
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
            ifelse(is.list(x) == TRUE, xet <- x[[1]], stop(paste("For ", 
                sQuote("EdgeTable"), " option ", sQuote(deparse(substitute(x))), 
                " should be a data frame in a list object.", 
                sep = "")))
        }
        ifelse(is.null(x$ord) == FALSE, n <- x$ord, n <- nrow(xet))
        lb <- rownames(xet)
        if (missing(gens) == TRUE) {
            ifelse(is.null(x$gens) == FALSE, gens <- x$gens, 
                gens <- colnames(xet))
        }
        else {
            invisible(NA)
        }
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
            stop(paste(sQuote(deparse(substitute(x))), " data should be at least an array of square matrices.", 
                sep = ""))
        x <- multiplex::semigroup(x, type = "symbolic")
        ifelse(is.null(x$ord) == FALSE, n <- x$ord, n <- dim(x$S)[1])
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
    if (isTRUE(n == 1L) == TRUE) {
        message("1-element semigroup.")
        return(x)
    }
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
            stop(paste("Values in ", sQuote("scope"), " should be in a list or a vector of lists.", 
                sep = ""))
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
    if (missing(bg) == TRUE) {
        bg <- graphics::par()$bg
    }
    else if (missing(bg) == FALSE && isTRUE(bg == 0) == TRUE) {
        bg <- "#FFFFFF"
    }
    else {
        NA
    }
    ifelse(missing(cex.main) == TRUE, cex.main <- graphics::par()$cex.main, 
        NA)
    ifelse(missing(col.main) == TRUE, col.main <- graphics::par()$col.main, 
        NA)
    ifelse(missing(font.main) == TRUE, font.main <- graphics::par()$font.main, 
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
    if (missing(ecol) == TRUE) {
        ecol <- rep(1, r)
    }
    else {
        ifelse(isTRUE(ecol == 0) == TRUE, ecol <- "transparent", 
            NA)
    }
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
            stop(paste(sQuote("clu"), " should be a vector.", 
                sep = ""))
        if (is.factor(clu) == TRUE) {
            tmpclu <- clu
            for (i in seq_len(nlevels(factor(clu)))) {
                levels(clu) <- c(levels(clu), i)
                clu[which(levels(factor(tmpclu))[i] == clu)] <- i
            }
            rm(i)
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
            rm(tmpclu)
        }
        else {
            NA
        }
        clu <- as.integer(as.vector(clu))
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
        cex <- 3L
    }
    if (isTRUE(length(cex) == 1L) == TRUE) {
        cex <- rep(cex, n)
    }
    else {
        if (is.vector(cex) == FALSE) 
            stop(paste(sQuote("cex"), " should be a vector.", 
                sep = ""))
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
        ifelse(missing(col) == TRUE, vcol <- rep("#FFFFFF", nclu), 
            vcol <- col)
        ifelse(missing(vcol0) == TRUE, vcol0 <- rep("#808080", 
            nclu), NA)
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
    }
    vcol[which(is.na(vcol))] <- graphics::par()$bg
    vcol[which(vcol == 0)] <- graphics::par()$bg
    if (isTRUE(any(pch %in% 21:25)) == TRUE) {
        if (missing(vcol0) == TRUE || isTRUE(vcol0 == 0) == TRUE) {
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
    if (missing(coord) == FALSE) {
        if (isTRUE(nrow(coord) == n) == FALSE) 
            stop(paste("Length of ", sQuote("coord"), " differ from number of vertices in ", 
                sQuote(deparse(substitute(x))), ".", sep = ""))
        flgcrd <- TRUE
        crd <- coord
    }
    else if (missing(coord) == TRUE) {
        flgcrd <- FALSE
        if (isTRUE(conc == TRUE) == FALSE) {
            crd <- frcd(netd, seed = seed, maxiter = maxiter)
        }
        else {
            if (missing(nr) == TRUE) {
                nr <- 1
            }
            else {
                if (isTRUE(dim(netd)[1] < length(nr)) == TRUE) {
                  message(paste("Number of radii in ", sQuote("nr"), 
                    " larger than number of vertices in ", sQuote(deparse(substitute(x))), 
                    ".", sep = ""))
                  nr <- nr[1:dim(netd)[1]]
                }
                else {
                  NA
                }
            }
            crd <- conc(netd, nr, ...)
        }
    }
    fds <- 140L
    ifelse(missing(scl) == TRUE, scl <- rep(1, 2), NA)
    ifelse(isTRUE(length(scl) == 1) == TRUE, scl <- rep(scl, 
        2), scl <- scl[1:2])
    if (isTRUE(max(scl) < 1) == TRUE) {
        fds <- fds - (1/(mean(scl)/30L))
    }
    else if (isTRUE(max(scl) > 1) == TRUE) {
        fds <- fds + (mean(scl) * 20L)
    }
    else {
        NA
    }
    if (missing(vedist) == TRUE) {
        vedist <- 0
    }
    else if (isTRUE(vedist > 10L) == TRUE) {
        vedist <- 10L
    }
    else if (isTRUE(vedist < (-10L)) == TRUE) {
        vedist <- -10L
    }
    else {
        NA
    }
    ifelse(isTRUE(vedist > n) == TRUE, fds <- fds - (n * 10L), 
        fds <- fds - (vedist * 10L))
    ifelse(isTRUE(flgcx == TRUE) == TRUE, fds <- fds - 25L, NA)
    if (missing(rot) == FALSE) {
        crd[, 1:2] <- xyrt(crd[, 1:2], as.numeric(rot))
        crd[, 1:2] <- crd[, 1:2] - min(crd[, 1:2])
        ifelse(isTRUE(n == 2) == TRUE && isTRUE(rot == -90) == 
            TRUE, rot <- -89.9, NA)
    }
    else {
        NA
    }
    ifelse(missing(flip) == FALSE && isTRUE(flip == TRUE) == 
        TRUE, flip <- TRUE, flip <- FALSE)
    if (isTRUE(conc == TRUE) == TRUE) {
        if (isTRUE(flip == TRUE) == TRUE || missing(pos) == TRUE) {
            ifelse(missing(rot) == FALSE, cnt <- 1L, cnt <- 0)
            angl <- vector(length = n)
            for (i in seq_len(n)) {
                ifelse((atan2((crd[i, 2] - cnt), (crd[i, 1] - 
                  cnt)) * (180L/pi)) >= 0, angl[i] <- atan2((crd[i, 
                  2] - cnt), (crd[i, 1] - cnt)) * (180L/pi), 
                  angl[i] <- ((atan2((crd[i, 2] - cnt), (crd[i, 
                    1] - cnt)) * (180L/pi))%%180L) + 180L)
            }
            rm(i)
        }
        if (isTRUE(flgcrd == FALSE) == TRUE && missing(pos) == 
            TRUE) {
            pos <- vector()
            for (i in seq_len(length(angl))) {
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
    if (isTRUE(loops == TRUE) == TRUE) {
        if (isTRUE(conc == TRUE) == TRUE) {
            if (isTRUE(flip == TRUE) == TRUE || missing(pos) == 
                TRUE) {
                ifelse(missing(rot) == FALSE, cnt <- 1L, cnt <- 0)
                angl <- vector(length = n)
                for (i in seq_len(n)) {
                  ifelse((atan2((crd[i, 2] - cnt), (crd[i, 1] - 
                    cnt)) * (180L/pi)) >= 0, angl[i] <- atan2((crd[i, 
                    2] - cnt), (crd[i, 1] - cnt)) * (180L/pi), 
                    angl[i] <- ((atan2((crd[i, 2] - cnt), (crd[i, 
                      1] - cnt)) * (180L/pi))%%180L) + 180L)
                }
                rm(i)
            }
            if (isTRUE(flgcrd == FALSE) == TRUE && missing(pos) == 
                TRUE) {
                pos <- vector()
                for (i in seq_len(length(angl))) {
                  if (isTRUE(65 < angl[i]) == TRUE && isTRUE(115 > 
                    angl[i]) == TRUE) {
                    pos <- append(pos, 3)
                  }
                  else if (isTRUE(115 <= angl[i]) == TRUE && 
                    isTRUE(245 >= angl[i]) == TRUE) {
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
        else {
            if (isTRUE(flip == TRUE) == TRUE) {
                crdf <- cbind(crd[, 1] - mean(crd[, 1]), crd[, 
                  2] - mean(crd[, 2]))
                angl <- vector(length = n)
                for (i in seq_len(n)) {
                  ifelse((atan2((crdf[i, 2]), (crdf[i, 1])) * 
                    (180L/pi)) >= 0, angl[i] <- atan2((crdf[i, 
                    2]), (crdf[i, 1])) * (180L/pi), angl[i] <- ((atan2((crdf[i, 
                    2]), (crdf[i, 1])) * (180L/pi))%%180L) + 
                    180L)
                }
                rm(i)
                if (isTRUE(conc == FALSE) == TRUE) {
                  cmp <- multiplex::comps(netd)$com
                  if (length(cmp) > 1) {
                    rownames(crdf) <- rownames(crd)
                    names(angl) <- rownames(crdf)
                    for (k in seq_len(length(cmp))) {
                      angl[which(names(angl) %in% names(which(crdf[which(names(angl) %in% 
                        cmp[[k]]), 2] < mean(crdf[which(names(angl) %in% 
                        cmp[[k]]), 2]))))] <- 181L
                      angl[which(names(angl) %in% names(which(crdf[which(names(angl) %in% 
                        cmp[[k]]), 2] > mean(crdf[which(names(angl) %in% 
                        cmp[[k]]), 2]))))] <- 179L
                    }
                    rm(k)
                  }
                  else {
                    invisible(NA)
                  }
                  rm(cmp)
                  req <- getElement(suppressWarnings(multiplex::green.rel(x)), 
                    "clu")$R
                  cmp <- vector(mode = "list", length = length(unique(req)))
                  for (i in unique(req)) {
                    cmp[[i]] <- x$st[which(req == i)]
                  }
                  rm(i)
                  if (length(cmp) > 1) {
                    rownames(crdf) <- rownames(crd)
                    names(angl) <- rownames(crdf)
                    for (k in seq_len(length(cmp))) {
                      if (isTRUE(length(cmp[[k]]) > 1L) == TRUE) {
                        angl[which(names(angl) %in% names(which(crdf[which(names(angl) %in% 
                          cmp[[k]]), 2] < mean(crdf[which(names(angl) %in% 
                          cmp[[k]]), 2]))))] <- 181L
                        angl[which(names(angl) %in% names(which(crdf[which(names(angl) %in% 
                          cmp[[k]]), 2] > mean(crdf[which(names(angl) %in% 
                          cmp[[k]]), 2]))))] <- 179L
                      }
                      else {
                        NA
                      }
                    }
                    rm(k)
                  }
                  else {
                    invisible(NA)
                  }
                  rm(cmp)
                }
            }
        }
    }
    if (missing(pos) == TRUE) {
        pos <- 4
    }
    else {
        if (isTRUE(pos < 0L) == TRUE | isTRUE(pos > 4L) == TRUE) 
            stop(paste("Values in ", sQuote("pos"), " must be between 0 and 4.", 
                sep = ""))
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
    if (isTRUE(loops == TRUE) == TRUE && isTRUE(flip == TRUE) == 
        TRUE) {
        ylim <- c(min(nds[, 2]) - (max(cex)/70L), max(nds[, 2]) + 
            (max(cex)/70L))
        xlim <- c(min(nds[, 1]) - (max(cex)/100L), max(nds[, 
            1]) + (max(cex)/100L))
    }
    else if (isTRUE(loops == TRUE) == TRUE) {
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
    ifelse(missing(font.main) == FALSE && isTRUE(font.main %in% 
        names(grDevices::postscriptFonts())) == TRUE, graphics::par(family = font.main), 
        NA)
    suppressWarnings(graphics::plot(nds, type = "n", axes = FALSE, 
        xlab = "", ylab = "", ylim = ylim, xlim = xlim, asp = asp, 
        main = main, cex.main = cex.main, col.main = col.main, 
        ...))
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
                }
                else if (isTRUE(collRecip == TRUE) == FALSE && 
                  isTRUE(tlbs[k] %in% c("recp")) == TRUE) {
                  ifelse(isTRUE(undRecip == TRUE) == TRUE, hd <- 0L, 
                    hd <- hds)
                }
                else {
                  bw <- bwd
                  hd <- hds
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
        netl <- net
        if (isTRUE(swp == TRUE) == TRUE) {
            bdlp <- bd$loop[rev(seq_len(length(bd$loop)))]
            ifelse(isTRUE(flip == TRUE) == FALSE && is.na(dim(net)[3]) == 
                TRUE, NA, netl <- net[, , rev(seq_len(dim(net)[3]))])
            Ltl <- Lt[rev(seq_len(length(Lt)))]
        }
        else {
            bdlp <- bd$loop
            Ltl <- Lt
        }
        ifelse(isTRUE(length(lty) == 1) == TRUE, Ltl <- rep(lty, 
            length(bdlp)), NA)
        if (missing(vldist) == FALSE) {
            if (isTRUE(abs(vldist) > cex[1]) == TRUE) {
                ifelse(isTRUE(vldist < 0) == TRUE, vldist <- cex[1] * 
                  -1L, vldist <- cex[1] * 1L)
            }
        }
        else {
            NA
        }
        ndss <- nds
        ndss[, 1] <- ndss[, 1] * scl[1]
        ndss[, 2] <- ndss[, 2] * scl[2]
        if (isTRUE(z == 1L) == TRUE) {
            lp <- as.vector(which(diag(net) > 0))
            if (isTRUE(length(lp) > 0) == TRUE) {
                if ((missing(coord) == FALSE || isTRUE(conc == 
                  TRUE) == TRUE) && isTRUE(flip == TRUE) == TRUE) {
                  flipl <- vector(length = n)
                  for (l in seq_len(length(angl))) {
                    ifelse(isTRUE(180 < angl[lp[l]]) == TRUE, 
                      flipl[lp[l]] <- TRUE, NA)
                  }
                  rm(l)
                }
                else {
                  flipl <- vector(length = length(bdlp))
                }
                dz <- (rng(1) + abs(min(rng(1))))/(10L)
                for (i in seq_len(length(lp))) {
                  if (isTRUE(n == 1L) == TRUE) {
                    dcx <- (cex[lp[i]]^2) * 0.00025
                    rds <- abs((cex[lp[i]]^2) * 0.00015)
                  }
                  else if (isTRUE(n < 3) == TRUE) {
                    dcx <- cex[lp[i]]/110L
                    rds <- abs((cex[lp[i]] * 0.007) - dz)
                  }
                  else {
                    dcx <- cex[lp[i]]/100L
                    rds <- abs((cex[lp[i]] * 0.0075) - dz)
                  }
                  ifelse(isTRUE(mirrorY == TRUE) == TRUE && isTRUE(flip == 
                    TRUE) == TRUE && isTRUE(n > 1) == TRUE, flipl <- !flipl, 
                    NA)
                  if (isTRUE(flipl[lp[i]] == TRUE) == TRUE) {
                    dcx <- dcx * (-1e-17)
                    if (missing(vldist) == FALSE) {
                      dcx <- dcx - (vldist * 0.008)
                    }
                    else {
                      ifelse(isTRUE(cex[lp[i]] > 5) == TRUE, 
                        dcx <- dcx - ((cex[lp[i]]/2L) * 0.002), 
                        NA)
                    }
                  }
                  else {
                    ifelse(missing(vldist) == FALSE, dcx <- dcx + 
                      (vldist * 0.01), NA)
                  }
                  hc(ndss[lp[i], 1], ndss[lp[i], 2] + (dcx), 
                    rds, col = grDevices::adjustcolor(vecol, 
                      alpha = alfa), lty = Ltl, lwd = lwd, flip = flipl[lp[i]])
                }
                rm(i)
            }
            else {
                NA
            }
        }
        else if (isTRUE(z > 1) == TRUE) {
            if (missing(bwd2) == TRUE) {
                bwd2 <- 1L
            }
            else if (missing(bwd2) == FALSE) {
                if (isTRUE(bwd2 < 1) == TRUE && isTRUE(bwd2 > 
                  0) == TRUE) {
                  NA
                }
                else if (isTRUE(bwd2 < 0) == TRUE) {
                  bwd2 <- 0
                }
                else if (isTRUE(bwd2 > 4L) == TRUE) {
                  ifelse(isTRUE(z > 2) == TRUE, bwd2 <- 2L, bwd2 <- 4L)
                }
            }
            if (isTRUE(bwd2 == 0) == TRUE) {
                dz <- rep(0, z)
            }
            else {
                if (isTRUE(max(cex) < 1.5) == TRUE) {
                  dz <- c(0, rep((cex/1000L), (z - 1L)))
                }
                else {
                  dz <- (bwd2) * (rng(z) + abs(min(rng(z))))/(20L)
                }
            }
            for (k in seq_len(length(bdlp))) {
                lp <- as.numeric(unique(multiplex::dhc(bdlp)[k][[1]]))
                if (isTRUE(length(lp) > 0) == TRUE) {
                  for (i in seq_len(length(lp))) {
                    if (isTRUE(flip == TRUE) == TRUE) {
                      flipl <- vector(length = n)
                      for (l in seq_len(length(angl))) {
                        ifelse(isTRUE(180 < angl[lp[i]]) == TRUE, 
                          flipl[k] <- TRUE, NA)
                      }
                      rm(l)
                    }
                    else {
                      flipl <- vector(length = length(bdlp))
                    }
                    ifelse(isTRUE(mirrorY == TRUE) == TRUE && 
                      isTRUE(flip == TRUE) == TRUE, flipl <- !flipl, 
                      NA)
                    if (isTRUE(n == 1L) == TRUE) {
                      dcx <- (cex[lp[i]]^2) * 0.00025
                      rds <- abs((cex[lp[i]]^2) * 0.00015 - (dz))
                    }
                    else if (isTRUE(n < 3) == TRUE) {
                      dcx <- cex[lp[i]]/110L
                      rds <- abs((cex[lp[i]] * 0.007) - dz[k])
                    }
                    else {
                      dcx <- cex[lp[i]]/100L
                      rds <- abs((cex[lp[i]] * 0.0075) - dz[k])
                    }
                    if (isTRUE(flipl[k] == TRUE) == TRUE) {
                      if (isTRUE(bwd2 > 0) == TRUE) {
                        if (isTRUE(k == 1L) == TRUE) {
                          dcx <- dcx - (cex[lp[i]] * 0.011)
                        }
                        else {
                          prds <- rds
                          if (isTRUE(cex[lp[i]] < 3L) == TRUE) {
                            rds <- rds * 0.8
                          }
                          else {
                            rds <- rds * 0.85
                          }
                          dcx <- dcx - 2L * (k/1000L) - (cex[lp[i]] * 
                            0.012)
                        }
                      }
                      else {
                        dcx <- dcx - (cex[lp[i]] * 0.0125)
                      }
                      if (missing(vldist) == FALSE) {
                        if (isTRUE(vldist < 0) == TRUE) {
                          dcx <- dcx - (vldist * 0.004)
                        }
                        else {
                          dcx <- dcx - (vldist * 0.01)
                        }
                      }
                    }
                    else if (isTRUE(flipl[k] == FALSE) == TRUE) {
                      if (isTRUE(k > 1L) == TRUE && isTRUE(bwd2 > 
                        0) == TRUE) {
                        if (isTRUE(cex[lp[i]] < 3L) == TRUE) {
                          rds <- rds * 0.8
                        }
                        else {
                          rds <- rds * 0.85
                        }
                      }
                      else {
                        NA
                      }
                      if (missing(vldist) == FALSE) {
                        if (isTRUE(vldist < 0) == TRUE) {
                          dcx <- dcx + (vldist * 0.004)
                        }
                        else {
                          dcx <- dcx + (vldist * 0.01)
                        }
                      }
                    }
                    hc(ndss[lp[i], 1], ndss[lp[i], 2] + (dcx), 
                      rds, col = grDevices::adjustcolor(vecol[k], 
                        alpha = alfa), lty = Ltl[k], lwd = lwd[k], 
                      flip = flipl[k])
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
                      function(z) bquote(italic(.(z))))), cex = fsize, 
                      adj = 0.5, col = fcol)
                  }
                  else if (isTRUE(fstyle == "bold") == TRUE) {
                    graphics::text(ndss, labels = as.expression(lapply(lbs, 
                      function(z) bquote(bold(.(z))))), cex = fsize, 
                      adj = 0.5, col = fcol)
                  }
                  else if (isTRUE(fstyle == "bolditalic") == 
                    TRUE) {
                    graphics::text(ndss, labels = as.expression(lapply(lbs, 
                      function(z) bquote(bolditalic(.(z))))), 
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
                      function(z) bquote(italic(.(z))))), cex = fsize, 
                      pos = pos, col = fcol, offset = (cex/4L), 
                      adj = c(0.5, 1))
                  }
                  else if (isTRUE(fstyle == "bold") == TRUE) {
                    graphics::text(ndss, as.expression(lapply(lbs, 
                      function(z) bquote(bold(.(z))))), cex = fsize, 
                      pos = pos, col = fcol, offset = (cex/4L), 
                      adj = c(0.5, 1))
                  }
                  else if (isTRUE(fstyle == "bolditalic") == 
                    TRUE) {
                    graphics::text(ndss, as.expression(lapply(lbs, 
                      function(z) bquote(bolditalic(.(z))))), 
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
                    function(z) bquote(italic(.(z))))), cex = fsize, 
                    pos = pos, col = fcol[1], offset = (cex/4L), 
                    adj = c(0.5, 1))
                }
                else if (isTRUE(fstyle == "bold") == TRUE) {
                  graphics::text(ndss, as.expression(lapply(lbs, 
                    function(z) bquote(bold(.(z))))), cex = fsize, 
                    pos = pos, col = fcol[1], offset = (cex/4L), 
                    adj = c(0.5, 1))
                }
                else if (isTRUE(fstyle == "bolditalic") == TRUE) {
                  graphics::text(ndss, as.expression(lapply(lbs, 
                    function(z) bquote(bolditalic(.(z))))), cex = fsize, 
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
                      function(z) bquote(italic(.(z))))), cex = fsize, 
                      adj = 0.5, col = fcol)
                  }
                  else if (isTRUE(fstyle == "bold") == TRUE) {
                    graphics::text(ndss, labels = as.expression(lapply(lbs, 
                      function(z) bquote(bold(.(z))))), cex = fsize, 
                      adj = 0.5, col = fcol)
                  }
                  else if (isTRUE(fstyle == "bolditalic") == 
                    TRUE) {
                    graphics::text(ndss, labels = as.expression(lapply(lbs, 
                      function(z) bquote(bolditalic(.(z))))), 
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
                      function(z) bquote(italic(.(z))))), cex = fsize, 
                      pos = pos[1], col = fcol, offset = (cex/4L), 
                      adj = c(0.5, 1))
                  }
                  else if (isTRUE(fstyle == "bold") == TRUE) {
                    graphics::text(ndss, as.expression(lapply(lbs, 
                      function(z) bquote(bold(.(z))))), cex = fsize, 
                      pos = pos[1], col = fcol, offset = (cex/4L), 
                      adj = c(0.5, 1))
                  }
                  else if (isTRUE(fstyle == "bolditalic") == 
                    TRUE) {
                    graphics::text(ndss, as.expression(lapply(lbs, 
                      function(z) bquote(bolditalic(.(z))))), 
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

