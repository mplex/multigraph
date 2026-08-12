## (GPL-3)  Antonio Rivero Ostoic, multiplex@post.com
mlgraph <-
function (net, layout = c("circ", "force", "stress", "rand", 
    "conc", "bip"), main = NULL, seed = NULL, maxiter = 100, 
    directed = TRUE, alpha = c(1, 1, 1), scope, collRecip, lbs, 
    showLbs, showAtts, cex.main, col.main, font.main, coord, 
    clu, cex, lwd, pch, lty, bwd, bwd2, att, bg, mar, pos, asp, 
    ecol, vcol, vcol0, col, lbat, loops, rot, mirrorX, mirrorY, 
    mirrorD, mirrorL, mirrorV, mirrorH, hds, scl, vedist, vldist, 
    flip, ffamily, fstyle, fsize, fcol, valued, modes, elv, lng, 
    nr, new, jitter, ...) 
{
    mlv <- net
    if (isTRUE("Multilevel" %in% attr(mlv, "class")) == TRUE) {
        ifelse(all(mlv$lbs$cdm %in% mlv$lbs$dm) == FALSE, invisible(NA), 
            warning(paste("Equal domain and codomain vertex labels in ", 
                sQuote(deparse(substitute(net))), ".", sep = ""), 
                call. = FALSE))
        if ("bpn" %in% attr(mlv, "class") || "cn2" %in% attr(mlv, 
            "class")) {
            flgcn <- FALSE
            net <- mlv$mlnet[, , which(mlv$modes == "1M")]
            met <- mlv$mlnet[, , which(mlv$modes == "2M")]
        }
        else if ("cn" %in% attr(mlv, "class")) {
            flgcn <- TRUE
            met <- multiplex::dichot(mlv$mlnet, c = max(mlv$mlnet) + 
                1L)
            net <- mlv$mlnet
        }
        else {
            vcn <- vector()
            for (i in seq_len(length(mlv$mlnet))) {
                vcn <- append(vcn, dimnames(mlv$mlnet[[i]])[[1]])
                vcn <- append(vcn, dimnames(mlv$mlnet[[i]])[[2]])
            }
            rm(i)
            bmlbs <- unique(vcn)
            bmat <- multiplex::transf(mlv$mlnet[[1]], type = "toarray2", 
                lbs = bmlbs)
            for (k in seq(from = 2, to = length(mlv$mlnet))) {
                bmat <- multiplex::zbnd(bmat, multiplex::transf(mlv$mlnet[[k]], 
                  type = "toarray2", lbs = bmlbs))
            }
            rm(k)
            dimnames(bmat)[[1]] <- dimnames(bmat)[[2]] <- bmlbs
            dimnames(bmat)[[3]] <- attr(mlv$mlnet, "names")
            for (i in which(mlv$modes == "2M")) {
                bmat[, , i] <- bmat[, , i] + t(bmat[, , i])
            }
            rm(i)
            net <- bmat[, , which(mlv$modes == "1M")]
            met <- bmat[, , which(mlv$modes == "2M")]
        }
    }
    else {
        if ((missing(modes) == FALSE && is.vector(modes) == TRUE) && 
            (missing(clu) == FALSE && is.vector(clu) == TRUE)) {
            ifelse(is.numeric(modes) == TRUE, modes <- paste(modes, 
                "M", sep = ""), NA)
            mlv <- list(mlnet = net, lbs = list(dm = dimnames(net)[[1]][which(clu == 
                1)], cdm = dimnames(net)[[1]][which(clu == 2)]), 
                modes = modes)
            net <- mlv$mlnet[, , which(mlv$modes == "1M")]
            met <- mlv$mlnet[, , which(mlv$modes == "2M")]
        }
        else {
            stop(paste(sQuote(deparse(substitute(net))), " is not a ", 
                dQuote("Multilevel"), " class object or array with clustering information.", 
                sep = ""))
        }
    }
    ifelse(isTRUE(dim(net)[3] == 1) == TRUE, net <- net[, , 1], 
        NA)
    ifelse(missing(valued) == FALSE && isTRUE(valued == TRUE) == 
        TRUE, valued <- TRUE, valued <- FALSE)
    ifelse(missing(loops) == FALSE && isTRUE(loops == FALSE) == 
        TRUE, loops <- FALSE, loops <- TRUE)
    ifelse(missing(collRecip) == FALSE && isTRUE(collRecip == 
        FALSE) == TRUE, collRecip <- FALSE, collRecip <- TRUE)
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
    ifelse(missing(showAtts) == FALSE && isTRUE(showAtts == FALSE) == 
        TRUE, showAtts <- FALSE, showAtts <- TRUE)
    ifelse(isTRUE(directed == FALSE) == TRUE, directed <- FALSE, 
        NA)
    if (missing(scope) == FALSE) {
        if (isTRUE(is.list(scope) == TRUE) == FALSE) 
            stop(paste(sQuote("scope"), " must be a list or a vector of lists.", 
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
    if (isTRUE("Multilevel" %in% attr(mlv, "class")) == TRUE) {
        if (missing(lbs) == TRUE) {
            if (missing(pch) == TRUE && all(mlv$lbs$cdm %in% 
                mlv$lbs$dm) == FALSE) {
                pch <- c(rep(21, length(mlv$lbs$dm)), rep(22, 
                  length(which(!(mlv$lbs$cdm %in% mlv$lbs$dm)))))
            }
            else {
                if (missing(pch) == FALSE) {
                  ifelse(isTRUE(length(pch) == 2) == TRUE, pch <- c(rep(pch[1], 
                    length(mlv$lbs$dm)), rep(pch[2], length(mlv$lbs$cdm))), 
                    pch <- pch[seq_len(length(mlv$lbs$dm) + length(mlv$lbs$cdm))])
                }
                else {
                  pch <- 21
                }
            }
        }
        else if (missing(lbs) == FALSE && is.list(lbs) == TRUE) {
            if (missing(pch) == TRUE) {
                pch <- rep(21, dim(mlv$mlnet)[1])
                pch[which(dimnames(mlv$mlnet)[[1]] %in% lbs[[1]])] <- 22
            }
            else {
                if (isTRUE(length(pch) == 2) == TRUE) {
                  pchtmp <- pch
                  pch <- rep(pchtmp[1], dim(mlv$mlnet)[1])
                  pch[which(dimnames(mlv$mlnet)[[1]] %in% lbs[[1]])] <- pchtmp[2]
                  rm(pchtmp)
                }
            }
        }
        else {
            NA
        }
    }
    ifelse(missing(lng) == TRUE, lng <- 50, NA)
    ifelse(isTRUE(lng <= 2) == TRUE, lng <- 3L, NA)
    if (match.arg(layout) == "bip") {
        ifelse(missing(asp) == TRUE, asp <- 2L, asp <- asp[1] * 
            2L)
    }
    else {
        ifelse(missing(asp) == TRUE, asp <- 1L, NA)
    }
    ifelse(missing(fcol) == TRUE, fcol <- c(1, 1), fcol <- rep(fcol, 
        2))
    if (missing(vcol) == TRUE) {
        ifelse(missing(col) == TRUE, vcol <- c(rep("#C0C0C0", 
            length(mlv$lbs$dm)), rep("#408080", length(which(!(mlv$lbs$cdm %in% 
            mlv$lbs$dm))))), vcol <- col)
        ifelse(missing(vcol0) == TRUE, vcol0 <- rep("#000000", 
            2), NA)
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
    ifelse(missing(bwd) == TRUE, bwd <- 1, NA)
    ifelse(isTRUE(bwd < 0L) == TRUE, bwd <- 0L, NA)
    obg <- graphics::par()$bg
    if (missing(bg) == TRUE) {
        bg <- graphics::par()$bg
    }
    else {
        obg <- graphics::par()$bg
    }
    ifelse(missing(cex.main) == TRUE, cex.main <- graphics::par()$cex.main, 
        NA)
    ifelse(missing(col.main) == TRUE, col.main <- graphics::par()$col.main, 
        NA)
    ifelse(missing(font.main) == TRUE, font.main <- graphics::par()$font.main, 
        NA)
    ifelse(missing(rot) == TRUE, rot <- 0, rot <- rot[1] * -1)
    ifelse(missing(vedist) == TRUE, vedist <- 0, NA)
    ifelse(isTRUE(vedist > 1L) == TRUE, vedist <- 1L, NA)
    ifelse(missing(scl) == TRUE, scl <- rep(1, 2), NA)
    ifelse(isTRUE(length(scl) == 1) == TRUE, scl <- rep(scl, 
        2), scl <- scl[1:2])
    if (!(missing(hds))) {
        if (isTRUE(hds > 1L) == TRUE) {
            hds <- (hds/1.5)
        }
        else if (isTRUE(hds < 1L) == TRUE) {
            hds <- (hds/(hds + 0.15))
        }
        else if (isTRUE(hds == 0L) == TRUE) {
            hds <- 0.01
        }
    }
    else {
        hds <- 1L * min(scl)
    }
    n <- dim(net)[1]
    ifelse(isTRUE(is.na(dim(net)[3]) == TRUE) == TRUE, z <- 1L, 
        z <- dim(net)[3])
    if (missing(lbs) == TRUE) {
        ifelse(is.null(dimnames(net)[[1]]) == TRUE, lbs <- as.character(seq_len(dim(net)[1])), 
            lbs <- dimnames(net)[[1]])
    }
    else {
        if (is.list(lbs) == FALSE) {
            NA
        }
        else {
            if ("cn2" %in% attr(mlv, "class")) {
                ifelse(isTRUE(length(lbs[[1]]) != length(mlv$lbs$dm)) == 
                  TRUE, lbs <- mlv$lbs$dm, NA)
            }
            else if (isTRUE(length(lbs) > 1) == TRUE) {
                if (all(mlv$lbs$dm %in% mlv$lbs$cdm) == TRUE && 
                  all(mlv$lbs$cdm %in% mlv$lbs$dm) == TRUE) {
                  ifelse(isTRUE(length(lbs[[1]]) != length(mlv$lbs$dm)) == 
                    TRUE, lbs[[1]] <- mlv$lbs$dm, NA)
                  ifelse(isTRUE(length(lbs[[2]]) != length(mlv$lbs$cdm)) == 
                    TRUE, lbs[[2]] <- mlv$lbs$cdm, NA)
                }
                else {
                  NA
                }
            }
            else {
                NA
            }
            lbs <- unlist(lbs, use.names = FALSE)
        }
    }
    if (isTRUE(n != length(unique(lbs))) == TRUE) 
        stop(paste(sQuote("lbs"), " differs number of vertices in ", 
            sQuote(deparse(substitute(net))), ".", sep = ""), 
            call. = FALSE)
    if (missing(att) == FALSE && is.array(att) == TRUE) {
        if (isTRUE(n != dim(att)[1]) == TRUE) {
            warning(paste("No attributes are shown due ", sQuote(deparse(substitute(net))), 
                " and ", sQuote(deparse(substitute(att))), " differ in dimensions.", 
                sep = ""), call. = FALSE)
            showAtts <- FALSE
        }
    }
    netd <- multiplex::dichot(net, c = 1L)
    if (isTRUE(directed == FALSE) == TRUE && isTRUE(collRecip == 
        TRUE) == TRUE && isTRUE(valued == TRUE) == TRUE) {
        if (isTRUE(z == 1L) == TRUE) {
            net <- net + t(net)
        }
        else {
            for (i in seq_len(z)) {
                net[, , i] <- net[, , i] + t(net[, , i])
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
            rcp[lower.tri(rcp, diag = TRUE)] <- 0L
        }
        else {
            nt <- array(0L, dim = c(n, n, z))
            dimnames(nt)[[1]] <- dimnames(nt)[[2]] <- unique(lbs)
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
    if (isTRUE(collRecip == TRUE) == TRUE) {
        bd <- multiplex::bundles(ucnet, loops = loops, lb2lb = FALSE, 
            collapse = FALSE)
        ifelse(isTRUE(directed == TRUE) == FALSE, NA, bd$recp <- multiplex::bundles(netd, 
            loops = loops, lb2lb = FALSE, collapse = FALSE)$recp)
    }
    else {
        bd <- multiplex::bundles(netd, loops = loops, lb2lb = FALSE, 
            collapse = FALSE)
    }
    ifelse(isTRUE(z == 1L) == TRUE, r <- 1L, r <- length(bd[[1]]))
    ifelse(isTRUE(sum(net) == 0) == TRUE && isTRUE(loops == TRUE) == 
        TRUE, bd$loop <- character(0), NA)
    bds <- multiplex::summaryBundles(bd, byties = TRUE)
    m <- dim(met)[1]
    ifelse(is.na(dim(met)[3]) == TRUE, zz <- 1L, zz <- dim(met)[3])
    metd <- multiplex::dichot(met, c = 1L)
    if (isTRUE(directed == FALSE) == TRUE && isTRUE(collRecip == 
        TRUE) == TRUE) {
        if (isTRUE(zz == 1L) == TRUE) {
            mt <- metd + t(metd)
            rcpm <- multiplex::dichot(mt, c = 1L)
            rcpm[lower.tri(rcpm, diag = TRUE)] <- 0L
        }
        else {
            mt <- array(0L, dim = c(m, m, zz))
            dimnames(mt)[[1]] <- dimnames(mt)[[2]] <- unique(lbs)
            dimnames(mt)[[3]] <- dimnames(met)[[3]]
            for (i in seq_len(zz)) {
                mt[, , i] <- metd[, , i] + t(metd[, , i])
            }
            rm(i)
            ifelse(isTRUE(max(mt) > 1L) == TRUE, rcpm <- multiplex::dichot(mt, 
                c = 2L), rcpm <- mt)
            for (i in seq_len(zz)) {
                rcpm[, , i][lower.tri(rcpm[, , i], diag = TRUE)] <- 0L
            }
            rm(i)
        }
    }
    else {
        NA
    }
    ucmet <- metd
    if (isTRUE(collRecip == TRUE) == TRUE) {
        bdm <- multiplex::bundles(ucmet, loops = loops, lb2lb = FALSE, 
            collapse = FALSE)
        ifelse(isTRUE(directed == TRUE) == FALSE, NA, bdm$recp <- multiplex::bundles(metd, 
            loops = loops, lb2lb = FALSE, collapse = FALSE)$recp)
    }
    else {
        bdm <- multiplex::bundles(metd, loops = loops, lb2lb = FALSE, 
            collapse = FALSE)
    }
    bdsm <- multiplex::summaryBundles(bdm, byties = TRUE)
    ifelse(isTRUE(zz == 1L) == TRUE, rr <- 1L, rr <- length(bdm[[1]]))
    ifelse(missing(ecol) == TRUE, ecol <- grDevices::gray.colors(r + 
        rr, start = 0.1, end = 0.8), NA)
    ifelse(missing(ecol) == FALSE && isTRUE(length(ecol) == 2) == 
        TRUE, ecol <- c(rep(ecol[1], length(which(mlv$modes == 
        "1M"))), rep(ecol[2], length(which(mlv$modes == "2M")))), 
        NA)
    ifelse(isTRUE(ecol == 0) == TRUE, ecol <- c(rep("#FFFFFF", 
        length(which(mlv$modes == "1M"))), rep("#FFFFFF", length(which(mlv$modes == 
        "2M")))), NA)
    ecol[which(ecol == 0)] <- "#FFFFFF"
    if (missing(lty) == TRUE) {
        lty <- c((seq_len(r) + 1L), rep(1, rr))
        ifelse(isTRUE(valued == TRUE) == TRUE, lty <- rep(1, 
            r + rr), lty <- seq_len(r + rr))
    }
    else if (missing(lty) == FALSE) {
        ifelse(isTRUE(valued == TRUE) == TRUE, lty <- abs(unlist(lty)), 
            lty <- rep(abs(unlist(lty)), (r + rr)))
    }
    ifelse(missing(lwd) == TRUE, lwd <- 1, NA)
    if (missing(hds) == FALSE) {
        if (isTRUE(hds > 2L) == TRUE) {
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
        hds <- 1L * scl
    }
    ifelse(isTRUE(dim(net)[1] > 8) == TRUE || isTRUE(valued == 
        TRUE) == TRUE || isTRUE(min(lwd) >= 4) == TRUE, hds <- hds * 
        0.75, NA)
    if (isTRUE((z + zz) == 1L) == TRUE) {
        Lt <- lty[1]
        vecol <- ecol[1]
    }
    else {
        ifelse(isTRUE("Multilevel" %in% attr(mlv, "class")) == 
            TRUE && isTRUE(flgcn == FALSE) == TRUE, lty <- rev(lty), 
            NA)
        if (isTRUE(zz > 1) == TRUE || isTRUE(zz == 1 && z == 
            1) == TRUE) {
            ifelse(isTRUE(length(ecol) == 1L) == TRUE, vecol <- rep(ecol, 
                z + zz), vecol <- rep(ecol, z + zz)[seq_len(z + 
                zz)])
            ifelse(isTRUE(length(lty) == 1L) == TRUE, Lt <- rep(lty, 
                r + rr), Lt <- rep(lty, r + rr)[seq_len(r + rr)])
        }
        else if (isTRUE(zz == 1) == TRUE) {
            ifelse(isTRUE(length(ecol) == 1L) == TRUE, vecol <- rep(ecol, 
                z), vecol <- rep(ecol, z)[seq_len(z)])
            ifelse(isTRUE(length(lty) == 1L) == TRUE, Lt <- seq_len(r), 
                Lt <- c(lty[1:r], lty[(r + 1):length(lty)]))
        }
        if (isTRUE(length(lty) == length(Lt)) == FALSE) {
            Ltc <- seq_along(vecol)
        }
        else {
            if (isTRUE(seq(lty) == lty) == TRUE) {
                Ltc <- Lt
            }
            else {
                if (isTRUE(zz > 1) == TRUE) {
                  Ltc <- seq_len(r + rr)
                }
                else if (isTRUE(zz == 1) == TRUE) {
                  Ltc <- seq_len(r)
                }
            }
        }
    }
    vltz <- Lt
    if (missing(clu) == FALSE) {
        if (any(unlist(lapply(clu, is.null))) == TRUE) {
            ifelse(unlist(lapply(clu, is.null))[1] == TRUE, clu[[1]] <- rep(1L, 
                nrow(net)), NA)
            ifelse(unlist(lapply(clu, is.null))[2] == TRUE, clu[[2]] <- rep(1L, 
                ncol(net)), NA)
        }
        else {
            NA
        }
        if (is.list(clu) == FALSE && (is.vector(clu) == TRUE || 
            is.factor(clu) == TRUE)) {
            ifelse(isTRUE(length(clu) == length(unlist(mlv$lbs))) == 
                TRUE, NA, clu <- list(clu, rep(1L, length(mlv$lbs$cdm))))
        }
        else if (is.list(clu) == TRUE) {
            if (any(unlist(lapply(clu, is.factor))) == TRUE) {
                uclu <- unlist(clu)
                uclu[(length(clu[[1]]) + 1):length(uclu)] <- max(uclu) + 
                  1L
                clu <- uclu
                rm(uclu)
            }
            else {
                clu <- unlist(clu)
            }
        }
        else {
            NA
        }
        ifelse(is.factor(clu[[1]]) == TRUE, uact <- levels(clu[[1]]), 
            uact <- unique(clu[[1]]))
        ifelse(is.factor(clu[[2]]) == TRUE, uevt <- levels(clu[[2]]), 
            uevt <- unique(clu[[2]]))
        if (is.character(uact) == TRUE) {
            tmpa <- as.vector(clu[[1]])
            for (i in seq_len(length(uact))) {
                tmpa[which(tmpa == uact[i])] <- i
            }
            rm(i)
            clu[[1]] <- as.numeric(tmpa)
            rm(tmpa)
        }
        if (is.character(uevt) == TRUE) {
            tmpe <- as.vector(clu[[2]])
            for (i in seq_len(length(uevt))) {
                tmpe[which(tmpe == uevt[i])] <- i
            }
            rm(i)
            clu[[2]] <- as.numeric(tmpe)
            rm(tmpe)
        }
        if (is.list(clu) == TRUE) {
            clutmp <- clu
            if (any(clutmp[[2]] %in% clutmp[[1]]) == TRUE) {
                k <- 0L
                tmp2 <- clutmp[[2]]
                while (any(tmp2 %in% clutmp[[1]]) == TRUE) {
                  tmp2 <- replace(tmp2, which(tmp2 == min(tmp2)), 
                    (max(clutmp[[1]]) + k))
                  k <- k + 1L
                }
                clutmp[[2]] <- tmp2
                rm(tmp2)
                clu <- as.vector(unlist(clutmp))
            }
            else if (any(clu[[2]] %in% clu[[1]]) == TRUE) {
                k <- 0L
                tmp2 <- clu[[2]]
                while (any(tmp2 %in% clu[[1]]) == TRUE) {
                  tmp2 <- replace(tmp2, which(tmp2 == min(tmp2)), 
                    (max(clu[[1]]) + k))
                  k <- k + 1L
                }
                clu[[2]] <- tmp2
                rm(tmp2)
                clu <- as.vector(unlist(clu))
            }
        }
        else {
            NA
        }
        nclu <- nlevels(factor(clu))
        flgclu <- TRUE
    }
    else {
        flgclu <- FALSE
        if (missing(vcol) == FALSE && (is.vector(vcol) == TRUE || 
            is.list(vcol) == TRUE)) {
            nclu <- c(rep(vcol[1], length(mlv$lbs$dm)), rep(vcol[2], 
                length(mlv$lbs$cdm)))
        }
        else {
            nclu <- 1L
        }
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
            stop(paste(sQuote("cex"), " must be a vector.", sep = ""), 
                call. = FALSE)
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
    else if (isTRUE(flgcx == FALSE) == TRUE && isTRUE(valued == 
        TRUE) == TRUE) {
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
        fsize <- fsize/10L
    }
    ifelse(isTRUE(valued == FALSE) == TRUE && isTRUE(bwd > 1L) == 
        TRUE, bwd <- 1L, NA)
    ifelse(isTRUE(max(cex) < 2) == TRUE, NA, bwd <- bwd * 0.75)
    if (isTRUE("Multilevel" %in% attr(mlv, "class")) == FALSE) {
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
    }
    if (isTRUE(length(vcol) == 1L) == TRUE) {
        vcol <- rep(vcol, n)
    }
    else if (isTRUE(flgclu == FALSE) == TRUE) {
        vcol <- c(rep(vcol[1], length(mlv$lbs$dm)), rep(vcol[2], 
            length(mlv$lbs$cdm)))
    }
    else if (isTRUE(flgclu == TRUE) == TRUE) {
        if (isTRUE(length(vcol) == nclu) == TRUE) {
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
        else {
            if (isTRUE(length(vcol) != nclu) == TRUE) {
                vvcol <- vector(length = length(clu))
                for (k in unique(clu)) vvcol[which(clu == k)] <- vcol[k]
                vvcol[(length(mlv$lbs$dm) + 1L):length(vvcol)] <- "#000000"
                ifelse(any(is.na(vvcol)) == TRUE, vvcol[which(is.na(vvcol))] <- "#000000", 
                  NA)
                vcol <- vvcol
                rm(k, vvcol)
            }
            else {
                NA
            }
        }
    }
    vcol[which(is.na(vcol))] <- graphics::par()$bg
    vcol[which(vcol == 0)] <- graphics::par()$bg
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
    if (missing(coord) == FALSE) {
        if (isTRUE(nrow(coord) == n) == FALSE) 
            stop(paste(sQuote("coord"), " does not match number of vertices in ", 
                sQuote(deparse(substitute(net))), ".", sep = ""))
        flgcrd <- TRUE
        crd <- coord
    }
    else if (missing(coord) == TRUE) {
        flgcrd <- FALSE
        switch(match.arg(layout), force = {
            crd <- frcd(zbnd(netd, met), seed = seed, maxiter = maxiter)
        }, circ = {
            crd <- data.frame(X = sin(2L * pi * ((0:(n - 1L))/n)), 
                Y = cos(2L * pi * ((0:(n - 1L))/n)))
        }, stress = {
            crd <- stsm(zbnd(netd, met), seed = seed, maxiter = maxiter, 
                ...)
            ifelse(missing(jitter) == FALSE, crd <- apply(crd, 
                2, function(z) {
                  jitter(z, amount = (min(10, jitter))/100L)
                }), NA)
        }, rand = {
            set.seed(seed)
            crd <- data.frame(X = round(stats::runif(n) * 1L, 
                5), Y = round(stats::runif(n) * 1L, 5))
        }, conc = {
            netdm <- zbnd(netd, met)
            if (missing(nr) == TRUE) {
                nr <- 1
            } else {
                if (isTRUE(dim(netdm)[1] < length(nr)) == TRUE) {
                  warning(paste(sQuote("nr"), " is larger than network dimensions.", 
                    sep = ""), call. = FALSE)
                  nr <- nr[1:dim(netdm)[1]]
                } else {
                  NA
                }
            }
            crd <- conc(netdm, nr = nr, ...)
        }, bip = {
            if ("bpn" %in% attr(mlv, "class") && identical(sort(mlv$lbs$dm), 
                sort(mlv$lbs$cdm)) == FALSE) {
                act <- nrm(rng(length(mlv$lbs$dm)))
                evt <- nrm(rng(length(mlv$lbs$cdm)))
                Act <- cbind(rep(0, length(mlv$lbs$dm)), act)
                Evt <- cbind(rep(1, length(mlv$lbs$cdm)), evt)
            } else {
                act <- nrm(rng(length(mlv$lbs$dm)/2L))
                evt <- nrm(rng(length(mlv$lbs$dm)/2L))
                Act <- cbind(rep(0, length(mlv$lbs$dm)/2L), act)
                Evt <- cbind(rep(1, length(mlv$lbs$dm)/2L), evt)
            }
            crd <- rbind(Act, Evt)
            crd[which(is.nan(crd))] <- 0.5
            crd[, 2] <- crd[, 2] * cos(pi) - crd[, 1] * sin(pi)
        })
    }
    ifelse(isTRUE(directed == TRUE) == TRUE, fds <- 130L, fds <- 140L)
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
    ifelse(isTRUE(valued == TRUE) == TRUE, fds <- fds + 25L, 
        NA)
    ifelse((missing(rot) == TRUE || (missing(rot) == FALSE && 
        isTRUE(as.numeric(rot) == 0) == TRUE)) && isTRUE(loops == 
        TRUE) == TRUE && match.arg(layout) == "conc", rot <- 0.01, 
        NA)
    if (missing(rot) == FALSE) {
        crd[, 1:2] <- xyrt(crd[, 1:2], as.numeric(rot))
        crd[, 1:2] <- crd[, 1:2] - min(crd[, 1:2])
    }
    else {
        NA
    }
    ifelse(missing(mirrorX) == FALSE && isTRUE(mirrorX == TRUE) == 
        TRUE || missing(mirrorV) == FALSE && isTRUE(mirrorV == 
        TRUE) == TRUE, crd[, 1] <- crd[, 1] * cos(pi) - crd[, 
        2] * sin(pi), mirrorX <- FALSE)
    ifelse(missing(mirrorY) == FALSE && isTRUE(mirrorY == TRUE) == 
        TRUE || missing(mirrorH) == FALSE && isTRUE(mirrorH == 
        TRUE) == TRUE, crd[, 2] <- crd[, 2] * cos(pi) - crd[, 
        1] * sin(pi), mirrorY <- FALSE)
    if (isTRUE(loops == TRUE) == TRUE) {
        ifelse(missing(flip) == FALSE && isTRUE(flip == FALSE) == 
            TRUE, flip <- FALSE, flip <- TRUE)
        if (match.arg(layout) %in% c("circ", "conc")) {
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
            if (match.arg(layout) == "circ") {
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
                    else if (isTRUE(245 < angl[i]) == TRUE && 
                      isTRUE(295 > angl[i]) == TRUE) {
                      pos <- append(pos, 1)
                    }
                    else {
                      pos <- append(pos, 4)
                    }
                  }
                  rm(i)
                }
                else if (missing(pos) == FALSE) {
                  ifelse(isTRUE(length(pos) > 1) == TRUE && isTRUE("cn2" %in% 
                    attr(mlv, "class")) == TRUE, warning(paste("Only first value of ", 
                    sQuote("pos"), " is taken.", sep = ""), call. = FALSE), 
                    NA)
                }
                else {
                  NA
                }
            }
        }
        if (match.arg(layout) %in% c("stress", "force", "rand") || 
            missing(coord) == FALSE) {
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
            }
        }
    }
    if (missing(pos) == TRUE) {
        if (isTRUE("cn2" %in% attr(mlv, "class")) == TRUE || 
            (all(mlv$lbs$dm %in% mlv$lbs$cdm) == TRUE && all(mlv$lbs$cdm %in% 
                mlv$lbs$dm) == TRUE)) {
            pos <- c(0, 0)
        }
        else {
            pos <- c(2, 4)
        }
        ifelse(match.arg(layout) == "bip" && isTRUE(mirrorX == 
            TRUE) == TRUE, pos <- rev(pos), NA)
    }
    else if (missing(pos) == FALSE) {
        if (isTRUE("cn2" %in% attr(mlv, "class")) == TRUE || 
            (all(mlv$lbs$dm %in% mlv$lbs$cdm) == TRUE && all(mlv$lbs$cdm %in% 
                mlv$lbs$dm) == TRUE)) {
            ifelse(isTRUE(length(pos) > 1) == TRUE && match.arg(layout) != 
                "circ", warning(paste("Only first value of ", 
                sQuote("pos"), " is taken.", sep = ""), call. = FALSE), 
                NA)
            pos <- rep(pos[1], 2)
        }
        ifelse(isTRUE(length(pos) == 1) == TRUE, pos <- rep(pos, 
            2), pos <- pos[1:2])
        if (isTRUE(min(pos) < 0L) == TRUE | isTRUE(max(pos) > 
            4L) == TRUE) 
            stop(paste("Values in ", sQuote("pos"), " must be between 0 and 4.", 
                sep = ""), call. = FALSE)
    }
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
    graphics::par(bg = grDevices::adjustcolor(bg, alpha = alpha[3]))
    if (isTRUE(loops == TRUE) == TRUE && isTRUE(flip == TRUE) == 
        TRUE) {
        ylim <- c(min(nds[, 2]) - (max(cex)/100L), max(nds[, 
            2]) + (max(cex)/100L))
        xlim <- c(min(nds[, 1]) - (max(cex)/200L), max(nds[, 
            1]) + (max(cex)/200L))
    }
    else if (isTRUE(loops == TRUE) == TRUE) {
        ylim <- c(min(nds[, 2]) - ((cex[1])/200L), max(nds[, 
            2]) + ((cex[1])/200L))
        xlim <- c(min(nds[, 1]) - ((cex[1])/200L), max(nds[, 
            1]) + ((cex[1])/200L))
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
    ifelse(missing(new) == FALSE && isTRUE(new == TRUE) == TRUE, 
        suppressWarnings(par(new = TRUE)), NA)
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
    tlbsm <- vector()
    if (isTRUE(length(bdsm) > 0) == TRUE) {
        for (i in seq_len(length(attr(bdsm, "names")))) {
            ifelse(isTRUE(length(multiplex::dhc(attr(bdsm, "names")[i], 
                sep = "")) > 4L) == TRUE, tlbsm <- append(tlbsm, 
                tolower(paste(multiplex::dhc(attr(bdsm, "names")[i], 
                  sep = "")[1:4], collapse = ""))), tlbsm <- append(tlbsm, 
                tolower(attr(bdsm, "names"))[i]))
        }
        rm(i)
    }
    if (isTRUE(directed == FALSE) == TRUE && isTRUE(collRecip == 
        TRUE) == TRUE) {
        if (is.na(dim(rcpm)[3]) == TRUE) {
            trcpm <- multiplex::trnf(rcpm, tolist = TRUE)
        }
        else {
            trcpm <- multiplex::transf(rcpm, type = "tolist")
        }
    }
    netdrp <- net
    ifelse(isTRUE(valued == TRUE) == TRUE && isTRUE(max(net) > 
        10L) == TRUE, fnnetdrp <- (norm(as.matrix(netdrp), type = "F")), 
        NA)
    cx <- cex
    if (isTRUE(loops == TRUE) == TRUE && isTRUE(valued == TRUE) == 
        TRUE && isTRUE(max(netdrp) > 10L) == TRUE) {
        if (isTRUE(z == 1L) == TRUE) {
            diag(netdrp) <- (diag(netdrp)/fnnetdrp) * (15L)
        }
        else {
            for (i in seq_len(z)) {
                diag(netdrp[, , i]) <- as.vector((diag(netdrp[, 
                  , i])/fnnetdrp) * (15L))
            }
        }
    }
    else {
        NA
    }
    if (isTRUE(loops == TRUE) == TRUE) {
        tlbslp <- tlbs
        tlbs <- tlbs[which(tlbs != "loop")]
    }
    else {
        NA
    }
    if (isTRUE(collRecip == TRUE) == TRUE && isTRUE(c("recp") %in% 
        attr(bds, "names")) == TRUE) {
        trcp <- multiplex::trnf(rcp, tolist = TRUE)
    }
    else {
        NA
    }
    if (isTRUE(length(c(tlbs, tlbsm)) > 0) == TRUE) {
        for (k in seq_len(length(tlbsm))) {
            prs <- as.numeric(multiplex::dhc(bdsm[[k]]))
            pars <- as.matrix(nds[as.numeric(levels(factor(multiplex::dhc(bdsm[[k]])))), 
                ])
            rbdsm <- length(bdsm[[k]])
            if (isTRUE(rbdsm > 0L) == TRUE) {
                qn <- which(tlbsm[k] == attr(bdm, "names"))
                if (isTRUE(zz == 1L) == TRUE) {
                  ifelse(isTRUE(length(lty) == 1) == TRUE, vlt <- rep(Lt, 
                    rbdsm), vlt <- rep(Lt[r + 1], rbdsm))
                  ifelse(isTRUE(length(ecol) == 1L) == TRUE, 
                    vecolm <- rep(ecol[1], rbdsm), vecolm <- rep(ecol[z + 
                      1], rbdsm))
                  tbnd <- as.vector(unlist(bdm[qn]))
                  if (isTRUE(length(tbnd) > 0L) == TRUE) {
                    ifelse(isTRUE(any(tbnd %in% bdsm[[k]])) == 
                      TRUE, vlt <- append(vlt, rep(Lt, qn)), 
                      NA)
                    ifelse(isTRUE(any(tbnd %in% bdsm[[k]])) == 
                      TRUE, vltz <- append(vltz, rep(Lt, qn)), 
                      NA)
                  }
                  vltc <- rep(Lt[zz + 1], rbdsm)[1]
                }
                else if (isTRUE(zz > 1L) == TRUE) {
                  vlt <- vector()
                  for (i in seq_along(Lt)) {
                    tbnd <- as.vector(unlist(bdm[[qn]][i]))
                    if (isTRUE(length(tbnd) > 0L) == TRUE) {
                      ifelse(isTRUE(any(tbnd %in% bdsm[[k]])) == 
                        TRUE, vlt <- append(vlt, rep(Lt[(z + 
                        1):length(Lt)][i], length(which(tbnd %in% 
                        bdsm[[k]])))), NA)
                      ifelse(isTRUE(any(tbnd %in% bdsm[[k]])) == 
                        TRUE, vltz <- append(vltz, rep(Lt[(z + 
                        1):length(Lt)][i], length(which(tbnd %in% 
                        bdsm[[k]])))), NA)
                    }
                  }
                  rm(i)
                  if (isTRUE(length(lty) == 1L) == TRUE) {
                    vlt1 <- rep(lty, length(vlt))
                    vltc <- vlt
                  }
                  else {
                    vltc <- vector()
                    if (isTRUE(Lt %in% Ltc) == FALSE) {
                      for (i in seq_along(Ltc)) {
                        tbnd <- as.vector(unlist(bdm[[qn]][i]))
                        if (isTRUE(length(tbnd) > 0L) == TRUE) {
                          ifelse(isTRUE(any(tbnd %in% bdsm[[k]])) == 
                            TRUE, vltc <- append(vltc, rep(Ltc[(z + 
                            1):length(Ltc)][i], length(which(tbnd %in% 
                            bdsm[[k]])))), NA)
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
                            Lt[i])], vlt[which(vlt == Lt[i])] != 
                            i, i))
                        }
                        rm(i)
                      }
                    }
                  }
                }
                if (isTRUE(valued == TRUE) == TRUE && identical(met, 
                  multiplex::dichot(met, c = 1)) == FALSE) {
                  Lw <- vector()
                  i <- 1
                  for (j in seq_len(length(bdsm[[k]]))) {
                    qn <- c(prs[i], prs[(i + 1)])
                    if (isTRUE(zz == 1L) == TRUE) {
                      Lw <- append(Lw, met[qn[1], qn[2]])
                    }
                    else {
                      for (l in seq_len(zz)) {
                        Lw <- append(Lw, met[qn[1], qn[2], l] + 
                          t(met[qn[1], qn[2], l]))
                      }
                      rm(l)
                    }
                    i <- i + 2L
                  }
                  rm(j)
                  rm(i)
                  if (isTRUE(max(met) > 10L) == TRUE) {
                    lw <- (Lw/fnnetdrp) * (10L * 5L)
                  }
                  else {
                    lw <- Lw
                  }
                }
                else if (isTRUE(valued == TRUE) == TRUE) {
                  lw <- rep(lwd[1], rbdsm)
                }
                else {
                  ifelse(isTRUE(length(bdsm) == 0) == TRUE, NA, 
                    lw <- rep(lwd[1], rbdsm))
                }
                lwdfct <- (lw + (1L/lw)) * mscl
                if (isTRUE(collRecip == TRUE) == TRUE && isTRUE(tlbsm[k] %in% 
                  c("recp")) == TRUE) {
                  bw <- 0L
                  lw <- lwdfct
                }
                else if (isTRUE(collRecip == FALSE) == TRUE && 
                  isTRUE(tlbsm[k] %in% c("recp")) == TRUE) {
                  bw <- bwd
                }
                else {
                  bw <- bwd
                  lw <- lwdfct
                }
                hd <- 0L
                if (isTRUE(collRecip == TRUE) == TRUE && isTRUE(tlbsm[k] %in% 
                  c("recp")) == FALSE) {
                  flgcr <- numeric()
                  sbdsm <- multiplex::swp(bdsm[[k]])
                  if (isTRUE(directed == FALSE) == TRUE && any(sbdsm %in% 
                    unlist(trcpm)) == TRUE) {
                    for (i in seq_len(zz)) {
                      ifelse(any(sbdsm %in% trcpm[[i]]) == TRUE, 
                        flgcr <- append(flgcr, as.numeric(i)), 
                        NA)
                    }
                    rm(i)
                  }
                }
                else {
                  flgcr <- rep(0L, zz)
                }
                pars[, 1] <- pars[, 1] * scl[1]
                pars[, 2] <- pars[, 2] * scl[2]
                if (isTRUE(zz == 1L) == TRUE) {
                  ccbnd(pars, rbdsm, bdsm[[k]], vlt, cx * mscl, 
                    lw, vecolm, bw, alfa, fds, flgcx, flgcr, 
                    hd, m)
                }
                else {
                  ifelse(isTRUE(length(lty) == 1L) == TRUE, ccbnd(pars, 
                    rbdsm, bdsm[[k]], vlt1, cx * mscl, lw, vecol[vltc], 
                    bw, alfa, fds, flgcx, flgcr, hd, m), ccbnd(pars, 
                    rbdsm, bdsm[[k]], vlt, cx * mscl, lw, vecol[vltc], 
                    bw, alfa, fds, flgcx, flgcr, hd, m))
                }
            }
            else {
                NA
            }
        }
        rm(k)
        for (k in seq_len(length(tlbs))) {
            prs <- as.numeric(multiplex::dhc(bds[[k]]))
            pars <- as.matrix(nds[as.numeric(levels(factor(multiplex::dhc(bds[[k]])))), 
                ])
            rbds <- length(bds[[k]])
            if (isTRUE(rbds > 0L) == TRUE) {
                qn <- which(tlbs[k] == attr(bd, "names"))
                if (isTRUE(z == 1L) == TRUE) {
                  vlt <- rep(Lt, rbds)
                  vecol <- rep(ecol[1], rbds)
                  tbnd <- as.vector(unlist(bd[qn]))
                  if (isTRUE(length(tbnd) > 0L) == TRUE) {
                    ifelse(isTRUE(any(tbnd %in% bds[[k]])) == 
                      TRUE, vlt <- append(vlt, rep(Lt, qn)), 
                      NA)
                    ifelse(isTRUE(any(tbnd %in% bds[[k]])) == 
                      TRUE, vltz <- append(vltz, rep(Lt, qn)), 
                      NA)
                  }
                  vltc <- vlt[1]
                }
                else if (isTRUE(z > 1L) == TRUE) {
                  vlt <- vector()
                  for (i in seq_along(Lt)) {
                    tbnd <- as.vector(unlist(bd[[qn]][i]))
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
                    if (isTRUE(Lt %in% Ltc) == FALSE) {
                      for (i in seq_along(Ltc)) {
                        tbnd <- as.vector(unlist(bd[[qn]][i]))
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
                if (isTRUE(flgcx == TRUE) == FALSE) {
                  cx <- rep(cex[1], 2)
                }
                if (isTRUE(valued == TRUE) == TRUE) {
                  Lw <- vector()
                  i <- 1
                  for (j in seq_len(length(bds[[k]]))) {
                    qn <- c(prs[i], prs[(i + 1)])
                    if (isTRUE(collRecip == TRUE) == TRUE) {
                      ifelse(isTRUE(z == 1L) == TRUE, Lw <- append(Lw, 
                        netdrp[qn[1], qn[2]]), Lw <- append(Lw, 
                        netdrp[qn[1], qn[2], vltc[j]] + t(netdrp[qn[1], 
                          qn[2], vltc[j]])))
                    }
                    else if (isTRUE(collRecip == FALSE) == TRUE) {
                      ifelse(isTRUE(z == 1L) == TRUE, Lw <- append(Lw, 
                        netdrp[qn[1], qn[2]]), Lw <- append(Lw, 
                        netdrp[qn[1], qn[2], vltc[j]]))
                    }
                    i <- i + 2L
                  }
                  if (isTRUE(max(netdrp) > 10L) == TRUE) {
                    lw <- (Lw/fnnetdrp) * (10L * 3L)
                  }
                  else {
                    lw <- Lw
                  }
                }
                else if (isTRUE(valued == TRUE) == FALSE) {
                  ifelse(isTRUE(length(bds) == 0) == TRUE, NA, 
                    lw <- rep(lwd[1], rbds))
                }
                if (isTRUE(collRecip == TRUE) == TRUE && isTRUE(tlbs[k] %in% 
                  c("recp")) == TRUE) {
                  bw <- 0L
                  lw <- lw * mscl
                }
                else {
                  bw <- bwd
                  lw <- lw * mscl
                }
                ifelse(isTRUE(directed == FALSE) == TRUE, hd <- 0L, 
                  hd <- hds * 0.75)
                if (isTRUE(collRecip == TRUE) == TRUE && isTRUE(tlbs[k] %in% 
                  c("recp")) == FALSE && isTRUE(valued == TRUE) == 
                  FALSE && isTRUE(c("recp") %in% attr(bds, "names")) == 
                  TRUE) {
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
                if (match.arg(layout) == "bip") {
                  if (missing(elv) == TRUE) {
                    elv <- 0.25
                  }
                  else {
                    ifelse(isTRUE(elv > 1L) == TRUE, elv <- 1L, 
                      NA)
                  }
                  bzrc((pars), cex = cx, lty = vlt, col = vecol[vltc], 
                    lwd = lw, elv = elv, ...)
                }
                else {
                  if (isTRUE(z == 1L) == TRUE) {
                    ccbnd(pars, rbds, bds[[k]], vlt, cx * mscl, 
                      lw, vecol, bw, alfa, fds, flgcx, flgcr, 
                      hd, n)
                  }
                  else {
                    ifelse(isTRUE(length(lty) == 1L) == TRUE, 
                      ccbnd(pars, rbds, bds[[k]], vlt1, cx * 
                        mscl, lw, vecol[vltc], bw, alfa, fds, 
                        flgcx, flgcr, hd, n), ccbnd(pars, rbds, 
                        bds[[k]], vlt, cx * mscl, lw, vecol[vltc], 
                        bw, alfa, fds, flgcx, flgcr, hd, n))
                  }
                }
            }
            else {
                NA
            }
        }
        rm(k)
    }
    else {
        ifelse(isTRUE(loops == TRUE) == TRUE && isTRUE(length(lty) == 
            1) == TRUE, vlt <- lty, NA)
    }
    if (isTRUE(loops == TRUE) == TRUE) {
        bdlp <- bd$loop
        if (missing(vldist) == FALSE) {
            if (isTRUE(abs(vldist) > cx[1]) == TRUE) {
                ifelse(isTRUE(vldist < 0) == TRUE, vldist <- cx[1] * 
                  -1L, vldist <- cx[1] * 1L)
            }
        }
        else {
            NA
        }
        dz <- (rng(z) + abs(min(rng(z))))/(10L)
        ndss <- nds
        ndss[, 1] <- ndss[, 1] * scl[1]
        ndss[, 2] <- ndss[, 2] * scl[2]
        if (isTRUE(length(lty) == 1) == TRUE) {
            Ltl <- rep(lty, length(bdlp))
        }
        else {
            Ltl <- Lt
        }
        if (isTRUE(z == 1L) == TRUE) {
            lp <- as.vector(which(diag(net) > 0))
            if (isTRUE(length(lp) > 0) == TRUE) {
                if ((missing(coord) == FALSE || match.arg(layout) %in% 
                  c("circ", "conc", "force", "stress", "rand")) && 
                  isTRUE(flip == TRUE) == TRUE) {
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
                for (i in seq_len(length(lp))) {
                  if (isTRUE(n == 1L) == TRUE) {
                    dcx <- (cex[lp[i]]^2) * 0.00025
                    rds <- (cex[lp[i]]^2) * 0.00015
                  }
                  else if (isTRUE(n < 3) == TRUE) {
                    dcx <- (cex[lp[i]] * 0.0075)
                    rds <- (cex[lp[i]] * 0.005) - (dz)
                  }
                  else {
                    dcx <- (cex[lp[i]] * 0.01)
                    rds <- (cex[lp[i]] * 0.0075) - (dz)
                  }
                  ifelse(isTRUE(mirrorY == TRUE) == TRUE && isTRUE(flip == 
                    TRUE) == TRUE, flipl <- !flipl, NA)
                  if (isTRUE(flipl[lp[i]] == TRUE) == TRUE) {
                    dcx <- dcx * (-1e-17)
                    ifelse(missing(vldist) == FALSE, dcx <- dcx - 
                      (vldist * 0.008), NA)
                  }
                  else {
                    ifelse(missing(vldist) == FALSE, dcx <- dcx + 
                      (vldist * 0.01), NA)
                  }
                  if (isTRUE(valued == TRUE) == TRUE) {
                    lwl <- diag(netdrp)[lp[i]]
                  }
                  else {
                    if (isTRUE(Ltl == 1) == TRUE) {
                      lwl <- lwd - 0.5
                    }
                    else {
                      lwl <- lwd
                    }
                  }
                  hc(x = ndss[lp[i], 1], y = ndss[lp[i], 2] + 
                    (dcx), r = rds, col = grDevices::adjustcolor(vecol, 
                    alpha = alfa), lty = Ltl, lwd = lwl, flip = flipl[lp[i]])
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
                  bwd2 <- 4L
                }
            }
            if (isTRUE(bwd2 == 0) == TRUE | isTRUE(max(cex) < 
                3) == TRUE) {
                dz <- rep(0, z)
            }
            else {
                ifelse(missing(bwd2) == TRUE, dz <- (rng(z) + 
                  abs(min(rng(z))))/(10L), dz <- (bwd2) * (rng(z) + 
                  abs(min(rng(z))))/(20L))
            }
            ifelse(isTRUE(length(lwd) == 1) == TRUE, lwd <- rep(lwd, 
                z), NA)
            for (k in seq_len(length(bdlp))) {
                lp <- as.numeric(unique(multiplex::dhc(bdlp)[k][[1]]))
                if (isTRUE(length(lp) > 0) == TRUE) {
                  for (i in seq_len(length(lp))) {
                    if ((missing(coord) == FALSE || match.arg(layout) %in% 
                      c("circ", "conc", "force", "stress", "rand")) && 
                      isTRUE(flip == TRUE) == TRUE) {
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
                    if (isTRUE(n == 1L) == TRUE) {
                      dcx <- (cex[lp[i]]^2) * 0.00025
                      rds <- abs((cex[lp[i]]^2) * 0.00015)
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
                        rds <- rds * 0.85
                        dcx <- dcx - (cex[lp[i]] * 0.011)
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
                    if (isTRUE(valued == TRUE) == TRUE) {
                      lwl <- netdrp[i, i, k]
                    }
                    else {
                      if (isTRUE(Ltl[k] == 1) == TRUE) {
                        lwl <- lwd[k] - k
                      }
                      else {
                        lwl <- lwd[k]
                      }
                    }
                    hc(x = ndss[lp[i], 1], y = ndss[lp[i], 2] + 
                      (dcx), r = rds, col = grDevices::adjustcolor(vecol[k], 
                      alpha = alfa), lty = Ltl[k], lwd = lwl, 
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
        nn <- length(mlv$lbs$dm)
        if (isTRUE(pos[1] == 0) == TRUE) {
            if (missing(fstyle) == TRUE || (missing(fstyle) == 
                FALSE && isTRUE(fstyle %in% c("italic", "bold", 
                "bolditalic") == FALSE))) {
                graphics::text(ndss[seq_len(nn), ], labels = lbs[seq_len(nn)], 
                  cex = fsize, adj = 0.5, col = fcol[1])
            }
            else if (missing(fstyle) == FALSE) {
                if (isTRUE(fstyle == "italic") == TRUE) {
                  graphics::text(ndss[seq_len(nn), ], labels = as.expression(lapply(lbs[seq_len(nn)], 
                    function(w) bquote(italic(.(w))))), cex = fsize, 
                    adj = 0.5, col = fcol[1])
                }
                else if (isTRUE(fstyle == "bold") == TRUE) {
                  graphics::text(ndss[seq_len(nn), ], labels = as.expression(lapply(lbs[seq_len(nn)], 
                    function(w) bquote(bold(.(w))))), cex = fsize, 
                    adj = 0.5, col = fcol[1])
                }
                else if (isTRUE(fstyle == "bolditalic") == TRUE) {
                  graphics::text(ndss[seq_len(nn), ], labels = as.expression(lapply(lbs[seq_len(nn)], 
                    function(w) bquote(bolditalic(.(w))))), cex = fsize, 
                    adj = 0.5, col = fcol[1])
                }
            }
        }
        else {
            if (missing(fstyle) == TRUE || (missing(fstyle) == 
                FALSE && isTRUE(fstyle %in% c("italic", "bold", 
                "bolditalic") == FALSE))) {
                graphics::text(ndss[seq_len(nn), ], lbs[seq_len(nn)], 
                  cex = fsize, pos = pos[1], col = fcol[1], offset = (cex/4L), 
                  adj = c(0.5, 1))
            }
            else if (missing(fstyle) == FALSE) {
                if (isTRUE(fstyle == "italic") == TRUE) {
                  graphics::text(ndss[seq_len(nn), ], as.expression(lapply(lbs[seq_len(nn)], 
                    function(w) bquote(italic(.(w))))), cex = fsize, 
                    pos = pos[1], col = fcol[1], offset = (cex/4L), 
                    adj = c(0.5, 1))
                }
                else if (isTRUE(fstyle == "bold") == TRUE) {
                  graphics::text(ndss[seq_len(nn), ], as.expression(lapply(lbs[seq_len(nn)], 
                    function(w) bquote(bold(.(w))))), cex = fsize, 
                    pos = pos[1], col = fcol[1], offset = (cex/4L), 
                    adj = c(0.5, 1))
                }
                else if (isTRUE(fstyle == "bolditalic") == TRUE) {
                  graphics::text(ndss[seq_len(nn), ], as.expression(lapply(lbs[seq_len(nn)], 
                    function(w) bquote(bolditalic(.(w))))), cex = fsize, 
                    pos = pos[1], col = fcol[1], offset = (cex/4L), 
                    adj = c(0.5, 1))
                }
            }
        }
        if (isTRUE(pos[2] == 0) == TRUE) {
            if (missing(fstyle) == TRUE || (missing(fstyle) == 
                FALSE && isTRUE(fstyle %in% c("italic", "bold", 
                "bolditalic") == FALSE))) {
                graphics::text(ndss[(nn + 1L):dim(mlv$mlnet)[1], 
                  ], labels = lbs[(nn + 1L):dim(mlv$mlnet)[1]], 
                  cex = fsize, adj = 0.5, col = fcol[2])
            }
            else if (missing(fstyle) == FALSE) {
                if (isTRUE(fstyle == "italic") == TRUE) {
                  graphics::text(ndss[(nn + 1L):dim(mlv$mlnet)[1], 
                    ], labels = as.expression(lapply(lbs[(nn + 
                    1L):dim(mlv$mlnet)[1]], function(w) bquote(italic(.(w))))), 
                    cex = fsize, adj = 0.5, col = fcol[2])
                }
                else if (isTRUE(fstyle == "bold") == TRUE) {
                  graphics::text(ndss[(nn + 1L):dim(mlv$mlnet)[1], 
                    ], labels = as.expression(lapply(lbs[(nn + 
                    1L):dim(mlv$mlnet)[1]], function(w) bquote(bold(.(w))))), 
                    cex = fsize, adj = 0.5, col = fcol[2])
                }
                else if (isTRUE(fstyle == "bolditalic") == TRUE) {
                  graphics::text(ndss[(nn + 1L):dim(mlv$mlnet)[1], 
                    ], labels = as.expression(lapply(lbs[(nn + 
                    1L):dim(mlv$mlnet)[1]], function(w) bquote(bolditalic(.(w))))), 
                    cex = fsize, adj = 0.5, col = fcol[2])
                }
            }
        }
        else {
            if (missing(fstyle) == TRUE || (missing(fstyle) == 
                FALSE && isTRUE(fstyle %in% c("italic", "bold", 
                "bolditalic") == FALSE))) {
                graphics::text(ndss[(nn + 1L):dim(mlv$mlnet)[1], 
                  ], lbs[(nn + 1L):dim(mlv$mlnet)[1]], cex = fsize, 
                  pos = pos[2], col = fcol[2], offset = (cex/4L), 
                  adj = c(0.5, 1))
            }
            else if (missing(fstyle) == FALSE) {
                if (isTRUE(fstyle == "italic") == TRUE) {
                  graphics::text(ndss[(nn + 1L):dim(mlv$mlnet)[1], 
                    ], as.expression(lapply(lbs[(nn + 1L):dim(mlv$mlnet)[1]], 
                    function(w) bquote(italic(.(w))))), cex = fsize, 
                    pos = pos[2], col = fcol[2], offset = (cex/4L), 
                    adj = c(0.5, 1))
                }
                else if (isTRUE(fstyle == "bold") == TRUE) {
                  graphics::text(ndss[(nn + 1L):dim(mlv$mlnet)[1], 
                    ], as.expression(lapply(lbs[(nn + 1L):dim(mlv$mlnet)[1]], 
                    function(w) bquote(bold(.(w))))), cex = fsize, 
                    pos = pos[2], col = fcol[2], offset = (cex/4L), 
                    adj = c(0.5, 1))
                }
                else if (isTRUE(fstyle == "bolditalic") == TRUE) {
                  graphics::text(ndss[(nn + 1L):dim(mlv$mlnet)[1], 
                    ], as.expression(lapply(lbs[(nn + 1L):dim(mlv$mlnet)[1]], 
                    function(w) bquote(bolditalic(.(w))))), cex = fsize, 
                    pos = pos[2], col = fcol[2], offset = (cex/4L), 
                    adj = c(0.5, 1))
                }
            }
        }
    }
    if (isTRUE(showAtts == TRUE) == TRUE) {
        ndss <- nds
        ndss[, 1] <- ndss[, 1] * scl[1]
        ndss[, 2] <- ndss[, 2] * scl[2]
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
                        diag = TRUE)) != 0)] <- lbat
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
            graphics::text(ndss, labels = atts, cex = fsize, 
                pos = pos%%4 + 1L, col = fcol, offset = (cex/4L), 
                adj = c(0.5, 1))
        }
        else if (isTRUE(flgcx == TRUE) == TRUE) {
            graphics::text(ndss, labels = atts, cex = fsize, 
                pos = pos%%4 + 1L, col = fcol, offset = (min(cex)/4L), 
                adj = c(0.5, 1))
        }
    }
    graphics::par(mar = omr)
    graphics::par(bg = obg)
    graphics::par(lend = 0)
    graphics::par(mai = omi)
}

