multigraph <-
function (net, layout = c("circ", "force", "stress", "conc", 
    "rand"), scope, directed = TRUE, main = NULL, opt, seed = NULL, 
    maxiter = 100, alpha = c(1, 1, 1), collRecip, showLbs, showAtts, 
    cex.main, weighted, weights, cex2, coord, clu, cex, lwd, 
    pch, lty, bwd, tcol, tcex, att, bg, mar, pos, asp, ecol, 
    vcol, vcol0, tcex2, tcol2, hds, vedist, rot, mirrorX, mirrorY, 
    col, lbat, drp, swp, loops, swp2, signed, scl, add, mirrorD, 
    mirrorL, lbs, mirrorV, mirrorH, ...) 
{
    if (is.null(dim(net)) == TRUE || (isTRUE(dim(net)[1] != dim(net)[2]) == 
        TRUE)) {
        stop("'net' must be matrix, array, or a list of a \"Signed\" or \"Multilevel\" class.")
    }
    else {
        NA
    }
    if (isTRUE(attr(net, "class") == "Multilevel") == TRUE) {
        mlvl <- net
        net <- mlvl$bmat
        if (missing(pch) == TRUE) {
            pch <- c(rep(1, length(mlvl$lbs$dm)), rep(0, length(mlvl$lbs$cdm)))
        }
        else {
            NA
        }
    }
    else if (isTRUE(attr(net, "class") == "Data Set") == TRUE) {
        att <- net$net[, , which(net$at == 1)]
        net <- net$net[, , which(net$at == 0)]
        ifelse(missing(pch) == TRUE, pch <- 21, NA)
    }
    else {
        ifelse(missing(pch) == TRUE, pch <- 21, NA)
    }
    if (isTRUE(is.list(net) == TRUE) == TRUE && isTRUE(attr(net, 
        "class") == "Signed") == FALSE) {
        net <- multiplex::transf(net, type = "toarray", lb2lb = TRUE, 
            lbs = sort(unique(multiplex::dhc(unlist(net)))))
    }
    else if (isTRUE(is.vector(net) == TRUE) == TRUE) {
        ifelse(missing(add) == FALSE && isTRUE(is.vector(add) == 
            TRUE) == TRUE, net <- multiplex::transf(net, type = "toarray", 
            lb2lb = TRUE, lbs = sort(unique(multiplex::dhc(c(net, 
                add))))), net <- multiplex::trnf(net, tolist = FALSE, 
            lb2lb = TRUE))
    }
    else if (isTRUE(is.null(net) == TRUE) == TRUE && (missing(add) == 
        FALSE && isTRUE(is.vector(add) == TRUE) == TRUE)) {
        net <- matrix(0L, nrow = length(add), ncol = length(add), 
            dimnames = list(add, add))
    }
    else {
        NA
    }
    ifelse(missing(signed) == FALSE && isTRUE(signed == TRUE) == 
        TRUE, signed <- TRUE, signed <- FALSE)
    ifelse(isTRUE(dim(net)[3] == 1) == TRUE, net <- net[, , 1], 
        NA)
    ifelse(missing(weighted) == FALSE && isTRUE(weighted == TRUE) == 
        TRUE, weighted <- TRUE, weighted <- FALSE)
    ifelse(missing(weights) == FALSE && isTRUE(weights == TRUE) == 
        TRUE, weights <- TRUE, weights <- FALSE)
    ifelse(missing(loops) == FALSE && isTRUE(loops == TRUE) == 
        TRUE, loops <- TRUE, loops <- FALSE)
    ifelse(missing(mirrorH) == FALSE && isTRUE(mirrorH == TRUE) == 
        TRUE, mirrorY <- TRUE, NA)
    ifelse(missing(mirrorV) == FALSE && isTRUE(mirrorV == TRUE) == 
        TRUE, mirrorX <- TRUE, NA)
    if (missing(scope) == FALSE && isTRUE(any(attr(scope, "names") == 
        "directed") == TRUE) == TRUE) {
        if (isTRUE(scope[[which(attr(scope, "names") == "directed")]] == 
            FALSE) == TRUE) {
            ifelse(missing(collRecip) == FALSE && isTRUE(collRecip == 
                TRUE) == FALSE, collRecip <- FALSE, collRecip <- TRUE)
        }
        else if (isTRUE(scope[[which(attr(scope, "names") == 
            "directed")]] == TRUE) == TRUE) {
            ifelse(missing(collRecip) == FALSE && isTRUE(collRecip == 
                TRUE) == TRUE, collRecip <- TRUE, collRecip <- FALSE)
        }
    }
    else {
        if (isTRUE(directed == FALSE) == TRUE) {
            ifelse(missing(collRecip) == FALSE && isTRUE(collRecip == 
                TRUE) == FALSE, collRecip <- FALSE, collRecip <- TRUE)
        }
        else if (isTRUE(directed == TRUE) == TRUE) {
            ifelse(missing(collRecip) == FALSE && isTRUE(collRecip == 
                TRUE) == TRUE, collRecip <- TRUE, collRecip <- FALSE)
        }
    }
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
    ifelse(missing(swp) == FALSE && isTRUE(swp == TRUE) == TRUE, 
        swp <- TRUE, swp <- FALSE)
    ifelse(missing(swp2) == FALSE && isTRUE(swp2 == TRUE) == 
        TRUE, swp2 <- TRUE, swp2 <- FALSE)
    ifelse(isTRUE(directed == FALSE) == TRUE, directed <- FALSE, 
        NA)
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
    ifelse(missing(tcol) == TRUE, tcol <- 1, NA)
    ifelse(missing(bwd) == TRUE, bwd <- 1, NA)
    ifelse(isTRUE(bwd < 0L) == TRUE, bwd <- 0L, NA)
    if (missing(bg) == TRUE) {
        bg <- graphics::par()$bg
        flgbg <- 0L
    }
    else {
        obg <- graphics::par()$bg
        flgbg <- 1L
    }
    ifelse(isTRUE(weights == TRUE) == TRUE && isTRUE(max(net) > 
        99L) == TRUE, flgbg <- 3L, NA)
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
    ifelse(isTRUE(dim(net)[1] > 8) == TRUE, hds <- hds * 0.75, 
        NA)
    ifelse(missing(scl) == TRUE, scl <- rep(1, 2), NA)
    ifelse(isTRUE(length(scl) == 1) == TRUE, scl <- rep(scl, 
        2), scl <- scl[1:2])
    ifelse(missing(vedist) == TRUE, vedist <- 0, NA)
    ifelse(isTRUE(vedist > 1L) == TRUE, vedist <- 1L, NA)
    if (isTRUE(signed == TRUE) == TRUE || isTRUE(attr(net, "class") == 
        "Signed") == TRUE) {
        if (isTRUE(attr(net, "class") == "Signed") == TRUE) {
            if (any(net$val %in% c(-1, 0, 1)) == TRUE) {
                net <- multiplex::zbind(multiplex::dichot(net$s, 
                  c = 1L), 1 - multiplex::dichot(net$s, c = 0))
            }
            else {
                nets <- multiplex::zbnd(net$s, net$s)
                net <- array(NA, dim = c(dim(nets)[1], dim(nets)[2], 
                  2), dimnames = list(dimnames(net$s)[[1]], dimnames(net$s)[[1]]))
                net <- replace(net, nets == "o", 0)
                net <- replace(net, nets == "a", 1)
                net[, , 1] <- replace(net[, , 1], nets[, , 1] == 
                  "p", 1)
                net[, , 1] <- replace(net[, , 1], nets[, , 1] == 
                  "n", 0)
                net[, , 2] <- replace(net[, , 2], nets[, , 2] == 
                  "n", 1)
                net[, , 2] <- replace(net[, , 2], nets[, , 2] == 
                  "p", 0)
                rm(nets)
            }
        }
        else {
            ifelse(isTRUE(is.na(dim(net)[3]) == TRUE) == TRUE, 
                NA, net <- net[, , seq_len(2)])
        }
        ifelse(is.null(dimnames(net)[[1]]) == TRUE | isTRUE(showLbs == 
            FALSE) == TRUE, showLbs <- FALSE, showLbs <- TRUE)
    }
    if (isTRUE(attr(net, "class") == "Signed") == FALSE) {
        if (missing(lbs) == FALSE) {
            NA
        }
        else {
            ifelse(is.null(dimnames(net)[[1]]) == TRUE, lbs <- as.character(seq_len(dim(net)[1])), 
                lbs <- dimnames(net)[[1]])
        }
    }
    n <- dim(net)[1]
    ifelse(isTRUE(is.na(dim(net)[3]) == TRUE) == TRUE, z <- 1L, 
        z <- dim(net)[3])
    ifelse(isTRUE(swp == TRUE) == TRUE && isTRUE(z > 1L) == TRUE, 
        net <- net[, , rev(seq_len(z))], NA)
    if (missing(att) == FALSE && is.array(att) == TRUE) {
        if (isTRUE(n != dim(att)[1]) == TRUE) {
            warning("Dimensions in \"net\" and \"att\" differ. No attributes are shown.")
            showAtts <- FALSE
        }
    }
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
            for (i in seq_len(z)) {
                netdrp[, , i] <- netdrp[, , i] + t(netdrp[, , 
                  i])
            }
            rm(i)
        }
    }
    else {
        NA
    }
    if (isTRUE(collRecip == TRUE) == TRUE && (isTRUE(weighted == 
        TRUE) == FALSE | isTRUE(directed == FALSE) == TRUE)) {
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
    if (missing(ecol) == TRUE) {
        ifelse(isTRUE(signed == TRUE) == TRUE, ecol <- rep(1, 
            2), ecol <- grDevices::gray.colors(r, start = 0.1, 
            end = 0.5))
    }
    else {
        ifelse(isTRUE(ecol == 0) == TRUE, ecol <- "#FFFFFF", 
            NA)
    }
    if (isTRUE(weighted == TRUE) == TRUE) {
        ifelse(missing(lty) == TRUE, lty <- rep(1, r), lty <- rep(lty, 
            r))
    }
    else if (isTRUE(signed == TRUE) == TRUE) {
        ifelse(missing(lty) == TRUE, lty <- c(1, 3), NA)
    }
    else {
        ifelse(missing(lty) == TRUE, lty <- seq_len(r), NA)
    }
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
                ifelse(isTRUE(swp == TRUE) == TRUE && isTRUE(weighted == 
                  TRUE) == FALSE, Ltc <- rev(seq_len(r)), Ltc <- seq_len(r))
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
    else if (isTRUE(flgcx == FALSE) == TRUE && isTRUE(weighted == 
        TRUE) == TRUE) {
        ifelse(isTRUE(max(cex) >= 21L) == TRUE, cex <- 20L, NA)
    }
    else {
        NA
    }
    if (missing(tcex) == TRUE) {
        ifelse(isTRUE(max(cex) < 2) == TRUE, tcex <- cex * 0.66, 
            tcex <- cex * 0.25)
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
    ifelse(isTRUE(directed == TRUE) == TRUE, fds <- 125L - (n * 
        ffds), fds <- 145L - (n * ffds))
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
        switch(match.arg(layout), force = {
            crd <- frcd(netd, seed = seed, maxiter = maxiter)
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
        ifelse(isTRUE(n == 2) == TRUE && isTRUE(rot == -90) == 
            TRUE, rot <- -89.9, NA)
    }
    else {
        cnt <- 0
    }
    if (isTRUE(flgcrd == FALSE) == TRUE) {
        if (match.arg(layout) == "circ" && missing(pos) == TRUE) {
            angl <- vector()
            length(angl) <- n
            for (i in seq_len(n)) {
                ifelse((atan2((crd[i, 2] - cnt), (crd[i, 1] - 
                  cnt)) * (180L/pi)) >= 0, angl[i] <- atan2((crd[i, 
                  2] - cnt), (crd[i, 1] - cnt)) * (180L/pi), 
                  angl[i] <- ((atan2((crd[i, 2] - cnt), (crd[i, 
                    1] - cnt)) * (180L/pi))%%180L) + 180L)
            }
            rm(i)
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
    if (missing(pos) == TRUE) {
        pos <- 4
    }
    else {
        if (isTRUE(pos < 0L) == TRUE | isTRUE(pos > 4L) == TRUE) 
            stop("Invalid \"pos\" value.")
    }
    ifelse(missing(mirrorX) == FALSE && isTRUE(mirrorX == TRUE) == 
        TRUE || missing(mirrorV) == FALSE && isTRUE(mirrorV == 
        TRUE) == TRUE, crd[, 1] <- crd[, 1] * cos(pi) - crd[, 
        2] * sin(pi), mirrorX <- FALSE)
    ifelse(missing(mirrorY) == FALSE && isTRUE(mirrorY == TRUE) == 
        TRUE || missing(mirrorH) == FALSE && isTRUE(mirrorH == 
        TRUE) == TRUE, crd[, 2] <- crd[, 2] * cos(pi) - crd[, 
        1] * sin(pi), mirrorY <- FALSE)
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
    tcex <- tcex * mscl
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
    if (missing(opt) == FALSE && isTRUE(opt == TRUE) == TRUE) {
        xlim <- xlim + c(-0.1, 0.1)
        ylim <- ylim + c(0, 0.1)
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
    ifelse(isTRUE(weighted == TRUE) == TRUE && isTRUE(max(netdrp) > 
        10L) == TRUE, fnnetdrp <- (norm(as.matrix(netdrp), type = "F")), 
        NA)
    if (isTRUE(swp == TRUE) == TRUE) {
        Lt <- Lt[rev(seq_len(length(Lt)))]
        lwd <- lwd[length(lwd):1]
        ifelse(isTRUE(weighted == TRUE) == TRUE, vecol <- vecol[rev(seq_len(length(vecol)))], 
            NA)
        alfa <- alfa[rev(seq_len(length(alfa)))]
    }
    if (isTRUE(loops == TRUE) == TRUE && isTRUE(weighted == TRUE) == 
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
    if (isTRUE(collRecip == TRUE) == TRUE && isTRUE(weighted == 
        TRUE) == FALSE) {
        trcp <- multiplex::transf(rcp, type = "tolist")
        utrcp <- unlist(trcp)
    }
    else {
        flgcr <- rep(0L, z)
    }
    if (isTRUE(weights == TRUE) == TRUE) {
        mdp <- data.frame(matrix(nrow = 0L, ncol = 3L))
    }
    else {
        NA
    }
    if (isTRUE(length(tlbs) > 0) == TRUE) {
        for (k in seq_len(length(tlbs))) {
            prs <- as.numeric(multiplex::dhc(bds[[k]]))
            if (isTRUE(directed == TRUE) == TRUE | isTRUE(collRecip == 
                FALSE) == TRUE | (isTRUE(directed == FALSE) == 
                TRUE && isTRUE(tlbs[k] %in% c("txch", "mixd", 
                "full")) == TRUE)) {
                pars <- as.matrix(nds[as.numeric(levels(factor(multiplex::dhc(bds[[k]])))), 
                  ])
            }
            else {
                pars <- as.matrix(nds[prs, ])
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
                      c("asym", "tent", "txch", "mixd")) == TRUE) | 
                    (isTRUE(directed == TRUE) == FALSE && isTRUE(tlbs[k] %in% 
                      c("txch")) == TRUE)) {
                    if (isTRUE(multiplex::dhc(bds[[k]][1])[1] > 
                      multiplex::dhc(bds[[k]][1])[2]) == TRUE) {
                      Prs <- vector()
                      length(Prs) <- length(prs)
                      Prs[which(seq_len(length(prs))%%2L == 1L)] <- prs[which(seq_len(length(prs))%%2L == 
                        0L)]
                      Prs[which(seq_len(length(prs))%%2L == 0L)] <- prs[which(seq_len(length(prs))%%2L == 
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
                  cx <- cex
                }
                if (isTRUE(weighted == TRUE) == TRUE) {
                  Lw <- vector()
                  i <- 1
                  for (j in seq_len(length(bds[[k]]))) {
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
                ifelse(isTRUE(swp2 == TRUE) == TRUE && isTRUE(tlbs[k] %in% 
                  c("recp")) == TRUE, bds[[k]] <- multiplex::swp(bds[[k]]), 
                  NA)
                ifelse(isTRUE(collRecip == TRUE) == TRUE && isTRUE(directed == 
                  TRUE) == TRUE && isTRUE(tlbs[k] %in% c("recp")) == 
                  TRUE, bw <- 0, bw <- bwd)
                if ((isTRUE(collRecip == TRUE) == TRUE && isTRUE(directed == 
                  TRUE) == TRUE) && isTRUE(tlbs[k] %in% c("recp")) == 
                  FALSE && isTRUE(weighted == TRUE) == FALSE) {
                  flgcr <- numeric()
                  sbds <- multiplex::swp(bds[[k]])
                  if (any(sbds %in% utrcp) == TRUE) {
                    for (i in seq_len(z)) {
                      ifelse(any(sbds %in% trcp[[i]]) == TRUE, 
                        flgcr <- append(flgcr, as.numeric(i)), 
                        NA)
                    }
                    rm(i)
                  }
                }
                else {
                  NA
                }
                pars[, 1] <- pars[, 1] * scl[1]
                pars[, 2] <- pars[, 2] * scl[2]
                if (isTRUE(z == 1L) == TRUE) {
                  mbnd(pars, rr, bds[[k]], vlt, cx * mscl, lwd * 
                    mscl, vecol, directed, bw, alfa, fds, flgcx, 
                    weighted, flgcr, hds, n)
                }
                else {
                  ifelse(isTRUE(length(lty) == 1L) == TRUE, mbnd(pars, 
                    rr, bds[[k]], vlt1, cx * mscl, lwd * mscl, 
                    vecol[vltc], directed, bw, alfa, fds, flgcx, 
                    weighted, flgcr, hds, n), mbnd(pars, rr, 
                    bds[[k]], vlt, cx * mscl, lwd * mscl, vecol[vltc], 
                    directed, bw, alfa, fds, flgcx, weighted, 
                    flgcr, hds, n))
                }
                if (isTRUE(weights == TRUE) == TRUE) {
                  if (isTRUE(z == 1L) == TRUE) {
                    ifelse(isTRUE(directed == FALSE) == TRUE, 
                      mdp[nrow(mdp) + 1, ] <- c((pars[1, 1] + 
                        pars[2, 1])/2L, (pars[1, 2] + pars[2, 
                        2])/2L, netdrp[as.numeric(multiplex::dhc(bds[[k]])[1]), 
                        as.numeric(multiplex::dhc(bds[[k]])[2])]/2L), 
                      mdp[nrow(mdp) + 1, ] <- c((pars[1, 1] + 
                        pars[2, 1])/2L, (pars[1, 2] + pars[2, 
                        2])/2L, netdrp[as.numeric(multiplex::dhc(bds[[k]])[1]), 
                        as.numeric(multiplex::dhc(bds[[k]])[2])]))
                  }
                  else {
                    ifelse(isTRUE(directed == FALSE) == TRUE, 
                      mdp[nrow(mdp) + 1, ] <- c((pars[1, 1] + 
                        pars[2, 1])/2L, (pars[1, 2] + pars[2, 
                        2])/2L, multiplex::mnplx(netdrp, dichot = FALSE, 
                        directed = TRUE)[as.numeric(multiplex::dhc(bds[[k]])[1]), 
                        as.numeric(multiplex::dhc(bds[[k]])[2])]), 
                      mdp[nrow(mdp) + 1, ] <- c((pars[1, 1] + 
                        pars[2, 1])/2L, (pars[1, 2] + pars[2, 
                        2])/2L, multiplex::mnplx(netdrp, dichot = FALSE, 
                        directed = FALSE)[as.numeric(multiplex::dhc(bds[[k]])[1]), 
                        as.numeric(multiplex::dhc(bds[[k]])[2])]))
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
        NA
    }
    if (isTRUE(loops == TRUE) == TRUE) {
        if (isTRUE(swp == TRUE) == TRUE) {
            bdlp <- bd$loop[rev(seq_len(length(bd$loop)))]
            if (isTRUE(weighted == TRUE) == FALSE) {
                NA
            }
            else {
                vecol <- vecol[rev(seq_len(length(vecol)))]
                netdrpl <- netdrp[, , rev(seq_len(dim(netdrp)[3]))]
            }
        }
        else {
            bdlp <- bd$loop
            ifelse(isTRUE(weighted == TRUE) == TRUE, netdrpl <- netdrp, 
                NA)
        }
        dz <- (rng(z) + abs(min(rng(z))))/(5L)
        ndss <- nds
        ndss[, 1] <- ndss[, 1] * scl[1]
        ndss[, 2] <- ndss[, 2] * scl[2]
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
                  ifelse(isTRUE(weighted == TRUE) == TRUE, hc(ndss[lp[i], 
                    1], ndss[lp[i], 2] + (dcx), lpsz, col = vecol, 
                    lty = lty, lwd = diag(netdrpl)[lp[i]]), hc(ndss[lp[i], 
                    1], ndss[lp[i], 2] + (dcx), lpsz, col = vecol, 
                    lty = lty, lwd = lwd))
                }
                rm(i)
            }
            else {
                NA
            }
        }
        else if (isTRUE(z > 1) == TRUE) {
            ifelse(isTRUE(bwd == 0) == TRUE, dz <- rep(0, z), 
                NA)
            for (k in seq_len(length(bdlp))) {
                lp <- as.numeric(unique(multiplex::dhc(bdlp)[k][[1]]))
                if (isTRUE(length(lp) > 0) == TRUE) {
                  for (i in seq_len(length(lp))) {
                    ifelse(isTRUE(cex[lp[i]] <= 3L) == TRUE | 
                      isTRUE(n < 3) == TRUE, dz <- dz * 0.75, 
                      NA)
                    if (isTRUE(n < 3) == TRUE) {
                      dcx <- cex[lp[i]]/110L
                      lpsz <- abs((cex[lp[i]] * 0.007) - dz[k])
                    }
                    else {
                      dcx <- cex[lp[i]]/100L
                      lpsz <- abs((cex[lp[i]] * 0.0075) - dz[k])
                    }
                    ifelse(isTRUE(weighted == TRUE) == TRUE, 
                      hc(ndss[lp[i], 1], ndss[lp[i], 2] + (dcx), 
                        lpsz, col = grDevices::adjustcolor(vecol[k], 
                          alpha = alfa), lty = lty[k], lwd = netdrpl[i, 
                          i, k]), hc(ndss[lp[i], 1], ndss[lp[i], 
                        2] + (dcx), lpsz, col = grDevices::adjustcolor(vecol[k], 
                        alpha = alfa), lty = lty[k], lwd = lwd[k]))
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
    if (isTRUE(weights == TRUE) == TRUE) {
        mdps <- mdp[, 1:2]
        mdps[, 1] <- mdps[, 1] * scl[1]
        mdps[, 2] <- mdps[, 2] * scl[2]
        ifelse(missing(tcex2) == TRUE, tcex2 <- tcex, NA)
        ifelse(missing(cex2) == TRUE, cex2 <- tcex2, NA)
        ifelse(missing(tcol2) == TRUE, tcol2 <- "#000000", NA)
        if (isTRUE(flgbg == 1) == TRUE) {
            graphics::points(mdps[, 1], mdps[, 2], pch = 15, 
                cex = tcex2 * 2, col = bg, bg = bg)
        }
        else if (isTRUE(flgbg == 0) == TRUE) {
            graphics::points(mdps[, 1], mdps[, 2], pch = 15, 
                cex = tcex2 * 2, col = "#FFFFFF", bg = obg)
        }
        else {
            NA
        }
        for (k in seq_along(tlbs)) {
            graphics::text(mdps[k, ], labels = mdp[k, 3], cex = tcex2, 
                col = tcol2)
        }
        rm(k)
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
        if (isTRUE(length(pos) == 1) == TRUE) {
            if (isTRUE(pos == 0) == TRUE) {
                graphics::text(ndss, labels = lbs, cex = tcex, 
                  adj = 0.5, col = tcol)
            }
            else {
                graphics::text(ndss, lbs, cex = tcex, pos = pos, 
                  col = tcol, offset = (cex/4L), adj = c(0.5, 
                    1))
            }
        }
        else if (isTRUE(length(pos) == n) == TRUE) {
            graphics::text(ndss, lbs, cex = tcex, pos = pos, 
                col = tcol[1], offset = (cex/4L), adj = c(0.5, 
                  1))
        }
        else {
            if (isTRUE(pos[1] == 0) == TRUE) {
                graphics::text(ndss, labels = lbs, cex = tcex, 
                  adj = 0.5, col = tcol)
            }
            else {
                graphics::text(ndss, lbs, cex = tcex, pos = pos[1], 
                  col = tcol, offset = (cex/4L), adj = c(0.5, 
                    1))
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
            graphics::text(ndss, labels = atts, cex = tcex, pos = pos%%4 + 
                1L, col = tcol, offset = (cex/4L), adj = c(0.5, 
                1))
        }
        else if (isTRUE(flgcx == TRUE) == TRUE) {
            graphics::text(ndss, labels = atts, cex = tcex, pos = pos%%4 + 
                1L, col = tcol, offset = (min(cex)/4L), adj = c(0.5, 
                1))
        }
    }
    graphics::par(mar = omr)
    graphics::par(bg = obg)
    graphics::par(lend = 0)
    graphics::par(mai = rep(0, 4))
}
