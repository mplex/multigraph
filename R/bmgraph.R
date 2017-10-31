bmgraph <-
function (net, layout = c("bip", "bip3", "bip3e", "bip4", "force", 
    "rand", "circ", "stress", "CA", "circ2"), scope, coord, tcex, 
    alpha = c(1, 1, 1), showLbs, showAtts, att = NULL, lbat = "1", 
    main = NULL, cex.main, bg, mar, directed, weighted, collRecip, 
    cex, pos, lwd, lty, ecol, vcol, vcol0, asp, seed = NULL, 
    maxiter = 100, bwd, clu, pch, tcol, rot, mirrorX, mirrorY, 
    col, hds, vedist, jitter, add, ...) 
{
    ifelse(is.data.frame(net) == TRUE, net <- as.matrix(net), 
        NA)
    if (isTRUE(is.list(net) == TRUE) == TRUE) {
        net <- multiplex::transf(net, type = "toarray", lb2lb = TRUE)
    }
    else if (isTRUE(is.vector(net) == TRUE) == TRUE) {
        ifelse(missing(add) == FALSE && isTRUE(is.list(add) == 
            TRUE) == TRUE, net <- multiplex::transf(net, type = "toarray", 
            lb2lb = TRUE, lbs = sort(unique(multiplex::dhc(c(net, 
                add))))), net <- multiplex::trnf(net, tolist = FALSE, 
            lb2lb = TRUE))
    }
    else {
        NA
    }
    ifelse(missing(weighted) == FALSE && isTRUE(weighted == TRUE) == 
        TRUE, weighted <- TRUE, weighted <- FALSE)
    ifelse(missing(collRecip) == FALSE && isTRUE(collRecip == 
        FALSE) == TRUE, collRecip <- FALSE, collRecip <- TRUE)
    ifelse(missing(showLbs) == FALSE && isTRUE(showLbs == FALSE) == 
        TRUE, showLbs <- FALSE, showLbs <- TRUE)
    ifelse(missing(showAtts) == FALSE && isTRUE(showAtts == FALSE) == 
        TRUE, showAtts <- FALSE, showAtts <- TRUE)
    ifelse(missing(directed) == FALSE && isTRUE(directed == TRUE) == 
        TRUE, directed <- TRUE, directed <- FALSE)
    if (missing(scope) == FALSE) {
        if (isTRUE(is.list(scope) == TRUE) == FALSE) 
            stop("\"scope\" should be a list or a vector of lists.")
        scope <- list(scope)
        ifelse(is.null(scope[[1]]) == TRUE, scope <- scope[2:length(scope)], 
            NA)
        if (isTRUE(length(scope) > 1) == TRUE && isTRUE(names(scope[1]) == 
            "coord") == TRUE) {
            scope <- scope[length(scope):1]
            flgrev <- TRUE
        }
        else {
            flgrev <- FALSE
        }
        tmp <- scope[[1]]
        if (isTRUE(length(scope) > 1) == TRUE && isTRUE(length(scope[[1]]) > 
            1) == TRUE) {
            for (k in 2:length(scope)) {
                tmp[length(tmp) + 1L] <- as.list(scope[k])
                names(tmp)[length(tmp)] <- attr(scope[k], "names")
            }
            rm(k)
        }
        else if (isTRUE(length(scope) > 1) == TRUE) {
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
        ifelse(isTRUE(flgrev == TRUE) == TRUE, scope <- tmp[length(tmp):1], 
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
    ifelse(isTRUE(dim(net)[3] == 1) == TRUE, net <- net[, , 1], 
        NA)
    nn <- dim(net)[1]
    mm <- dim(net)[2]
    ifelse(isTRUE(is.na(dim(net)[3]) == TRUE) == TRUE, z <- 1, 
        z <- dim(net)[3])
    ifelse(isTRUE(weighted == FALSE) == TRUE, net <- multiplex::dichot(net, 
        c = 1L), NA)
    ifelse(missing(tcol) == TRUE, tcol <- c(1, 1), NA)
    ifelse(missing(bwd) == TRUE, bwd <- 1, NA)
    ifelse(isTRUE(bwd < 0L) == TRUE, bwd <- 0L, NA)
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
    if (isTRUE(z == 1) == TRUE) {
        if (isTRUE(directed == TRUE) == FALSE) {
            tnet <- t(net)
            dimnames(tnet)[[2]] <- dimnames(net)[[1]]
            bmnet <- (rbind(cbind(matrix(0, ncol = nn, nrow = nn, 
                dimnames = list(rownames(net), rownames(net))), 
                net), cbind(tnet, matrix(0, ncol = mm, nrow = mm, 
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
        for (k in seq_len(dim(net)[3])) {
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
        for (i in seq_len(nrow(bmnetdf))) {
            bmnet[, , i][seq_len((nn + mm)^2)] <- as.numeric(bmnetdf[i, 
                ])
        }
        rm(i)
    }
    n <- dim(bmnet)[1]
    ifelse(is.null(dimnames(net)[[1]]) == TRUE, nlbs <- as.character(seq_len(nn)), 
        nlbs <- dimnames(net)[[1]])
    ifelse(is.null(dimnames(net)[[2]]) == TRUE, mlbs <- as.character(nn + 
        seq_len(mm)), mlbs <- dimnames(net)[[2]])
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
        else {
            NA
        }
    }
    else {
        hds <- 1L
    }
    ifelse(missing(vedist) == TRUE, vedist <- 0, NA)
    ifelse(isTRUE(vedist > 1L) == TRUE, vedist <- 1L, NA)
    ifelse(missing(rot) == TRUE, NA, rot <- rot * -1)
    if (isTRUE(directed == FALSE) == TRUE) {
        if (isTRUE(z == 1L) == TRUE) {
            ucbmnet <- bmnet
            ucbmnet[upper.tri(ucbmnet)] <- 0L
        }
        else {
            if (isTRUE(weighted == FALSE) == TRUE) {
                ucbmnet <- multiplex::dichot(bmnet, c = 1L)
                for (k in seq_len(dim(bmnet)[3])) {
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
            else {
                ucbmnet <- bmnet
                if (isTRUE(collRecip == TRUE) == TRUE) {
                  for (k in seq_len(dim(bmnet)[3])) {
                    ucbmnet[upper.tri(ucbmnet[, , k])] <- 0L
                  }
                  rm(k)
                }
            }
        }
    }
    else {
        ucbmnet <- bmnet
    }
    if (isTRUE(directed == FALSE) == TRUE) {
        ifelse(isTRUE(z == 1L) == TRUE, bd <- multiplex::bundles(as.array(as.matrix(ucbmnet)), 
            loops = FALSE, lb2lb = FALSE, collapse = FALSE), 
            bd <- multiplex::bundles((ucbmnet), loops = FALSE, 
                lb2lb = FALSE, collapse = FALSE))
    }
    else {
        ifelse(isTRUE(z == 1L) == TRUE, bd <- multiplex::bundles(as.array(as.matrix(bmnet)), 
            loops = FALSE, lb2lb = FALSE, collapse = FALSE), 
            bd <- multiplex::bundles((bmnet), loops = FALSE, 
                lb2lb = FALSE, collapse = FALSE))
    }
    ifelse(isTRUE(z == 1L) == TRUE, r <- 1L, r <- length(bd[[1]]))
    bds <- multiplex::summaryBundles(bd, byties = TRUE)
    ifelse(isTRUE(length(bds) == 0) == TRUE, showAtts <- FALSE, 
        NA)
    ifelse(missing(ecol) == TRUE, ecol <- grDevices::gray.colors(r, 
        start = 0.1, end = 0.5), NA)
    if (isTRUE(weighted == TRUE) == TRUE) {
        ifelse(missing(lty) == TRUE, lty <- rep(1, r), lty <- rep(lty, 
            r))
    }
    else {
        ifelse(missing(lty) == TRUE, lty <- seq_len(r), NA)
    }
    if (isTRUE(z == 1L) == TRUE) {
        Lt <- lty[1]
        vecol <- ecol[1]
    }
    else {
        ifelse(isTRUE(length(lty) == 1L) == TRUE, Lt <- seq_len(r), 
            Lt <- rep(lty, r)[seq_len(r)])
        ifelse(isTRUE(length(ecol) == 1L) == TRUE, vecol <- rep(ecol, 
            z), vecol <- rep(ecol, z)[seq_len(z)])
        if (isTRUE(length(lty) == length(Lt)) == FALSE) {
            Ltc <- seq_along(vecol)
        }
        else {
            ifelse(isTRUE(seq(lty) == lty) == TRUE, Ltc <- Lt, 
                Ltc <- seq_len(r))
        }
    }
    vltz <- Lt
    if (isTRUE(weighted == TRUE) == TRUE) {
        NA
    }
    else {
        ifelse(isTRUE(bwd > 1L) == TRUE, bwd <- 1L, NA)
    }
    if (!(missing(clu))) {
        flgclu <- TRUE
        if (is.vector(clu) == FALSE) 
            stop("'clu' must be a vector")
        if (is.character(clu) == TRUE) {
            tmpclu <- clu
            for (i in seq_len(nlevels(factor(clu)))) {
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
        if (isTRUE(match.arg(layout) == "force") == TRUE | isTRUE(match.arg(layout) == 
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
    if (missing(tcex) == TRUE) {
        ifelse(isTRUE(max(cex) < 2) == TRUE, tcex <- cex * 0.66, 
            tcex <- cex * 0.33)
    }
    else {
        NA
    }
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
            for (i in seq_len(nclu)) {
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
                for (i in seq_len(nclu)) {
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
    ifelse(isTRUE(directed == TRUE) == TRUE, fds <- 90L, fds <- 100L)
    if (missing(coord) == FALSE) {
        if (isTRUE(nrow(coord) == n) == FALSE) 
            stop("Length of 'coord' does not match network order.")
        flgcrd <- TRUE
        crd <- coord
    }
    else if (missing(coord) == TRUE) {
        flgcrd <- FALSE
        switch(match.arg(layout), force = {
            crd <- frcd(as.matrix(multiplex::mnplx(bmnet)), seed = seed, 
                maxiter = maxiter)
            fds <- fds + 20L
        }, bip = {
            act <- nrm(rng(nn))
            evt <- nrm(rng(mm))
            Act <- cbind(rep(0, nrow(net)), act)
            Evt <- cbind(rep(1, ncol(net)), evt)
            crd <- rbind(Act, Evt)
            crd[which(is.nan(crd))] <- 0.5
            crd[, 2] <- crd[, 2] * cos(pi) - crd[, 1] * sin(pi)
            rownames(crd) <- lbs
            fds <- fds - 30L
        }, bip3 = {
            act1 <- nrm(rng(ceiling(nn/2)))
            act2 <- nrm(rng(floor(nn/2)))
            evt <- nrm(rng(mm))
            Act1 <- cbind(rep(0, ceiling(nrow(net)/2)), act1)
            Act2 <- cbind(rep(2, floor(nrow(net)/2)), act2)
            Evt <- cbind(rep(1, ncol(net)), evt)
            crd <- rbind(Act1, Act2, Evt)
            crd[which(is.nan(crd))] <- 0.5
            crd[, 2] <- crd[, 2] * cos(pi) - crd[, 1] * sin(pi)
            rownames(crd) <- lbs
        }, bip3e = {
            act <- nrm(rng(nn))
            Act <- cbind(rep(1, nrow(net)), act)
            evt1 <- nrm(rng(ceiling(mm/2)))
            evt2 <- nrm(rng(floor(mm/2)))
            Evt1 <- cbind(rep(0, ceiling(ncol(net)/2)), evt1)
            Evt2 <- cbind(rep(2, floor(ncol(net)/2)), evt2)
            crd <- rbind(Act, Evt1, Evt2)
            crd[which(is.nan(crd))] <- 0.5
            crd[, 2] <- crd[, 2] * cos(pi) - crd[, 1] * sin(pi)
            rownames(crd) <- lbs
        }, bip4 = {
            qdn <- nn/(nn + 3)
            qdm <- mm/(mm + 3)
            act1 <- nrm(rng(ceiling(nn/2)))
            act2 <- nrm(rng(floor(nn/2)))
            evt1 <- nrm(rng(ceiling(mm/2)))
            evt2 <- nrm(rng(floor(mm/2)))
            Act2 <- cbind((act2 * qdn) + (0.1/qdn), rep(0, floor(nn/2)))
            Evt2 <- cbind((evt2 * qdm) + (0.1/qdm), rep(1, floor(mm/2)))
            Act1 <- cbind(rep(0, ceiling(nn/2)), (act1 * qdn) + 
                (0.1/qdn))
            Evt1 <- cbind(rep(1, ceiling(mm/2)), (evt1 * qdm) + 
                (0.1/qdm))
            crd <- as.data.frame(rbind(Act1, Act2, Evt1, Evt2))
            crd[which(is.na(crd))] <- 0.5
            crd[, 2] <- crd[, 2] * cos(pi) - crd[, 1] * sin(pi)
        }, rand = {
            set.seed(seed)
            crd <- data.frame(X = round(stats::runif(n) * 1L, 
                5), Y = round(stats::runif(n) * 1L, 5))
        }, circ = {
            crd <- data.frame(X = sin(2L * pi * ((0:(n - 1L))/n)), 
                Y = cos(2L * pi * ((0:(n - 1L))/n)))
            ifelse(isTRUE(flgcx == TRUE) == TRUE, fds <- fds - 
                10L, NA)
        }, stress = {
            crd <- stsm(as.matrix(multiplex::mnplx(bmnet)), seed = seed, 
                maxiter = maxiter, ...)
        }, CA = {
            ifelse(is.na(dim(net)[3]) == FALSE, nt <- as.matrix(multiplex::mnplx(net)), 
                nt <- net)
            Wn <- drop(as.matrix(nt) %*% (rep(1/sum(nt), ncol(nt))))
            Wm <- drop((rep(1/sum(nt), nn)) %*% as.matrix(nt))
            outWW <- t(t((as.matrix(nt)/sum(nt) - outer(Wn, Wm)) * 
                1/sqrt(Wn)) * 1/sqrt(Wm))
            outWW[which(is.nan(outWW))] <- mean(outWW, na.rm = TRUE)
            M <- svd(outWW)
            crd <- rbind(M$u[, 1L:2L] * 1/sqrt(Wn), M$v[, 1L:2L] * 
                1/sqrt(Wm))
            crd[which(crd[, 1] == -Inf), 1] <- min(crd[, 1])
            crd[which(crd[, 2] == -Inf), 2] <- min(crd[, 2])
            crd[which(crd[, 1] == Inf), 1] <- max(crd[, 1])
            crd[which(crd[, 2] == Inf), 2] <- max(crd[, 2])
            crd[which(is.nan(crd))] <- 0L
            ifelse(missing(jitter) == FALSE, crd <- jitter(crd, 
                amount = jitter), NA)
        }, circ2 = {
            act <- seq(1.25, pi * 1.2, len = nn)
            evt <- seq(min(act), max(act), len = mm)
            act1 <- rbind(cos(act), sin(act))
            evt1 <- rbind(cos(evt), sin(evt))
            evt1[1, ] <- evt1[1, ] * cos(pi) - evt1[2, ] * sin(pi)
            evt1[2, ] <- evt1[2, ] * cos(pi) - evt1[1, ] * sin(pi)
            crd <- as.data.frame(t(cbind(act1, evt1)))
            crd[which(is.na(crd))] <- 0.5
            rt <- 10
            crd[, 1] <- (crd[, 1] * cos(rt * (pi/180L)) - crd[, 
                2] * sin(rt * (pi/180L)))/0.4
            crd[, 2] <- (crd[, 2] * cos(rt * (pi/180L)) + crd[, 
                1] * sin(rt * (pi/180L)))/1
            crd[, 1:2] <- crd[, 1:2] - min(crd[, 1:2])
            ifelse(isTRUE(flgcx == TRUE) == TRUE, fds <- fds - 
                50L, NA)
        })
    }
    if (missing(rot) == FALSE) {
        crd[, 1:2] <- xyrt(crd[, 1:2], as.numeric(rot))
        crd[, 1:2] <- crd[, 1:2] - min(crd[, 1:2])
    }
    else {
        NA
    }
    if (missing(pos) == TRUE) {
        flgpos <- TRUE
        ifelse(match.arg(layout) == "bip3" | match.arg(layout) == 
            "bip3e", pos <- c(2, 1, 4), pos <- c(2, 4))
        ifelse(match.arg(layout) == "bip4", pos <- c(2, 4, 3, 
            1), NA)
        ifelse(match.arg(layout) == "circ" | match.arg(layout) == 
            "stress", pos <- c(4, 2), NA)
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
    rat <- (max(crd[, 1]) - min(crd[, 1]))/(max(crd[, 2]) - min(crd[, 
        2]))
    crd[, 1] <- (crd[, 1] - min(crd[, 1]))/(max(crd[, 1]) - min(crd[, 
        1]))
    ifelse(isTRUE(rat > 0) == TRUE, crd[, 2] <- ((crd[, 2] - 
        min(crd[, 2]))/(max(crd[, 2]) - min(crd[, 2]))) * (1L/rat), 
        crd[, 2] <- ((crd[, 2] - min(crd[, 2]))/(max(crd[, 2]) - 
            min(crd[, 2]))) * (rat))
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
    nds <- (2L/max(nds)) * nds
    ifelse(match.arg(layout) == "stress", nds[, 1] <- nds[, 1] * 
        1L, NA)
    if (isTRUE(flgcx == TRUE) == TRUE && isTRUE(sqrt(((max(nds[, 
        1]) - min(nds[, 1])) * (max(nds[, 2]) - min(nds[, 2])))/nrow(nds)) < 
        (1/3)) == TRUE) {
        nds <- nds * (2.223 - (4.45 * (sqrt(((max(nds[, 1]) - 
            min(nds[, 1])) * (max(nds[, 2]) - min(nds[, 2])))/n))))
    }
    else {
        nds <- nds * (0.5)
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
    if (match.arg(layout) == "force" | match.arg(layout) == "rand" | 
        match.arg(layout) == "circ" | match.arg(layout) == "stress") {
        xlim <- c(min(nds[, 1]) - (max(cex)/50L) - (0), max(nds[, 
            1]) + (max(cex)/50L) + (0))
        ylim <- c(min(nds[, 2]) - (max(cex)/100L), max(nds[, 
            2]) + (max(cex)/100L))
        ifelse(missing(asp) == TRUE, asp <- 1, NA)
        fds <- fds + 20L
    }
    else if (match.arg(layout) == "bip4") {
        xlim <- c(min(nds[, 1]) - (max(cex)/100L) - (0), max(nds[, 
            1]) + (max(cex)/100L) + (0))
        ylim <- c(min(nds[, 2]) - (max(cex)/100L), max(nds[, 
            2]) + (max(cex)/100L))
        ifelse(missing(asp) == TRUE, asp <- 1.2, NA)
    }
    else {
        mx <- ceiling(length(multiplex::dhc(lbs[seq_len(nn)], 
            ""))/nn) * 0.02
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
    for (i in seq_len(length(attr(bds, "names")))) {
        ifelse(isTRUE(length(multiplex::dhc(attr(bds, "names")[i], 
            sep = "")) > 4L) == TRUE, tlbs <- append(tlbs, tolower(paste(multiplex::dhc(attr(bds, 
            "names")[i], sep = "")[1:4], collapse = ""))), tlbs <- append(tlbs, 
            tolower(attr(bds, "names"))[i]))
    }
    rm(i)
    ifelse(isTRUE(weighted == TRUE) == TRUE && isTRUE(max(ucbmnet) > 
        10L) == TRUE, fnucbmnet <- (norm(as.matrix(ucbmnet), 
        type = "F")), NA)
    if (isTRUE(length(bds) > 0) == TRUE) {
        for (k in seq_len(length(bds))) {
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
                  if (isTRUE(directed == TRUE) == TRUE && isTRUE(tlbs[k] %in% 
                    c("asym", "tent", "txch", "mixd")) == TRUE) {
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
                  cx <- rep(cex[1], 2)
                }
                ifelse(missing(vedist) == TRUE, vedist <- 0, 
                  NA)
                if ((match.arg(layout) == "bip" | match.arg(layout) == 
                  "bip3" | match.arg(layout) == "bip3e")) {
                  ifelse(isTRUE(vedist < 0) == TRUE, fds <- fds + 
                    (vedist * 0.2), fds <- fds + (vedist * -0.2))
                }
                else if (match.arg(layout) == "bip4") {
                  ifelse(isTRUE(vedist < 0) == TRUE, fds <- fds - 
                    (vedist * 0.2), fds <- fds - (vedist * -0.2))
                }
                else {
                  NA
                }
                if (isTRUE(weighted == TRUE) == TRUE) {
                  Lw <- vector()
                  i <- 1
                  for (j in seq_len(length(bds[[k]]))) {
                    qn <- c(prs[i], prs[(i + 1)])
                    if (isTRUE(directed == FALSE) == TRUE && 
                      isTRUE(collRecip == TRUE) == TRUE) {
                      ifelse(isTRUE(z == 1L) == TRUE, Lw <- append(Lw, 
                        ucbmnet[qn[1], qn[2]]), Lw <- append(Lw, 
                        ucbmnet[qn[1], qn[2], vltc[j]] + t(ucbmnet[qn[1], 
                          qn[2], vltc[j]])))
                    }
                    else if (isTRUE(directed == TRUE) == TRUE | 
                      isTRUE(collRecip == FALSE) == TRUE) {
                      ifelse(isTRUE(z == 1L) == TRUE, Lw <- append(Lw, 
                        ucbmnet[qn[1], qn[2]] + t(ucbmnet)[qn[1], 
                          qn[2]]), Lw <- append(Lw, ucbmnet[qn[1], 
                        qn[2], vltc[j]]))
                    }
                    i <- i + 2L
                  }
                  rm(j)
                  rm(i)
                  if (isTRUE(max(ucbmnet) > 10L) == TRUE) {
                    lwd <- (Lw/fnucbmnet) * (10L * 5L)
                  }
                  else {
                    lwd <- Lw
                  }
                }
                else {
                  lwd <- rep(lwd[1], rr)
                }
                flgcr <- rep(0L, z)
                if (isTRUE(z == 1L) == TRUE) {
                  mbnd(pars, rr, bds[[k]], vlt, cx, lwd, vecol, 
                    directed, bwd, alfa, fds, flgcx, weighted, 
                    flgcr, hds, n)
                }
                else {
                  ifelse(isTRUE(length(lty) == 1L) == TRUE, mbnd(pars, 
                    rr, bds[[k]], vlt1, cx, lwd, vecol[vltc], 
                    directed, bwd, alfa, fds, flgcx, weighted, 
                    flgcr, hds, n), mbnd(pars, rr, bds[[k]], 
                    vlt, cx, lwd, vecol[vltc], directed, bwd, 
                    alfa, fds, flgcx, weighted, flgcr, hds, n))
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
        ndss <- nds
        if (match.arg(layout) == "bip3") {
            if (isTRUE(pos[1] == 0) == TRUE) {
                graphics::text(ndss[seq_len(length(act1)), ], 
                  labels = lbs[seq_len(length(act1))], cex = tcex, 
                  adj = 0.5, col = tcol[1])
            }
            else {
                graphics::text(ndss[seq_len(length(act1)), ], 
                  lbs[seq_len(length(act1))], cex = tcex, pos = pos[1], 
                  col = tcol[1], offset = (cex/4L), adj = c(0.5, 
                    1))
            }
            if (isTRUE(pos[2] == 0) == TRUE) {
                graphics::text(ndss[(nn + 1L):dim(bmnet)[1], 
                  ], labels = lbs[(nn + 1L):dim(bmnet)[1]], cex = tcex, 
                  adj = 0.5, col = tcol[2])
            }
            else {
                graphics::text(ndss[(nn + 1L):dim(bmnet)[1], 
                  ], lbs[(nn + 1L):dim(bmnet)[1]], cex = tcex, 
                  pos = pos[2], col = tcol[2], offset = (cex/4L), 
                  adj = c(0.5, 1))
            }
            if (isTRUE(pos[3] == 0) == TRUE) {
                graphics::text(ndss[(length(act1) + 1):nn, ], 
                  labels = lbs[(length(act1) + 1):nn], cex = tcex, 
                  adj = 0.5, col = tcol[1])
            }
            else {
                graphics::text(ndss[(length(act1) + 1):nn, ], 
                  lbs[(length(act1) + 1):nn], cex = tcex, pos = pos[3], 
                  col = tcol[1], offset = (cex/4L), adj = c(0.5, 
                    1))
            }
        }
        else if (match.arg(layout) == "bip3e") {
            if (isTRUE(pos[2] == 0) == TRUE) {
                graphics::text(ndss[seq_len(length(act)), ], 
                  labels = lbs[seq_len(length(act))], cex = tcex, 
                  adj = 0.5, col = tcol[2])
            }
            else {
                graphics::text(ndss[seq_len(length(act)), ], 
                  lbs[seq_len(length(act))], cex = tcex, pos = pos[2], 
                  col = tcol[2], offset = (cex/4L), adj = c(0.5, 
                    1))
            }
            if (isTRUE(pos[1] == 0) == TRUE) {
                graphics::text(ndss[(nn + 1L):(nn + length(evt1)), 
                  ], labels = lbs[(nn + 1L):(nn + length(evt1))], 
                  cex = tcex, adj = 0.5, col = tcol[1])
            }
            else {
                graphics::text(ndss[(nn + 1L):(nn + length(evt1)), 
                  ], lbs[(nn + 1L):(nn + length(evt1))], cex = tcex, 
                  pos = pos[1], col = tcol[1], offset = (cex/4L), 
                  adj = c(0.5, 1))
            }
            if (isTRUE(pos[3] == 0) == TRUE) {
                graphics::text(ndss[(nn + 1 + length(evt1)):dim(bmnet)[1], 
                  ], labels = lbs[(nn + 1 + length(evt1)):dim(bmnet)[1]], 
                  cex = tcex, adj = 0.5, col = tcol[1])
            }
            else {
                graphics::text(ndss[(nn + 1 + length(evt1)):dim(bmnet)[1], 
                  ], lbs[(nn + 1 + length(evt1)):dim(bmnet)[1]], 
                  cex = tcex, pos = pos[3], col = tcol[1], offset = (cex/4L), 
                  adj = c(0.5, 1))
            }
        }
        else if (match.arg(layout) == "bip4") {
            if (isTRUE(pos[1] == 0) == TRUE) {
                graphics::text(ndss[seq_len(length(act1)), ], 
                  labels = lbs[seq_len(length(act1))], cex = tcex, 
                  adj = 0.5, col = tcol[1])
                graphics::text(ndss[(length(act1) + 1):nn, ], 
                  labels = lbs[(length(act1) + 1):nn], cex = tcex, 
                  adj = 0.5, col = tcol[1])
            }
            else {
                graphics::text(ndss[seq_len(length(act1)), ], 
                  lbs[seq_len(length(act1))], cex = tcex, pos = pos[1], 
                  col = tcol[1], offset = (cex/4L), adj = c(0.5, 
                    1))
                graphics::text(ndss[(length(act1) + 1):nn, ], 
                  labels = lbs[(length(act1) + 1):nn], cex = tcex, 
                  pos = pos[3], col = tcol[1], offset = (cex/4L), 
                  adj = c(0.5, 1))
            }
            if (isTRUE(pos[2] == 0) == TRUE) {
                graphics::text(ndss[(nn + 1L):(nn + length(evt1)), 
                  ], labels = lbs[(nn + 1L):(nn + length(evt1))], 
                  cex = tcex, adj = 0.5, col = tcol[2])
                graphics::text(ndss[(nn + 1 + length(evt1)):dim(bmnet)[1], 
                  ], labels = lbs[(nn + 1 + length(evt1)):dim(bmnet)[1]], 
                  cex = tcex, adj = 0.5, col = tcol[2])
            }
            else {
                graphics::text(ndss[(nn + 1L):(nn + length(evt1)), 
                  ], labels = lbs[(nn + 1L):(nn + length(evt1))], 
                  cex = tcex, pos = pos[2], col = tcol[2], offset = (cex/4L), 
                  adj = c(0.5, 1))
                graphics::text(ndss[(nn + 1 + length(evt1)):dim(bmnet)[1], 
                  ], lbs[(nn + 1 + length(evt1)):dim(bmnet)[1]], 
                  cex = tcex, pos = pos[4], col = tcol[2], offset = (cex/4L), 
                  adj = c(0.5, 1))
            }
        }
        else {
            if (isTRUE(length(pos) == 2) == TRUE) {
                if (isTRUE(pos[1] == 0) == TRUE) {
                  graphics::text(ndss[seq_len(nn), ], labels = lbs[seq_len(nn)], 
                    cex = tcex, adj = 0.5, col = tcol[1])
                }
                else {
                  graphics::text(ndss[seq_len(nn), ], lbs[seq_len(nn)], 
                    cex = tcex, pos = pos[1], col = tcol[1], 
                    offset = (cex/4L), adj = c(0.5, 1))
                }
                if (isTRUE(pos[2] == 0) == TRUE) {
                  graphics::text(ndss[(nn + 1L):dim(bmnet)[1], 
                    ], labels = lbs[(nn + 1L):dim(bmnet)[1]], 
                    cex = tcex, adj = 0.5, col = tcol[2])
                }
                else {
                  graphics::text(ndss[(nn + 1L):dim(bmnet)[1], 
                    ], lbs[(nn + 1L):dim(bmnet)[1]], cex = tcex, 
                    pos = pos[2], col = tcol[2], offset = (cex/4L), 
                    adj = c(0.5, 1))
                }
            }
            else if (isTRUE(length(pos) == 1) == TRUE) {
                if (isTRUE(pos == 0) == TRUE) {
                  graphics::text(ndss, labels = lbs, cex = tcex, 
                    adj = 0.5, col = tcol[1])
                }
                else {
                  graphics::text(ndss, lbs, cex = tcex, pos = pos, 
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
    graphics::par(lend = 0)
    x <- NULL
    rm(x)
}
