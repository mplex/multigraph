## (GPL-3)  Antonio Rivero Ostoic, multiplex@post.com
multigraph <-
function (net, layout = c("circ", "force", "stress", "conc", 
    "rand"), scope, directed = TRUE, loops, signed, valued, values, 
    lbs, showLbs, att, lbat, showAtts, main = NULL, cex.main, 
    col.main, font.main, coord, collRecip, undRecip, seed = NULL, 
    maxiter = 100, clu, cex, cex2, pch, lwd, lty, vcol, vcol0, 
    col, ecol, bwd, bwd2, pos, bg, bg2, asp, drp, add, swp, swp2, 
    flip, alpha = c(1, 1, 1, 1), rot, mirrorX, mirrorY, mirrorD, 
    mirrorL, mirrorV, mirrorH, scl, hds, vedist, vldist, mar, 
    ffamily, fstyle, fsize, fsize2, fcol, fcol2, lclu, sel, new, 
    mai, lscl, rm.isol, nr, jitter, ...) 
{
    cnet <- attr(net, "class")
    flgmlvl <- FALSE
    flgcn2 <- FALSE
    flgpf <- FALSE
    ifelse(missing(valued) == FALSE && isTRUE(valued == TRUE) == 
        TRUE, valued <- TRUE, valued <- FALSE)
    ifelse(missing(values) == FALSE && isTRUE(values == TRUE) == 
        TRUE, values <- TRUE, values <- FALSE)
    ifelse(missing(showAtts) == FALSE && isTRUE(showAtts == FALSE) == 
        TRUE, showAtts <- FALSE, showAtts <- TRUE)
    if (isTRUE("Rel.System" %in% cnet) == TRUE) {
        net <- net$Ties
    }
    else if (isTRUE("Multilevel" %in% cnet) == TRUE) {
        mlv <- net
        if (isTRUE("bpn" %in% cnet) == TRUE) {
            flgmlvl <- TRUE
        }
        else if (isTRUE("cn2" %in% cnet) == TRUE) {
            flgcn2 <- TRUE
        }
        else {
            NA
        }
        net <- mlv$mlnet
    }
    else if (isTRUE("pathfinder" %in% cnet) == TRUE) {
        flgpf <- TRUE
        maxw <- net$max
        net <- net$Q
    }
    else if (isTRUE(tolower(cnet) == "dataset") == TRUE || any(tolower(cnet) %in% 
        "data.set") == TRUE) {
        if (is.array(net$net) == TRUE) {
            if (isTRUE(showAtts == TRUE) == TRUE) {
                att <- net$net[, , which(net$atnet[[1]] == 1)]
                net <- net$net[, , which(net$atnet[[1]] == 0)]
            }
            else {
                net <- net$net
            }
        }
        else if (is.data.frame(net$net) == TRUE) {
            net <- multiplex::edgel(net$net, toarray = TRUE, 
                header = TRUE, valued = valued, sep = ",")
        }
        else {
            NA
        }
    }
    else if (is.vector(net) == TRUE && is.list(net) == FALSE) {
        ifelse(missing(add) == FALSE && isTRUE(is.vector(add) == 
            TRUE) == TRUE, net <- multiplex::transf(net, type = "toarray", 
            lb2lb = TRUE, lbs = sort(unique(multiplex::dhc(c(unlist(net), 
                as.character(add)))))), net <- multiplex::trnf(net, 
            tolist = FALSE, lb2lb = TRUE))
    }
    else if (is.data.frame(net) == TRUE) {
        net <- multiplex::edgel(net, toarray = TRUE, header = TRUE, 
            valued = valued, sep = ",")
    }
    else {
        NA
    }
    ifelse(missing(add) == FALSE && is.vector(add) == TRUE, flgad <- TRUE, 
        flgad <- FALSE)
    if (isTRUE(is.data.frame(net) == TRUE) == FALSE) {
        if (is.list(net) == TRUE && isTRUE(cnet == "Signed") == 
            FALSE) {
            if (isTRUE("data.set" %in% tolower(class(net))) == 
                TRUE) {
                if (isTRUE(typeof(net) == "list") == TRUE) {
                  net <- lapply(net, as.vector)
                }
                else if (isTRUE(typeof(net) == "character") == 
                  TRUE) {
                  net <- as.vector(net)
                }
            }
            else {
                NA
            }
            if (isTRUE(length(unique(unlist(lapply(net, dim)))) > 
                1L) == TRUE) 
                stop(paste(sQuote("net"), " should be a square array.", 
                  sep = ""), call. = FALSE)
            ifelse(isTRUE(flgad == TRUE) == TRUE, net <- multiplex::transf(net, 
                type = "toarray", lb2lb = TRUE, lbs = sort(unique(unlist(multiplex::dhc(c(unlist(net), 
                  as.character(add))))))), net <- multiplex::transf(net, 
                type = "toarray", lb2lb = TRUE))
        }
        else if (is.null(net) == TRUE && isTRUE(flgad == TRUE) == 
            TRUE) {
            net <- matrix(0L, nrow = length(add), ncol = length(add), 
                dimnames = list(add, add))
        }
        else {
            ifelse(isTRUE(flgad == TRUE) == TRUE, net <- multiplex::transf(net, 
                type = "toarray", lbs = c(dimnames(net)[[1]], 
                  as.character(add))), NA)
        }
    }
    else {
        net <- as.matrix(net)
    }
    if (isTRUE(cnet == "Signed") == FALSE) {
        ifelse(is.array(net) == TRUE || is.matrix(net) == TRUE, 
            NA, stop(paste(sQuote("net"), " should be a matrix or array.", 
                sep = ""), call. = FALSE))
    }
    ifelse(isTRUE(dim(net)[3] == 1) == TRUE, net <- net[, , 1], 
        NA)
    ifelse(missing(undRecip) == FALSE && isTRUE(undRecip == TRUE) == 
        TRUE, undRecip <- TRUE, undRecip <- FALSE)
    if (isTRUE(flgmlvl == TRUE) == TRUE) {
        ifelse(missing(loops) == FALSE && isTRUE(loops == FALSE) == 
            TRUE, loops <- FALSE, loops <- TRUE)
    }
    else {
        ifelse(missing(loops) == FALSE && isTRUE(loops == TRUE) == 
            TRUE, loops <- TRUE, loops <- FALSE)
    }
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
        if (isTRUE(directed == FALSE) == TRUE || isTRUE(undRecip == 
            TRUE) == TRUE) {
            ifelse(missing(collRecip) == FALSE && isTRUE(collRecip == 
                TRUE) == FALSE, collRecip <- FALSE, collRecip <- TRUE)
        }
        else if (isTRUE(directed == TRUE) == TRUE) {
            ifelse(missing(collRecip) == FALSE && isTRUE(collRecip == 
                TRUE) == TRUE, collRecip <- TRUE, collRecip <- FALSE)
        }
        else {
            ifelse(missing(collRecip) == FALSE && isTRUE(collRecip == 
                FALSE) == TRUE, collRecip <- FALSE, collRecip <- TRUE)
        }
    }
    if ((missing(showLbs) == FALSE && isTRUE(showLbs == TRUE) == 
        TRUE) || (missing(lbs) == FALSE)) {
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
    if (missing(swp) == FALSE && isTRUE(swp == TRUE) == TRUE) {
        if (isTRUE(valued == TRUE) == TRUE) {
            warning(paste("Option ", sQuote("swp"), " not supported for ", 
                sQuote("valued"), " graphs.", sep = ""), call. = FALSE)
            swp <- FALSE
        }
        else {
            swp <- TRUE
        }
    }
    else {
        swp <- FALSE
    }
    ifelse(missing(swp2) == FALSE && isTRUE(swp2 == TRUE) == 
        TRUE, swp2 <- TRUE, swp2 <- FALSE)
    ifelse(isTRUE(directed == FALSE) == TRUE, directed <- FALSE, 
        NA)
    if (missing(scope) == FALSE) {
        if (isTRUE(is.list(scope) == TRUE) == FALSE) 
            stop(paste("Values in ", sQuote("scope"), " should be in a list or a vector of lists.", 
                sep = ""), call. = FALSE)
        for (k in seq_len(length(scope))) {
            if (is.factor(scope[[k]]) == TRUE) {
                warning(paste("Vectors with class ", dQuote("factor"), 
                  " must be outside ", sQuote("scope"), " since levels are ignored.", 
                  sep = ""), call. = FALSE)
                break
            }
        }
        rm(k)
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
    ifelse(missing(lwd) == TRUE, lwd <- 1, lwd <- abs(lwd))
    if (isTRUE(flgmlvl == TRUE) == TRUE) {
        ifelse(missing(pch) == TRUE && all(mlv$lbs$cdm %in% mlv$lbs$dm) == 
            FALSE, pch <- c(rep(21, length(mlv$lbs$dm)), rep(22, 
            length(which(!(mlv$lbs$cdm %in% mlv$lbs$dm))))), 
            NA)
        ifelse(missing(pch) == TRUE, pch <- c(21, 22), NA)
        ifelse(missing(ecol) == FALSE && isTRUE(length(ecol) == 
            2) == TRUE, ecol <- c(rep(ecol[1], length(which(mlv$modes == 
            "1M"))), rep(ecol[2], length(which(mlv$modes == "2M")))), 
            NA)
        ifelse(missing(fcol) == TRUE, fcol <- c(1, 1), fcol <- rep(fcol, 
            2))
    }
    else {
        ifelse(missing(pch) == TRUE, pch <- 21, NA)
        ifelse(missing(fcol) == TRUE, fcol <- 1, NA)
    }
    ifelse(missing(bwd) == TRUE, bwd <- 1, NA)
    ifelse(isTRUE(bwd < 0L) == TRUE, bwd <- 0, NA)
    if (missing(bg) == TRUE) {
        bg <- graphics::par()$bg
        flgbg <- FALSE
    }
    else {
        obg <- graphics::par()$bg
        ifelse(isTRUE(bg == 0) == TRUE, bg <- "#FFFFFF", NA)
        flgbg <- TRUE
    }
    ifelse(missing(cex.main) == TRUE, cex.main <- graphics::par()$cex.main, 
        NA)
    ifelse(missing(col.main) == TRUE, col.main <- graphics::par()$col.main, 
        NA)
    ifelse(missing(font.main) == TRUE, font.main <- graphics::par()$font.main, 
        NA)
    if (missing(rot) == FALSE) {
        ifelse(isTRUE(rot == -90L) == TRUE, rot <- -89.99 * -1, 
            rot <- rot[1] * -1)
    }
    else {
        NA
    }
    if (isTRUE(length(alpha) < 2L) == TRUE) {
        alfa <- 1L
        alpha <- rep(alpha, 4L)
    }
    else {
        alfa <- alpha[2]
    }
    if (isTRUE(length(alpha) < 3L) == TRUE) 
        alpha <- append(alpha, c(0.1, 0.5))
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
    }
    else {
        ifelse(missing(scl) == TRUE, hds <- 1L, hds <- 1L * scl)
    }
    ifelse(isTRUE(dim(net)[1] > 8) == TRUE || isTRUE(valued == 
        TRUE) == TRUE || isTRUE(flgmlvl == TRUE) == TRUE, hds <- hds * 
        0.75, NA)
    ifelse(missing(signed) == FALSE && isTRUE(signed == TRUE) == 
        TRUE, signed <- TRUE, signed <- FALSE)
    if (isTRUE(signed == TRUE) == TRUE || isTRUE(cnet == "Signed") == 
        TRUE) {
        if (isTRUE(cnet == "Signed") == TRUE) {
            if (any(net$val %in% c(-1, 0, 1)) == TRUE) {
                net <- multiplex::zbind(multiplex::dichot(net$s, 
                  c = 1L), 1 - multiplex::dichot(net$s, c = 0))
            }
            else {
                nets <- multiplex::zbnd(net$s, net$s)
                net <- array(NA, dim = c(dim(nets)[1], dim(nets)[2], 
                  2), dimnames = list(dimnames(net$s)[[1]], dimnames(net$s)[[1]]))
                net <- replace(net, nets == "o", 0L)
                net <- replace(net, nets == "a", 1L)
                net[, , 1] <- replace(net[, , 1], nets[, , 1] == 
                  "p", 1L)
                net[, , 1] <- replace(net[, , 1], nets[, , 1] == 
                  "n", 0L)
                net[, , 2] <- replace(net[, , 2], nets[, , 2] == 
                  "n", 1L)
                net[, , 2] <- replace(net[, , 2], nets[, , 2] == 
                  "p", 0L)
                rm(nets)
            }
        }
        else {
            ifelse(isTRUE(is.na(dim(net)[3]) == TRUE) == TRUE, 
                NA, net <- net[, , seq_len(2)])
        }
    }
    else {
        net[is.na(net)] <- 0L
        net <- replace(net, net == Inf, 0L)
    }
    ifelse(isTRUE(is.na(dim(net)[3]) == TRUE) == TRUE, z <- 1L, 
        z <- dim(net)[3])
    if (missing(lscl) == FALSE) {
        if (isTRUE(valued == TRUE) == TRUE) {
            ifelse(isTRUE(lscl > 2L) == TRUE || isTRUE(lscl < 
                0) == TRUE, lscl <- 1, NA)
            if (isTRUE(max(net) == max(diag(multiplex::mnplx(net, 
                valued = TRUE)))) == TRUE) {
                if (isTRUE(z > 1L) == TRUE) {
                  for (k in z) {
                    net[which(net[, , z] == max(net[, , z]))] <- max(net[, 
                      , z]) * lscl
                  }
                  rm(k)
                }
                else {
                  net[which(net == max(net))] <- max(net) * lscl
                }
            }
            else {
                invisible(NA)
            }
        }
        else if (isTRUE(valued == FALSE) == TRUE) {
            ifelse(isTRUE(lscl > 1L) == TRUE || isTRUE(lscl < 
                0) == TRUE, lscl <- 1, NA)
        }
    }
    else {
        NA
    }
    ifelse(missing(rm.isol) == FALSE && isTRUE(rm.isol == TRUE) == 
        TRUE, rm.isol <- TRUE, rm.isol <- FALSE)
    if (isTRUE(rm.isol == TRUE) == TRUE && isTRUE(dim(net)[1] > 
        1L) == TRUE) {
        inc <- multiplex::rel.sys(net, bonds = "entire")$incl
        if (missing(clu) == FALSE && isTRUE(length(clu) == dim(net)[1]) == 
            TRUE && is.vector(clu) == TRUE) {
            clu <- clu[which(dimnames(net)[[1]] %in% inc)]
        }
        else {
            clu <- rep(1, dim(net)[1])
        }
        ifelse(missing(lbs) == FALSE, lbs <- lbs[inc], NA)
        ifelse(missing(lclu) == TRUE, lclu <- seq(0, length(unique(clu))), 
            NA)
        net <- multiplex::rm.isol(net)
    }
    else {
        ifelse(missing(clu) == FALSE, NA, clu <- rep(1, dim(net)[1]))
    }
    n <- dim(net)[1]
    if (missing(lbs) == TRUE) {
        ifelse(is.null(dimnames(net)[[1]]) == TRUE, lbs <- as.character(seq_len(dim(net)[1])), 
            lbs <- dimnames(net)[[1]])
    }
    else if (missing(lbs) == FALSE) {
        ifelse(isTRUE(flgmlvl == TRUE) == TRUE, lbs <- unlist(lbs), 
            NA)
        ifelse(isTRUE(is.character(lbs) == TRUE) == FALSE, lbs <- as.character(lbs), 
            NA)
        if (identical(lbs, make.unique(lbs)) == FALSE) {
            warning(paste("Make unique ", sQuote("lbs"), " because of duplicated entries.", 
                sep = ""), call. = FALSE)
            lbs <- make.unique(lbs)
        }
        else {
            NA
        }
        if (isTRUE(length(lbs) == dim(net)[1]) == TRUE) {
            dimnames(net)[[1]] <- dimnames(net)[[2]] <- lbs
        }
        else {
            warning(paste("Length of ", sQuote("lbs"), " not equal to number of vertices in ", 
                sQuote("net"), ".", sep = ""), call. = FALSE)
            dimnames(net)[[1]] <- dimnames(net)[[2]] <- lbs[seq_len(dim(net)[1])]
        }
    }
    else {
        NA
    }
    if (isTRUE(swp == TRUE) == TRUE && isTRUE(z > 1L) == TRUE) {
        net <- net[, , rev(seq_len(z))]
    }
    else {
        NA
    }
    if (missing(att) == FALSE && is.array(att) == TRUE) {
        if (isTRUE(n != dim(att)[1]) == TRUE && isTRUE(flgad == 
            FALSE) == TRUE) {
            warning(paste("Attributes are not shown since dimensions in ", 
                sQuote("net"), " and ", sQuote("att"), " differ.", 
                sep = ""), call. = FALSE)
            showAtts <- FALSE
        }
    }
    if (any(is.na(net)) == TRUE) {
        net[is.na(net)] <- 0
        warning(paste("Missing values in ", sQuote("net"), " are ignored.", 
            sep = ""), call. = FALSE)
    }
    if (missing(drp) == FALSE && is.numeric(drp) == TRUE) {
        netdrp <- replace(net, net <= drp, 0)
        netd <- multiplex::dichot(netdrp, c = 1L)
    }
    else {
        ifelse(isTRUE(sum(net) == 0) == TRUE, netd <- net, netd <- multiplex::dichot(net, 
            c = min(net[net > 0])))
        netdrp <- net
    }
    if (isTRUE(directed == FALSE) == TRUE && isTRUE(collRecip == 
        TRUE) == TRUE && isTRUE(valued == TRUE) == TRUE) {
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
    if (isTRUE(directed == FALSE) == TRUE && (isTRUE(collRecip == 
        TRUE) == TRUE && isTRUE(undRecip == TRUE) == FALSE) || 
        isTRUE(directed == FALSE) == TRUE && isTRUE(undRecip == 
            TRUE) == TRUE) {
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
    ifelse(isTRUE(sum(net) == 0) == TRUE && isTRUE(loops == TRUE) == 
        TRUE, bd$loop <- character(0), NA)
    bds <- multiplex::summaryBundles(bd, byties = TRUE)
    ifelse(isTRUE(z == 1L) == TRUE, r <- 1L, r <- length(bd[[1]]))
    if (missing(ecol) == TRUE) {
        ifelse(isTRUE(signed == TRUE) == TRUE, ecol <- rep(1, 
            2), ecol <- grDevices::gray.colors(r, start = 0.1, 
            end = 0.5))
    }
    else {
        ifelse(isTRUE(ecol == 0) == TRUE, ecol <- "transparent", 
            NA)
    }
    if (isTRUE(valued == TRUE) == TRUE) {
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
                ifelse(isTRUE(swp == TRUE) == TRUE, Ltc <- rev(seq_len(r)), 
                  Ltc <- seq_len(r))
            }
        }
    }
    ifelse(isTRUE(flgmlvl == FALSE) == TRUE, Lt <- abs(Lt), Lt <- unlist(Lt))
    vltz <- Lt
    if (missing(clu) == FALSE) {
        flgclu <- TRUE
        ifelse(is.numeric(clu) == TRUE && isTRUE(min(clu) < 0) == 
            TRUE, clu <- abs(clu), NA)
        if (isTRUE(flgmlvl == FALSE) == TRUE) {
            if (isTRUE(length(clu) == 1) == TRUE) {
                clu <- rep(clu, n)
                nclu <- 1L
            }
            else if (isTRUE(length(clu) > 1L) == TRUE) {
                if (is.vector(clu) == TRUE & is.factor(clu) == 
                  TRUE) {
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
                  for (i in seq_len(nlevels(factor(tmpclu)) - 
                    1L)) {
                    clu[which(clu == clu[which((clu %in% tmpclu) == 
                      TRUE)[(i - 0)]])] <- (i + 1L)
                  }
                  rm(i)
                  clu[which((clu %in% tmpclu) == TRUE)] <- nlevels(factor(tmpclu))
                  rm(tmpclu)
                }
                clu <- as.integer(as.vector(clu))
                if (missing(lclu) == FALSE && is.vector(lclu) == 
                  TRUE) {
                  nclu <- length(lclu)
                }
                else {
                  nclu <- nlevels(factor(clu))
                }
            }
        }
        else if (isTRUE(flgmlvl == TRUE) == TRUE) {
            if (any(unlist(lapply(clu, is.null))) == TRUE) {
                ifelse(unlist(lapply(clu, is.null))[1] == TRUE, 
                  clu[[1]] <- rep(1L, nrow(net)), NA)
                ifelse(unlist(lapply(clu, is.null))[2] == TRUE, 
                  clu[[2]] <- rep(1L, ncol(net)), NA)
            }
            else {
                NA
            }
            ifelse("cn2" %in% attr(mlv, "class"), clu <- clu[[1]], 
                NA)
            if (is.list(clu) == TRUE) {
                ifelse(is.factor(clu[[1]]) == TRUE, uact <- levels(clu[[1]]), 
                  uact <- unique(clu[[1]]))
                ifelse(is.factor(clu[[2]]) == TRUE, uevt <- levels(clu[[2]]), 
                  uevt <- unique(clu[[2]]))
                clutmp <- clu
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
        }
    }
    else {
        flgclu <- FALSE
        if (isTRUE(flgmlvl == TRUE) == TRUE && missing(vcol) == 
            FALSE && (is.vector(vcol) == TRUE || is.list(vcol) == 
            TRUE)) {
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
    else if (isTRUE(flgcx == FALSE) == TRUE && isTRUE(valued == 
        TRUE) == TRUE) {
        ifelse(isTRUE(max(cex) >= 21L) == TRUE, cex <- 20L, NA)
    }
    else {
        NA
    }
    if (missing(fsize) == TRUE) {
        ifelse(isTRUE(max(cex) < 2) == TRUE, fsize <- cex * 0.66, 
            fsize <- cex * 0.25)
    }
    else {
        fsize <- fsize/10L
    }
    ifelse(isTRUE(valued == FALSE) == TRUE && isTRUE(bwd > 1L) == 
        TRUE, bwd <- 1L, NA)
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
        if (isTRUE(flgmlvl == TRUE) == TRUE) {
            ifelse(missing(col) == TRUE, vcol <- c(rep("#C0C0C0", 
                length(mlv$lbs$dm)), rep("#408080", length(which(!(mlv$lbs$cdm %in% 
                mlv$lbs$dm))))), vcol <- col)
            ifelse(missing(vcol0) == TRUE, vcol0 <- rep(vcol, 
                2), NA)
        }
        else {
            ifelse(missing(col) == TRUE, vcol <- rep("#FFFFFF", 
                nclu), vcol <- col)
            ifelse(missing(vcol0) == TRUE, vcol0 <- rep("#808080", 
                nclu), NA)
        }
    }
    else {
        ifelse(isTRUE(length(vcol) == 1L) == TRUE, vcol <- rep(vcol, 
            n), NA)
    }
    if (isTRUE(flgmlvl == TRUE) == TRUE) {
        if (isTRUE(flgclu == FALSE) == TRUE) {
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
        }
    }
    else {
        if (isTRUE(length(vcol) != n) == TRUE & isTRUE(nclu == 
            1) == TRUE) {
            vcol <- rep(vcol[1], n)
        }
        else if ((isTRUE(length(vcol) == nclu) == TRUE) && (isTRUE(nclu < 
            length(vcol)) == TRUE || identical(vcol, clu) == 
            FALSE || missing(lclu) == FALSE)) {
            if (isTRUE(n > 1) == TRUE) {
                if ((is.factor(clu) == FALSE && isTRUE(min(clu) == 
                  0) == TRUE) && missing(lclu) == TRUE) {
                  lclu <- as.numeric(levels(as.factor(clu)))
                }
                else if (missing(lclu) == FALSE) {
                  NA
                }
                else {
                  lclu <- NULL
                }
                tmpvcol <- rep(0, n)
                if (is.null(lclu) == FALSE) {
                  if (isTRUE(min(lclu) == 0L) == TRUE) {
                    lclu[which(lclu == 0L)] <- max(lclu) + 1L
                    lclu <- sort(lclu)
                    clu[which(clu == 0L)] <- max(lclu)
                  }
                  for (i in seq_len(nclu)) {
                    tmpvcol[which(clu == lclu[i])] <- vcol[i]
                  }
                  rm(i)
                }
                else {
                  for (i in seq_len(nclu)) {
                    tmpvcol[which(clu == (levels(factor(clu))[i]))] <- vcol[i]
                  }
                  rm(i)
                }
                vcol <- tmpvcol
                rm(tmpvcol)
            }
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
        flgcrd <- TRUE
        if (missing(rm.isol) == FALSE && isTRUE(rm.isol == TRUE) == 
            TRUE) {
            ifelse(all(lbs %in% rownames(coord)) == FALSE, stop(paste("Rownames in ", 
                sQuote("coord"), " do not match network labels in ", 
                sQuote("net"), ".", sep = ""), call. = FALSE), 
                crd <- coord[which(rownames(coord) %in% lbs), 
                  ])
        }
        else {
            if (isTRUE(nrow(coord) == n) == FALSE) 
                stop(paste("Length of ", sQuote("coord"), " does not match number of vertices in ", 
                  sQuote("net"), ".", sep = ""))
            crd <- coord
        }
    }
    else if (missing(coord) == TRUE) {
        flgcrd <- FALSE
        switch(match.arg(layout), force = {
            ifelse(isTRUE(n < 2L) == TRUE, crd <- data.frame(X = sin(2L * 
                pi * ((0:(n - 1L))/n)), Y = cos(2L * pi * ((0:(n - 
                1L))/n))), crd <- frcd(netd, seed = seed, maxiter = maxiter))
            ifelse(missing(jitter) == FALSE, crd <- apply(crd, 
                2, function(z) {
                  jitter(z, amount = (min(10, jitter))/100L)
                }), NA)
        }, stress = {
            ifelse(isTRUE(n < 2L) == TRUE, crd <- data.frame(X = sin(2L * 
                pi * ((0:(n - 1L))/n)), Y = cos(2L * pi * ((0:(n - 
                1L))/n))), crd <- stsm(netd, seed = seed, maxiter = maxiter, 
                ...))
            ifelse(missing(jitter) == FALSE, crd <- apply(crd, 
                2, function(z) {
                  jitter(z, amount = (min(10, jitter))/100L)
                }), NA)
        }, circ = {
            crd <- data.frame(X = sin(2L * pi * ((0:(n - 1L))/n)), 
                Y = cos(2L * pi * ((0:(n - 1L))/n)))
        }, conc = {
            if (missing(lbs) == FALSE && is.null(dimnames(netd)[[1]]) == 
                TRUE) {
                if (identical(lbs, make.unique(lbs)) == FALSE) {
                  warning(paste("Make unique ", sQuote("lbs"), 
                    " because of duplicated entries.", sep = ""), 
                    call. = FALSE)
                  dimnames(netd)[[1]] <- dimnames(netd)[[2]] <- make.unique(lbs)
                } else {
                  dimnames(netd)[[1]] <- dimnames(netd)[[2]] <- lbs
                }
            }
            if (missing(nr) == TRUE) {
                nr <- 1
            } else {
                if (isTRUE(dim(netd)[1] < length(nr)) == TRUE) {
                  warning(paste("Number of radii in ", sQuote("nr"), 
                    " larger than number of vertices in ", sQuote("net"), 
                    ".", sep = ""), call. = FALSE)
                  nr <- nr[1:dim(netd)[1]]
                } else {
                  NA
                }
            }
            crd <- conc(netd, nr = nr, ...)
        }, rand = {
            set.seed(seed)
            crd <- data.frame(X = round(stats::runif(n) * 1L, 
                5), Y = round(stats::runif(n) * 1L, 5))
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
    else {
        if (isTRUE(directed == TRUE) == TRUE) {
            if (isTRUE(vedist > 5L) == TRUE) {
                vedist <- 5L
            }
            else if (isTRUE(vedist < (-5L)) == TRUE) {
                vedist <- -5L
            }
        }
        else {
            if (isTRUE(vedist > 10L) == TRUE) {
                vedist <- 10L
            }
            else if (isTRUE(vedist < (-10L)) == TRUE) {
                vedist <- -10L
            }
        }
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
    ifelse(missing(flip) == FALSE && isTRUE(flip == FALSE) == 
        TRUE, flip <- FALSE, flip <- TRUE)
    if (match.arg(layout) %in% c("circ", "conc")) {
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
    if ((isTRUE(flip == TRUE) == TRUE && isTRUE(loops == TRUE) == 
        TRUE) & (match.arg(layout) %in% c("stress", "force", 
        "rand", "conc") || missing(coord) == FALSE)) {
        crdf <- cbind(crd[, 1] - mean(crd[, 1]), crd[, 2] - mean(crd[, 
            2]))
        angl <- vector(length = nrow(crd))
        for (i in seq_len(nrow(crd))) {
            ifelse((atan2((crdf[i, 2]), (crdf[i, 1])) * (180L/pi)) >= 
                0, angl[i] <- atan2((crdf[i, 2]), (crdf[i, 1])) * 
                (180L/pi), angl[i] <- ((atan2((crdf[i, 2]), (crdf[i, 
                1])) * (180L/pi))%%180L) + 180L)
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
    if (isTRUE(flgmlvl == TRUE) == TRUE) {
        if (missing(pos) == TRUE) {
            if (isTRUE("cn2" %in% attr(mlv, "class")) == TRUE || 
                (all(mlv$lbs$dm %in% mlv$lbs$cdm) == TRUE && 
                  all(mlv$lbs$cdm %in% mlv$lbs$dm) == TRUE)) {
                pos <- c(0, 0)
            }
            else {
                pos <- c(2, 4)
            }
        }
        else if (missing(pos) == FALSE) {
            if (isTRUE("cn2" %in% attr(mlv, "class")) == TRUE || 
                (all(mlv$lbs$dm %in% mlv$lbs$cdm) == TRUE && 
                  all(mlv$lbs$cdm %in% mlv$lbs$dm) == TRUE)) {
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
    }
    else {
        if (missing(pos) == TRUE) {
            pos <- 4
        }
        else {
            if (isTRUE(pos < 0L) == TRUE | isTRUE(pos > 4L) == 
                TRUE) {
                stop(paste("Values in ", sQuote("pos"), " must be between 0 and 4.", 
                  sep = ""), call. = FALSE)
            }
            else {
                invisible(NA)
            }
        }
    }
    ifelse(missing(mirrorX) == FALSE && isTRUE(mirrorX == TRUE) == 
        TRUE || missing(mirrorV) == FALSE && isTRUE(mirrorV == 
        TRUE) == TRUE, crd[, 1] <- crd[, 1] * cos(pi) - crd[, 
        2] * sin(pi), mirrorX <- FALSE)
    if (isTRUE(n > 1) == TRUE) {
        ifelse(missing(mirrorY) == FALSE && isTRUE(mirrorY == 
            TRUE) == TRUE || missing(mirrorH) == FALSE && isTRUE(mirrorH == 
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
    }
    else {
        mirrorY <- FALSE
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
        NA
    }
    if (missing(mai) == TRUE) {
        ifelse(isTRUE(loops == TRUE) == TRUE, mai <- c(0, 0, 
            0.2, 0), mai <- c(0, 0, 0, 0))
    }
    else {
        NA
    }
    graphics::par(mar = mar)
    ifelse(is.null(main) == TRUE, graphics::par(mai = mai), graphics::par(mai = mai + 
        c(0, 0, cex.main/10L, 0)))
    obg <- graphics::par()$bg
    graphics::par(bg = grDevices::adjustcolor(bg, alpha = alpha[3]))
    if (isTRUE(loops == TRUE) == TRUE && isTRUE(flip == TRUE) == 
        TRUE && missing(new) == TRUE) {
        ylim <- c(min(nds[, 2]) - (max(cex)/50L), max(nds[, 2]) + 
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
    ifelse(missing(new) == FALSE && isTRUE(new == TRUE) == TRUE, 
        suppressWarnings(par(new = TRUE)), NA)
    suppressWarnings(graphics::plot(nds, type = "n", axes = FALSE, 
        xlab = "", ylab = "", ylim = ylim, xlim = xlim, asp = asp, 
        main = main, cex.main = cex.main, col.main = col.main, 
        ...))
    if (isTRUE(showLbs == TRUE) == TRUE || isTRUE(values == TRUE) == 
        TRUE) {
        ifelse(missing(ffamily) == FALSE && isTRUE(ffamily %in% 
            names(grDevices::postscriptFonts())) == TRUE, graphics::par(family = ffamily), 
            NA)
    }
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
        tlbs <- tlbs[which(tlbs != "loop")]
    }
    else {
        NA
    }
    if (isTRUE(flgcn2 == TRUE) == TRUE) {
        tlbs <- rev(tlbs)
        bds <- rev(bds)
        bd <- rev(bd)
    }
    else {
        NA
    }
    ifelse(isTRUE(flgpf == TRUE) == TRUE, fnnetdrp <- maxw, fnnetdrp <- (norm(as.matrix(netdrp), 
        type = "F")))
    if (isTRUE(swp == TRUE) == TRUE) {
        Lt <- Lt[rev(seq_len(length(Lt)))]
        lwd <- lwd[length(lwd):1]
        alfa <- alfa[rev(seq_len(length(alfa)))]
    }
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
    if ((isTRUE(collRecip == TRUE) == TRUE && isTRUE(undRecip == 
        TRUE) == FALSE) && isTRUE(valued == TRUE) == FALSE && 
        isTRUE(c("recp") %in% attr(bds, "names")) == TRUE) {
        trcp <- multiplex::trnf(rcp, tolist = TRUE)
    }
    else {
        NA
    }
    if (isTRUE(values == TRUE) == TRUE) {
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
                else {
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
                    if (isTRUE(Lt == Ltc) == FALSE) {
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
                if (isTRUE(flgcx == TRUE) == TRUE) {
                  if ((isTRUE(directed == TRUE) == TRUE && isTRUE(tlbs[k] %in% 
                    c("asym", "tent", "txch", "mixd")) == TRUE) | 
                    (isTRUE(collRecip == FALSE) == TRUE && isTRUE(tlbs[k] %in% 
                      c("asym", "tent", "txch", "mixd")) == TRUE) | 
                    (isTRUE(directed == TRUE) == FALSE && isTRUE(tlbs[k] %in% 
                      c("txch")) == TRUE)) {
                    if (isTRUE(multiplex::dhc(bds[[k]][1])[1] > 
                      multiplex::dhc(bds[[k]][1])[2]) == TRUE) {
                      Prs <- vector(length = length(prs))
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
                if (isTRUE(valued == TRUE) == TRUE) {
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
                    lw <- (Lw/fnnetdrp) * (10L * 5L)
                  }
                  else {
                    lw <- Lw
                  }
                }
                else {
                  lw <- rep(lwd[1], rbds)
                }
                ifelse((isTRUE(swp2 == TRUE) == TRUE && isTRUE(tlbs[k] %in% 
                  c("recp")) == TRUE) && isTRUE(valued == FALSE) == 
                  TRUE, bds[[k]] <- multiplex::swp(bds[[k]]), 
                  NA)
                if (isTRUE(collRecip == TRUE) == TRUE && isTRUE(tlbs[k] %in% 
                  c("recp")) == TRUE) {
                  bw <- 0L
                  if (isTRUE(undRecip == TRUE) == TRUE) {
                    hd <- 0
                  }
                  else {
                    hd <- hds
                  }
                }
                else if (isTRUE(collRecip == TRUE) == FALSE && 
                  isTRUE(tlbs[k] %in% c("recp")) == TRUE) {
                  ifelse(isTRUE(undRecip == TRUE) == TRUE, hd <- 0L, 
                    hd <- hds)
                  bw <- bwd
                }
                else {
                  bw <- bwd
                  hd <- hds
                }
                if (isTRUE((collRecip == TRUE) == TRUE && isTRUE(undRecip == 
                  TRUE) == FALSE) && isTRUE(tlbs[k] %in% c("recp")) == 
                  FALSE && isTRUE(valued == TRUE) == FALSE && 
                  isTRUE(c("recp") %in% attr(bds, "names")) == 
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
                if (isTRUE(z == 1L) == TRUE) {
                  mbnd(pars, rbds, bds[[k]], vlt, cx * mscl, 
                    lw, vecol, directed, bw, alfa, fds, flgcx, 
                    valued, flgcr, hd, n)
                }
                else {
                  ifelse(isTRUE(length(lty) == 1L) == TRUE, mbnd(pars, 
                    rbds, bds[[k]], vlt1, cx * mscl, lw, vecol[vltc], 
                    directed, bw, alfa, fds, flgcx, valued, flgcr, 
                    hd, n), mbnd(pars, rbds, bds[[k]], vlt, cx * 
                    mscl, lw, vecol[vltc], directed, bw, alfa, 
                    fds, flgcx, valued, flgcr, hd, n))
                }
                if (isTRUE(values == TRUE) == TRUE) {
                  if (isTRUE(z == 1L) == TRUE) {
                    ifelse(isTRUE(isSymmetric(net) == TRUE) == 
                      TRUE, fsym <- 2L, fsym <- 1L)
                    ifelse(isTRUE(directed == FALSE) == TRUE, 
                      mdp[nrow(mdp) + 1, ] <- c((pars[1, 1] + 
                        pars[2, 1])/2L, (pars[1, 2] + pars[2, 
                        2])/2L, as.numeric(net[as.numeric(multiplex::dhc(bds[[k]])[1]), 
                        as.numeric(multiplex::dhc(bds[[k]])[2])]/fsym + 
                        net[as.numeric(multiplex::dhc(bds[[k]])[2]), 
                          as.numeric(multiplex::dhc(bds[[k]])[1])]/fsym)), 
                      mdp[nrow(mdp) + 1, ] <- c((pars[1, 1] + 
                        pars[2, 1])/2L, (pars[1, 2] + pars[2, 
                        2])/2L, as.numeric(net[as.numeric(multiplex::dhc(bds[[k]])[1]), 
                        as.numeric(multiplex::dhc(bds[[k]])[2])]) + 
                        net[as.numeric(multiplex::dhc(bds[[k]])[2]), 
                          as.numeric(multiplex::dhc(bds[[k]])[1])]))
                  }
                  else {
                    mpnet <- multiplex::mnplx(net, directed = TRUE, 
                      valued = TRUE)
                    ifelse(isTRUE(directed == FALSE) == TRUE, 
                      mdp[nrow(mdp) + 1, ] <- c((pars[1, 1] + 
                        pars[2, 1])/2L, (pars[1, 2] + pars[2, 
                        2])/2L, as.numeric(mpnet[as.numeric(multiplex::dhc(bds[[k]])[1]), 
                        as.numeric(multiplex::dhc(bds[[k]])[2])] + 
                        mpnet[as.numeric(multiplex::dhc(bds[[k]])[2]), 
                          as.numeric(multiplex::dhc(bds[[k]])[1])])), 
                      mdp[nrow(mdp) + 1, ] <- c((pars[1, 1] + 
                        pars[2, 1])/2L, (pars[1, 2] + pars[2, 
                        2])/2L, multiplex::mnplx(netdrp, directed = FALSE, 
                        valued = TRUE)[as.numeric(multiplex::dhc(bds[[k]])[1]), 
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
        ifelse(isTRUE(loops == TRUE) == TRUE && isTRUE(length(lty) == 
            1) == TRUE, vlt <- lty, NA)
    }
    if (isTRUE(loops == TRUE) == TRUE) {
        netdrpl <- netdrp
        if (isTRUE(swp == TRUE) == TRUE) {
            bdlp <- bd$loop[rev(seq_len(length(bd$loop)))]
            ifelse(isTRUE(flip == TRUE) == FALSE && is.na(dim(netdrp)[3]) == 
                TRUE, NA, netdrpl <- netdrp[, , rev(seq_len(dim(netdrp)[3]))])
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
                    dcx <- dcx - (cex[lp[i]] * 0.011)
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
                  if (isTRUE(valued == TRUE) == TRUE) {
                    lwl <- diag(netdrpl)[lp[i]]
                  }
                  else {
                    if (isTRUE(Ltl == 1) == TRUE) {
                      lwl <- lwd - 0.5
                    }
                    else {
                      lwl <- lwd
                    }
                  }
                  hc(ndss[lp[i], 1], ndss[lp[i], 2] + (dcx), 
                    rds, col = grDevices::adjustcolor(vecol, 
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
                  ifelse(isTRUE(z > 2) == TRUE, bwd2 <- 2L, bwd2 <- 4L)
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
                          rds <- rds * 0.85
                          if (isTRUE(cex[lp[i]] < 3L) == TRUE) {
                            dcx <- dcx - 0.5 * (k/1000L) - (cex[lp[i]] * 
                              0.012)
                          }
                          else {
                            dcx <- dcx - 2L * (k/1000L) - (cex[lp[i]] * 
                              0.012)
                          }
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
                        rds <- rds * 0.85
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
                      lwl <- netdrpl[i, i, k]
                    }
                    else {
                      if (isTRUE(Ltl[k] == 1) == TRUE) {
                        lwl <- lwd[k] - k
                      }
                      else {
                        lwl <- lwd[k]
                      }
                    }
                    hc(ndss[lp[i], 1], ndss[lp[i], 2] + (dcx), 
                      rds, col = grDevices::adjustcolor(vecol[k], 
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
    if (isTRUE(values == TRUE) == TRUE) {
        mdps <- mdp[, 1:2]
        mdps[, 1] <- mdps[, 1] * scl[1]
        mdps[, 2] <- mdps[, 2] * scl[2]
        ifelse(missing(fsize2) == TRUE, fsize2 <- fsize, fsize2 <- fsize2/10L)
        ifelse(missing(cex2) == TRUE, cex2 <- fsize2, NA)
        ifelse(missing(fcol2) == TRUE, fcol2 <- "#000000", NA)
        ifelse(missing(bg2) == TRUE, bg2 <- "transparent", NA)
        if (isTRUE(flgbg == TRUE) == TRUE) {
            graphics::points(mdps[, 1], mdps[, 2], pch = 16, 
                cex = fsize2 * 2, col = bg2, bg = bg)
        }
        else if (isTRUE(flgbg == FALSE) == TRUE) {
            graphics::points(mdps[, 1], mdps[, 2], pch = 16, 
                cex = fsize2 * 2, col = grDevices::adjustcolor(bg2, 
                  alpha = alpha[4]), bg = obg)
        }
        else {
            NA
        }
        if (isTRUE(directed == FALSE) == TRUE) {
            for (k in seq_along(tlbs)) {
                graphics::text(mdps[k, ], labels = mdp[k, 3], 
                  cex = fsize2, col = fcol2)
            }
            rm(k)
        }
        else {
            for (k in seq_along(tlbs)) {
                graphics::text(mdps[k, ], labels = mdp[k, 3], 
                  cex = fsize2, col = fcol2)
            }
            rm(k)
        }
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
        if (missing(sel) == TRUE) {
            lbs[which(is.na(lbs))] <- ""
        }
        else {
            ifelse(is.data.frame(sel) == TRUE, sel <- as.vector(unlist(sel, 
                use.names = FALSE)), NA)
            lbb <- lbs
            lbs <- rep("", length(lbs))
            vec <- vector()
            if (all(sel %in% lbb) == FALSE) {
                warning(paste("Values of ", sQuote("sel"), " are not found in ", 
                  sQuote("net"), ".", sep = ""), call. = FALSE)
                sel0 <- sel[which(sel %in% lbb)]
                for (i in seq_len(length(sel0))) {
                  vec <- append(vec, which(lbb %in% sel0[i]))
                }
                rm(i)
                suppressWarnings(lbs[vec] <- sel0)
            }
            else {
                for (i in seq_len(length(sel))) {
                  vec <- append(vec, which(lbb %in% sel[i]))
                }
                rm(i)
                suppressWarnings(lbs[vec] <- sel)
            }
        }
        if (isTRUE(flgmlvl == TRUE) == TRUE) {
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
                  else if (isTRUE(fstyle == "bolditalic") == 
                    TRUE) {
                    graphics::text(ndss[seq_len(nn), ], labels = as.expression(lapply(lbs[seq_len(nn)], 
                      function(w) bquote(bolditalic(.(w))))), 
                      cex = fsize, adj = 0.5, col = fcol[1])
                  }
                }
            }
            else {
                if (missing(fstyle) == TRUE || (missing(fstyle) == 
                  FALSE && isTRUE(fstyle %in% c("italic", "bold", 
                  "bolditalic") == FALSE))) {
                  graphics::text(ndss[seq_len(nn), ], lbs[seq_len(nn)], 
                    cex = fsize, pos = pos[1], col = fcol[1], 
                    offset = (cex/4L), adj = c(0.5, 1))
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
                  else if (isTRUE(fstyle == "bolditalic") == 
                    TRUE) {
                    graphics::text(ndss[seq_len(nn), ], as.expression(lapply(lbs[seq_len(nn)], 
                      function(w) bquote(bolditalic(.(w))))), 
                      cex = fsize, pos = pos[1], col = fcol[1], 
                      offset = (cex/4L), adj = c(0.5, 1))
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
                  else if (isTRUE(fstyle == "bolditalic") == 
                    TRUE) {
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
                  else if (isTRUE(fstyle == "bolditalic") == 
                    TRUE) {
                    graphics::text(ndss[(nn + 1L):dim(mlv$mlnet)[1], 
                      ], as.expression(lapply(lbs[(nn + 1L):dim(mlv$mlnet)[1]], 
                      function(w) bquote(bolditalic(.(w))))), 
                      cex = fsize, pos = pos[2], col = fcol[2], 
                      offset = (cex/4L), adj = c(0.5, 1))
                  }
                }
            }
        }
        else {
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
                  else if (isTRUE(fstyle == "bolditalic") == 
                    TRUE) {
                    graphics::text(ndss, as.expression(lapply(lbs, 
                      function(z) bquote(bolditalic(.(z))))), 
                      cex = fsize, pos = pos, col = fcol[1], 
                      offset = (cex/4L), adj = c(0.5, 1))
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
                      ifelse(isTRUE(flgad == TRUE) == TRUE, neta <- suppressMessages(multiplex::zbind(netd, 
                        att, force = TRUE)), neta <- multiplex::zbind(netd, 
                        att))
                      dimnames(neta)[[3]] <- c(dimnames(net)[[3]], 
                        dimnames(att)[[3]])
                      clss <- multiplex::expos(multiplex::rel.sys(neta, 
                        att = (z + 1L):dim(neta)[3]), classes = TRUE)$Classes
                      attr(clss, "names")[which(attr(clss, "names") == 
                        "ALL")] <- multiplex::jnt(dimnames(att)[[3]], 
                        sep = "")
                      ifelse(isTRUE(attr(clss, "names")[1] == 
                        "NONE") == TRUE, nnone <- 2, nnone <- 1)
                      for (i in nnone:length(clss)) {
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
    options(warn = 0)
    graphics::par(mar = omr)
    graphics::par(bg = obg)
    graphics::par(lend = 0)
    graphics::par(mai = omi)
}

