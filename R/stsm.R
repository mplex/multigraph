stsm <-
function (net, seed = seed, maxiter = 40, drp, jitter, method, 
    ...) 
{
    n <- dim(net)[1]
    ifelse(missing(drp) == FALSE && is.numeric(drp) == TRUE, 
        NA, drp <- 0)
    netd <- multiplex::mnplx(replace(net, net <= drp, 0), directed = FALSE, 
        dichot = TRUE)
    ifelse(is.null(rownames(netd)) == TRUE, lbs <- seq_len(n), 
        lbs <- rownames(netd))
    ifelse(missing(method) == TRUE, delta <- stats::dist(netd, 
        upper = TRUE, diag = TRUE, method = "binary"), delta <- stats::dist(netd, 
        upper = TRUE, diag = TRUE, method = method))
    ifelse(isTRUE(sum(delta) == 0) == TRUE, return(popl(dim(net)[1], 
        seed = seed)), NA)
    mwd <- delta^(-2)
    mwd[which(mwd == Inf)] <- 0
    lpmwd <- as.matrix(mwd)
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            ifelse(isTRUE(as.matrix(mwd)[i, j] != 0) == TRUE, 
                lpmwd[i, j] <- -1 * as.matrix(mwd)[i, j], NA)
        }
        rm(j)
    }
    rm(i)
    diag(lpmwd) <- apply(as.matrix(mwd), 1, sum)
    s <- svd(lpmwd)
    p <- (s$d > max(.Machine$double.eps^(2/3) * s$d[1], 0))
    if (all(p)) {
        pilpmwd <- s$v %*% (1/s$d * t(s$u))
    }
    else if (any(p)) {
        pilpmwd <- s$v[, p, drop = FALSE] %*% (1/s$d[p] * t(s$u[, 
            p, drop = FALSE]))
    }
    else {
        pilpmwd <- matrix(0, nrow = ncol(lpmwd), ncol = nrow(lpmwd))
    }
    set.seed(seed)
    vec <- stats::rnorm(n * 2L)
    X0 <- cbind(vec[seq_len(n)], vec[(n + 1):length(vec)])
    Xs <- as.matrix(X0)
    newstss <- sts(X0, delta = delta, mwd = mwd)
    for (iter in seq_len(maxiter)) {
        lzx0 <- lz(X0, delta, mwd)
        if (any(is.nan(lzx0)) == FALSE) {
            X <- pilpmwd %*% (lzx0 %*% X0)
            X[which(X == Inf)] <- 0
            oldstss <- newstss
            newstss <- sts(X, delta = delta, mwd = mwd)
            Xs <- X
            abstols <- abstolx <- reltols <- sqrt(.Machine$double.eps)
            ifelse(isTRUE(abs(newstss - oldstss) < (reltols * 
                newstss)) == TRUE, break, NA)
            ifelse(isTRUE(abs(newstss - oldstss) < abstols) == 
                TRUE, break, NA)
            ifelse(isTRUE(norm(X - X0, type = "F") < abstolx) == 
                TRUE, break, NA)
            X0 <- X
        }
    }
    rm(iter)
    cmps <- multiplex::comps(netd)
    nds <- Xs
    rownames(nds) <- lbs
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
            2) * length(cmps$isol)
        if (is.null(cmps$com) == FALSE) {
            locx <- ((tmpi[, 1]/3L) - (min(ndst[, 1])) - 0)
            ifelse(isTRUE(rat > 0) == TRUE, locy <- ((min(ndst[, 
                2])) - (tmpi[, 2]/3L) - 0), locy <- ((max(ndst[, 
                2])) + (tmpi[, 2]/3L) + 0))
            ndst.chull <- grDevices::chull(ndst)
            ndst.chull <- ndst[ndst.chull, ]
            ifelse(isTRUE(length(which(ndst.chull[, 1] < mean(ndst.chull[, 
                1]))) > length(which(ndst.chull[, 1] > mean(ndst.chull[, 
                1])))) == TRUE, locx <- locx + (1/n), locx <- locx + 
                ((1/n) * -1))
            ifelse(isTRUE(length(which(ndst.chull[, 2] < mean(ndst.chull[, 
                2]))) > length(which(ndst.chull[, 2] > mean(ndst.chull[, 
                2])))) == TRUE, locy <- locy - (1/n), locy <- locy - 
                ((1/n) * -1))
        }
        else {
            locx <- (tmpi[, 1])
            locy <- (tmpi[, 2])
        }
        nds[which(lbs %in% cmps$isol), ] <- (cbind(locx, locy))
    }
    else if (isTRUE(length(cmps$isol) == 1L) == TRUE) {
        ndst <- nds[which(nds[, 1] != 0), ]
        locx <- max(ndst[, 1]) + (1/n)
        ifelse(isTRUE(rat < 0) == TRUE, locy <- max(ndst[, 2]) + 
            (1/n), locy <- min(ndst[, 2]) - (1/n))
        nds[which(lbs %in% cmps$isol), ] <- (cbind(locx, locy))
    }
    nds[, 1] <- (nds[, 1] - min(nds[, 1]))/(max(nds[, 1]) - min(nds[, 
        1])) + (1/n)
    ifelse(isTRUE(rat > 0) == TRUE, nds[, 2] <- ((nds[, 2] - 
        min(nds[, 2]))/(max(nds[, 2]) - min(nds[, 2]))) * (1L/rat), 
        nds[, 2] <- ((nds[, 2] - min(nds[, 2]))/(max(nds[, 2]) - 
            min(nds[, 2]))) * (rat) + (1/n))
    nds[, 2] <- nds[, 2] * -1
    Xs <- as.data.frame(nds)
    if (missing(jitter) == FALSE && isTRUE(jitter) == TRUE) {
        jitter(Xs, amount = (n/100))
    }
    else {
        return(as.data.frame(Xs))
    }
}
