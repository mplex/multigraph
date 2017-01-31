stsm <-
function (net, seed = seed, maxiter = 40, drp, jitter, method, 
    ...) 
{
    n <- dim(net)[1]
    ifelse(missing(drp) == FALSE && is.numeric(drp) == TRUE, 
        NA, drp <- 0)
    netd <- multiplex::mnplx(replace(net, net <= drp, 0), directed = FALSE, 
        dichot = TRUE)
    ifelse(missing(method) == TRUE, delta <- stats::dist(netd, 
        upper = TRUE, diag = TRUE, method = "binary"), delta <- stats::dist(netd, 
        upper = TRUE, diag = TRUE, method = method))
    ifelse(isTRUE(sum(delta) == 0) == TRUE, return(popl(dim(net)[1], 
        seed = seed)), NA)
    mwd <- delta^(-2)
    mwd[which(mwd == Inf)] <- 0
    lpmwd <- as.matrix(mwd)
    for (i in 1:n) {
        for (j in 1:n) {
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
    X0 <- cbind(vec[1:n], vec[(n + 1):length(vec)])
    Xs <- as.matrix(X0)
    newstss <- sts(X0, delta = delta, mwd = mwd)
    for (iter in 1:maxiter) {
        if (any(is.nan(lz(X0, delta, mwd))) == FALSE) {
            X <- pilpmwd %*% (lz(X0, delta, mwd) %*% X0)
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
    if (isTRUE(length(multiplex::comps(netd)$isol) > 2) == TRUE) {
        isol <- which(duplicated(round(Xs, 3)) | duplicated(round(Xs, 
            3), fromLast = TRUE))
        tmpi <- popl(length(isol), seed = seed) * (1/length(isol) * 
            length(isol))
        locx <- ((tmpi[, 1]/3L) - (min(Xs[, 1])) - 0)
        ifelse(isTRUE((max(Xs[, 1]) - min(Xs[, 1]))/(max(Xs[, 
            2]) - min(Xs[, 2])) > 0) == TRUE, locy <- ((min(Xs[, 
            2])) - (tmpi[, 2]/3L) - 0), locy <- ((max(Xs[, 2])) + 
            (tmpi[, 2]/3L) + 0))
        Xs.chull <- grDevices::chull(Xs)
        Xs.chull <- Xs[Xs.chull, ]
        ifelse(isTRUE(length(which(Xs.chull[, 1] < mean(Xs.chull[, 
            1]))) > length(which(Xs.chull[, 1] > mean(Xs.chull[, 
            1])))) == TRUE, locx <- locx + (2L/10L), locx <- locx + 
            ((2L/10L) * -1))
        ifelse(isTRUE(length(which(Xs.chull[, 2] < mean(Xs.chull[, 
            2]))) > length(which(Xs.chull[, 2] > mean(Xs.chull[, 
            2])))) == TRUE, locy <- locy - (2L/10L), locy <- locy - 
            ((2L/10L) * -1))
        Xs[isol, ] <- (cbind(locx, locy))
    }
    else {
        NA
    }
    if (missing(jitter) == FALSE && isTRUE(jitter) == TRUE) {
        jitter(Xs, amount = (n/100))
    }
    else {
        return(as.data.frame(Xs))
    }
}
