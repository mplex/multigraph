circn <-
function (net, clu, inv, flip, irot, jitter, ...) 
{
    n <- dim(net)[1]
    ifelse(missing(clu) == TRUE, clu <- rep(1, n), NA)
    if (isTRUE(is.null(dimnames(net)[1]) == TRUE | is.null(dimnames(net)[1][[1]])) == 
        FALSE) {
        lbs <- dimnames(net)[[1]]
    }
    else {
        lbs <- 1:n
    }
    if (isTRUE(length(clu) == 1) == TRUE && is.numeric(clu) == 
        TRUE) {
        if (isTRUE(clu >= n) == TRUE | isTRUE(clu <= 0) == TRUE) 
            stop("Value of 'clu' must be greater than zero and lower than network order.")
        r <- clu
        clu <- vector()
        for (i in 1:r) {
            clu <- append(clu, rep(i, ceiling(n/r)))
        }
        rm(i)
    }
    else if (is.factor(clu) == TRUE) {
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
            clu[which(clu == clu[which((clu %in% tmpclu) == TRUE)[(i - 
                0)]])] <- (i + 1)
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
    ifelse(isTRUE(0 %in% (clu)) == TRUE, clu <- clu + 1L, NA)
    ifelse(isTRUE(length(clu) > n) == TRUE, clu <- clu[1:n], 
        NA)
    r <- nlevels(factor(clu))
    if (missing(inv) == FALSE && isTRUE(inv == TRUE) == TRUE) {
        oclu <- clu
        k <- 1L
        for (i in nlevels(factor(oclu)):1) {
            clu[which(levels(factor(oclu))[i] == oclu)] <- as.numeric(levels(factor(oclu))[k])
            k <- k + 1L
        }
        rm(i)
    }
    else {
        inv <- FALSE
    }
    nlst <- list()
    length(nlst) <- r
    for (i in 1:r) {
        nlst[[i]] <- lbs[which(clu == i)]
    }
    rm(i)
    rad <- 1
    nds <- data.frame(matrix(ncol = 2, nrow = 0))
    for (i in 1:length(nlst)) {
        nl <- length(nlst[[i]])
        if (isTRUE(i == 1) == TRUE && isTRUE(nl == 1) == TRUE) {
            x <- data.frame(X = 0, Y = 0)
            rad <- rad - 1L
        }
        else {
            x <- data.frame(X = cos(2L * pi * ((0:((nl) - 1L))/nl)) * 
                rad, Y = sin(2L * pi * ((0:(nl - 1L))/nl)) * 
                rad)
        }
        if (missing(flip) == FALSE && isTRUE(flip == TRUE) == 
            TRUE) {
            ifelse(isTRUE((i%%2L) == 0) == TRUE, x[, 1] <- x[, 
                1] * cos(pi) - x[, 2] * sin(pi), NA)
        }
        else {
            NA
        }
        if (isTRUE((i%%2L) == 0) == TRUE) {
            ifelse(missing(irot) == FALSE && is.numeric(irot) == 
                TRUE, x[, 2:1] <- xyrt(x[, 2:1], (irot * -1L)), 
                x[, 2:1] <- xyrt(x[, 2:1], -0L))
        }
        else {
            NA
        }
        nds <- rbind(nds, x[, 2:1])
        rad <- rad + 1L
    }
    rm(i)
    if (isTRUE(r > 1) == TRUE) {
        rownames(nds) <- dimnames(net)[[1]]
        nnds <- data.frame(matrix(ncol = ncol(nds), nrow = nrow(nds)))
        nnds[which(rownames(nds) %in% nlst[[1]]), ] <- nds[1:length(nlst[[1]]), 
            ]
        yp <- tabulate(clu)[1]
        for (i in 1:(r - 1)) {
            k <- 1L
            for (j in which(rownames(nds) %in% nlst[[i + 1L]])) {
                nnds[j, ] <- nds[(yp + 1L):(yp + tabulate(clu)[i + 
                  1L]), ][k, ]
                k <- k + 1L
            }
            rm(j)
            yp <- yp + tabulate(clu)[i + 1L]
        }
        rm(i)
        ifelse(missing(jitter) == FALSE, nnds <- jitter(as.matrix(nnds), 
            amount = jitter), NA)
        nnds
    }
    else {
        ifelse(missing(jitter) == FALSE, nds <- jitter(as.matrix(nds), 
            amount = jitter), NA)
        nds
    }
}
