## (GPL-3)  Antonio Rivero Ostoic, multiplex@post.com
conc <-
function (net, nr, irot, inv, iflip, mirror = c("N", "X", "Y", 
    "D", "L"), ...) 
{
    n <- dim(net)[1]
    if (isTRUE((is.null(dimnames(net)[1]) == TRUE | is.null(dimnames(net)[1][[1]])) == 
        FALSE) && any(duplicated(dimnames(net)[[1]])) == FALSE) {
        lbs <- dimnames(net)[[1]]
    }
    else {
        lbs <- seq_len(n)
    }
    if (missing(nr) == TRUE) {
        clu <- rep(1, n)
    }
    else if (isTRUE(length(nr) == 1L) == TRUE && is.numeric(nr) == 
        TRUE) {
        if (isTRUE(nr >= n) == TRUE | isTRUE(nr <= 0) == TRUE) {
            message(paste("Number of radii in ", sQuote(deparse(substitute(nr))), 
                " must be positive and lower than network order.", 
                sep = ""))
            nr <- 1L
        }
        else if (isTRUE(is.integer(nr) == FALSE) == TRUE) {
            nr <- round(nr)
        }
        else {
            NA
        }
        clu <- vector()
        for (i in seq_len(nr)) {
            clu <- append(clu, rep(i, ceiling(n/nr)))
        }
        rm(i)
        clu <- clu[seq_len(n)]
    }
    else if ((isTRUE(length(nr) > 1L) == TRUE && is.vector(nr) == 
        TRUE) || is.factor(nr) == TRUE) {
        tmpnr <- nr[seq_len(n)]
        if (any(is.na(tmpnr)) == TRUE) 
            stop(paste("Missing data found in ", sQuote(deparse(substitute(nr))), 
                ".", sep = ""))
        for (i in seq_len(nlevels(factor(nr)))) {
            levels(nr) <- c(levels(nr), i)
            nr[which(levels(factor(tmpnr))[i] == nr)] <- i
        }
        rm(i)
        clu <- as.integer(as.vector(nr))
        rm(tmpnr)
    }
    else if (is.character(nr) == TRUE) {
        tmpnr <- nr[seq_len(n)]
        if (any(is.na(tmpnr)) == TRUE) 
            stop(paste("Missing data found in ", sQuote(deparse(substitute(nr))), 
                ".", sep = ""))
        nr[which(nr == nr[1])] <- 1
        for (i in seq_len(nlevels(factor(tmpnr)) - 1L)) {
            nr[which(nr == nr[which((nr %in% tmpnr) == TRUE)[(i - 
                0)]])] <- (i + 1L)
        }
        rm(i)
        nr[which((nr %in% tmpnr) == TRUE)] <- nlevels(factor(tmpnr))
        clu <- as.integer(as.vector(nr))
        rm(tmpnr)
    }
    else {
        NA
    }
    r <- nlevels(factor(clu))
    if (missing(irot) == FALSE && is.numeric(irot) == TRUE) {
        flgirot <- TRUE
        irot <- irot[seq_len(r)]
        irot[which(is.na(irot))] <- 0L
    }
    else {
        flgirot <- FALSE
    }
    if (missing(iflip) == FALSE && (isTRUE(class(iflip) %in% 
        c("numeric", "integer")) == TRUE)) {
        flgiflp <- TRUE
    }
    else {
        flgiflp <- FALSE
    }
    if (missing(inv) == FALSE && isTRUE(inv == TRUE) == TRUE) {
        oclu <- clu
        k <- 1L
        for (i in rev(seq_len(nlevels(factor(oclu))))) {
            clu[which(levels(factor(oclu))[i] == oclu)] <- as.numeric(levels(factor(oclu))[k])
            k <- k + 1L
        }
        rm(i)
    }
    else {
        inv <- FALSE
    }
    nlst <- vector(mode = "list", length = r)
    for (i in seq_len(r)) {
        nlst[[i]] <- lbs[which(clu == i)]
    }
    rm(i)
    rad <- 1L
    nds <- data.frame(matrix(ncol = 2, nrow = 0))
    for (i in seq_len(length(nlst))) {
        nl <- length(nlst[[i]])
        if (isTRUE(i == 1L) == TRUE && isTRUE(nl == 1L) == TRUE) {
            x <- data.frame(X = 0L, Y = 0L)
            rad <- rad - 1L
        }
        else {
            x <- data.frame(X = cos(2L * pi * ((0:((nl) - 1L))/nl)) * 
                rad, Y = sin(2L * pi * ((0:(nl - 1L))/nl)) * 
                rad)
        }
        ifelse(isTRUE(flgirot == TRUE) == TRUE, x[, 2:1] <- xyrt(x[, 
            2:1], (irot[i] * -1L)), NA)
        ifelse(isTRUE(flgiflp == TRUE) == TRUE && isTRUE(i %in% 
            iflip) == TRUE, x[, 1] <- x[, 1] * cos(pi) - x[, 
            2] * sin(pi), NA)
        nds <- rbind(nds, x[, 2:1])
        rad <- rad + 1L
    }
    rm(i)
    if (isTRUE(r > 1L) == TRUE && isTRUE(inv == TRUE) == TRUE) {
        nnds <- data.frame(matrix(ncol = ncol(nds), nrow = nrow(nds)))
        nnds[which(lbs %in% nlst[[1]]), ] <- nds[seq_len(length(nlst[[1]])), 
            ]
        yp <- tabulate(clu)[1]
        for (i in seq_len(r - 1L)) {
            k <- 1L
            for (j in which(lbs %in% nlst[[i + 1L]])) {
                nnds[j, ] <- nds[(yp + 1L):(yp + tabulate(clu)[i + 
                  1L]), ][k, ]
                k <- k + 1L
            }
            rm(j)
            yp <- yp + tabulate(clu)[i + 1L]
        }
        rm(i)
        nds <- nnds
    }
    else {
        NA
    }
    switch(match.arg(mirror), N = {
        NA
    }, X = {
        nds[, 1] <- nds[, 1] * cos(pi) - nds[, 2] * sin(pi)
    }, Y = {
        nds[, 2] <- nds[, 2] * cos(pi) - nds[, 1] * sin(pi)
    }, D = {
        nds[, 1:2] <- xyrt(nds[, 1:2], as.numeric(-45))
        nds[, 1:2] <- nds[, 1:2] - min(nds[, 1:2])
        nds[, 2] <- nds[, 2] * cos(pi) - nds[, 1] * sin(pi)
        nds[, 1:2] <- xyrt(nds[, 1:2], as.numeric(45))
        nds[, 1:2] <- nds[, 1:2] - min(nds[, 1:2])
    }, L = {
        nds[, 1:2] <- xyrt(nds[, 1:2], as.numeric(45))
        nds[, 1:2] <- nds[, 1:2] - min(nds[, 1:2])
        nds[, 2] <- nds[, 2] * cos(pi) - nds[, 1] * sin(pi)
        nds[, 1:2] <- xyrt(nds[, 1:2], as.numeric(-45))
        nds[, 1:2] <- nds[, 1:2] - min(nds[, 1:2])
    })
    nds
}

