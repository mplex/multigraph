conc <-
function (net, nr, irot, inv, flip, ...) 
{
    n <- dim(net)[1]
    if (isTRUE(is.null(dimnames(net)[1]) == TRUE | is.null(dimnames(net)[1][[1]])) == 
        FALSE) {
        lbs <- dimnames(net)[[1]]
    }
    else {
        lbs <- 1:n
    }
    if (missing(nr) == TRUE) {
        clu <- rep(1, n)
    }
    else if (isTRUE(length(nr) == 1) == TRUE && is.numeric(nr) == 
        TRUE) {
        if (isTRUE(nr >= n) == TRUE | isTRUE(nr <= 0) == TRUE) 
            stop("Value of 'nr' must be greater than zero and lower than network order.")
        clu <- vector()
        for (i in 1:nr) {
            clu <- append(clu, rep(i, ceiling(n/nr)))
        }
        rm(i)
        clu <- clu[1:n]
    }
    else if (is.factor(nr) == TRUE) {
        tmpnr <- nr[1:n]
        if (any(is.na(tmpnr)) == TRUE) 
            stop("Insufficient length of 'nr'.")
        for (i in 1:nlevels(factor(nr))) {
            levels(nr) <- c(levels(nr), i)
            nr[which(levels(factor(tmpnr))[i] == nr)] <- i
        }
        rm(i)
        clu <- methods::as(as.vector(nr), "numeric")
        rm(tmpnr)
    }
    else if (is.character(nr) == TRUE) {
        tmpnr <- nr[1:n]
        if (any(is.na(tmpnr)) == TRUE) 
            stop("Insufficient length of 'nr'.")
        nr[which(nr == nr[1])] <- 1
        for (i in 1:(nlevels(factor(tmpnr)) - 1L)) {
            nr[which(nr == nr[which((nr %in% tmpnr) == TRUE)[(i - 
                0)]])] <- (i + 1)
        }
        rm(i)
        nr[which((nr %in% tmpnr) == TRUE)] <- nlevels(factor(tmpnr))
        clu <- methods::as(as.vector(nr), "numeric")
        rm(tmpnr)
    }
    else if (is.numeric(nr) == TRUE) {
        ifelse(is.integer(nr) == TRUE, clu <- nr[1:n], clu <- methods::as(as.vector(nr[1:n]), 
            "integer"))
        clu[which(is.na(clu))] <- 0
        ifelse(isTRUE(0 %in% (clu)) == TRUE, clu <- clu + 1L, 
            NA)
        tmpclu <- clu
        for (i in 1:nlevels(factor(clu))) {
            clu[which(levels(factor(tmpclu))[i] == clu)] <- i
        }
        rm(i)
        rm(tmpclu)
    }
    else {
        NA
    }
    r <- nlevels(factor(clu))
    if (missing(irot) == FALSE && is.numeric(irot) == TRUE) {
        flgirot <- TRUE
        irot <- irot[1:r]
        irot[which(is.na(irot))] <- 0
    }
    else {
        flgirot <- FALSE
    }
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
        ifelse(isTRUE(flgirot == TRUE) == TRUE, x[, 2:1] <- xyrt(x[, 
            2:1], (irot[i] * -1L)), NA)
        if (missing(flip) == FALSE && isTRUE(flip == TRUE) == 
            TRUE) {
            ifelse(isTRUE((i%%2L) == 0) == TRUE, x[, 1] <- x[, 
                1] * cos(pi) - x[, 2] * sin(pi), NA)
        }
        else {
            NA
        }
        nds <- rbind(nds, x[, 2:1])
        rad <- rad + 1L
    }
    rm(i)
    if (isTRUE(r > 1) == TRUE) {
        nnds <- data.frame(matrix(ncol = ncol(nds), nrow = nrow(nds)))
        nnds[which(lbs %in% nlst[[1]]), ] <- nds[1:length(nlst[[1]]), 
            ]
        yp <- tabulate(clu)[1]
        for (i in 1:(r - 1)) {
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
        nnds
    }
    else {
        nds
    }
}
