mbnd <-
function (pares, r, b, vlt, cx, lwd, ecol, directed, asp, bwd, 
    alfa, fds, flgcx) 
{
    ifelse(isTRUE(nrow(pares) > 2L) == TRUE, pares <- pares[c(1, 
        nrow(pares)), ], NA)
    ifelse(isTRUE(ncol(pares) > 2L) == TRUE, pares <- pares[, 
        1:2], NA)
    angp <- atan2((pares[2, 2] - pares[1, 2]), (pares[2, 1] - 
        pares[1, 1])) * (180L/pi)
    if (isTRUE(pares[1, 1] != 0L | pares[1, 2] != 0L) == TRUE) {
        xo <- 0L - pares[1, 1]
        yo <- 0L - pares[1, 2]
        orig <- pares
        orig[, 1] <- pares[, 1] + xo
        orig[, 2] <- pares[, 2] + yo
        orot <- xyrtb(orig, (0L - angp))
    }
    else if (isTRUE(pares[1, 1] != 0L | pares[1, 2] != 0L) == 
        FALSE) {
        orot <- xyrtb(pares, (0L - angp))
    }
    angx <- abs(angp)
    if (isTRUE(angx >= 360L) == TRUE) {
        agx <- angx%%360L
    }
    else {
        if (isTRUE(angx >= 270L) == TRUE) {
            gx <- abs(90L + (angx%%270L))
        }
        else if (isTRUE(angx > 180L) == TRUE) {
            ifelse(isTRUE(angx > 180L) == TRUE, lx <- abs(180L - 
                (angx%%180L)), lx <- angx%%180L)
            gx <- abs(90L - (lx%%90L))
        }
        else {
            ifelse(isTRUE(angx > 180L) == TRUE, gx <- abs(180L - 
                (x%%180L)), gx <- angx%%180L)
        }
        ifelse(isTRUE(gx >= 90L) == TRUE, agx <- abs(90L - (gx%%90L)), 
            agx <- gx%%90L)
    }
    qx <- ((30L + (agx * (6L/45L) * -1L))/100L)
    d <- (rng(r) * ((bwd * 1000L) * ((2L^(abs(sin(angp * (pi/180L)))))/1200L)) * 
        (mean(cx)/2L)) * qx
    ifelse(isTRUE(lwd > 1) == TRUE, fds <- fds - (lwd * 3L), 
        NA)
    orott <- orot
    orott[1, 1] <- (cx[1]/fds) - orot[1, 1]
    ifelse(isTRUE(flgcx == TRUE) == TRUE, orott[2, 1] <- orot[2, 
        1] - (cx[2]/(fds + 20)), orott[2, 1] <- orot[2, 1] - 
        (cx[1]/fds))
    lst <- array(0L, dim = c(2, 2, r))
    dat <- data.frame(matrix(nrow = 0L, ncol = 2L))
    for (i in 1:r) {
        lst[, 1, i] <- orott[, 1]
        lst[, 2, i] <- d[i]
        dat[(nrow(dat) + 1L):(nrow(dat) + 2L), ] <- lst[, , i]
    }
    rm(i)
    rrot <- xyrt(dat, as.numeric(angp))
    if (isTRUE(pares[1, 1] != 0L | pares[1, 2] != 0L) == TRUE) {
        rrot[, 1] <- rrot[, 1] - xo
        rrot[, 2] <- rrot[, 2] - yo
    }
    for (i in 1:r) {
        graphics::segments(rrot[which(seq(1:nrow(dat))%%2L == 
            1L)[i], 1], rrot[which(seq(1:nrow(dat))%%2L == 1L)[i], 
            2], rrot[which(seq(1:nrow(dat))%%2L == 0L)[i], 1], 
            rrot[which(seq(1:nrow(dat))%%2L == 0L)[i], 2], lty = vlt[i], 
            lwd = lwd, col = grDevices::adjustcolor(ecol[i], 
                alpha = alfa))
        if (isTRUE(directed == TRUE) == TRUE) {
            Hd <- data.frame(x = c(-0.8, -0.55, -0.8, 0.2), y = c(-0.5, 
                0, 0.5, 0))
            Hd <- Hd * (as.numeric(lwd))/60L
            if (isTRUE(b[i] %in% multiplex::men(b)[1]) == FALSE) {
                prx <- rrot[which(seq(1:nrow(dat))%%2L == 1L)[i], 
                  1]
                pry <- rrot[which(seq(1:nrow(dat))%%2L == 1L)[i], 
                  2]
                hd <- xyrt((Hd), (as.numeric(angp) - 180L))
            }
            else if (isTRUE(b[i] %in% multiplex::men(b)[1]) == 
                TRUE) {
                prx <- rrot[which(seq(1:nrow(dat))%%2L == 0L)[i], 
                  1]
                pry <- rrot[which(seq(1:nrow(dat))%%2L == 0L)[i], 
                  2]
                hd <- xyrt(Hd, (as.numeric(angp) - 0L))
            }
            hd[, 1] <- hd[, 1] + prx
            hd[, 2] <- hd[, 2] + pry
            graphics::polygon((hd), col = grDevices::adjustcolor(ecol[i], 
                alpha = alfa), border = NA)
        }
        else {
            NA
        }
    }
    rm(i)
    graphics::par(new = FALSE)
    x <- NULL
    rm(x)
}
