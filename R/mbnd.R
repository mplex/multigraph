## (GPL-3)  Antonio Rivero Ostoic, multiplex@post.com
mbnd <-
function (pars, r, b, vlt, cx, lw, ecol, directed, bw, alfa, 
    fds, flgcx, valued, flgcr, hds, n) 
{
    ifelse(isTRUE(nrow(pars) > 2L) == TRUE, pars <- pars[c(1, 
        nrow(pars)), ], NA)
    ifelse(isTRUE(ncol(pars) > 2L) == TRUE, pars <- pars[, 1:2], 
        NA)
    angp <- atan2((pars[2, 2] - pars[1, 2]), (pars[2, 1] - pars[1, 
        1])) * (180L/pi)
    if (isTRUE(pars[1, 1] != 0L | pars[1, 2] != 0L) == TRUE) {
        xo <- 0L - pars[1, 1]
        yo <- 0L - pars[1, 2]
        orig <- pars
        orig[, 1] <- pars[, 1] + xo
        orig[, 2] <- pars[, 2] + yo
        orot <- xyrtb(orig, (0L - angp))
    }
    else if (isTRUE(pars[1, 1] != 0L | pars[1, 2] != 0L) == FALSE) {
        orot <- xyrtb(pars, (0L - angp))
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
    if (isTRUE(valued == FALSE) == TRUE) {
        if (isTRUE(mean(lw) < 6L) == TRUE) {
            fds <- fds - (mean(lw) * 4)
        }
        else if (isTRUE(mean(lw) > 15L) == TRUE) {
            fds <- fds - (mean(lw) * 2)
        }
        else {
            fds <- fds - (mean(lw) * 2.5)
        }
    }
    ifelse(isTRUE(valued == TRUE) == TRUE && isTRUE(mean(cx) < 
        2L) == TRUE, cx <- cx + (2L - mean(cx)), NA)
    d <- (rng(r) * ((bw * 1000L) * ((2L^(abs(sin(angp * (pi/180L)))))/1200L)) * 
        (mean(cx)/2L)) * ((30L + (agx * (6L/45L) * -1L))/100L)
    orott <- orot
    if (isTRUE(flgcx == TRUE) == TRUE) {
        orott[1, 1] <- (cx[1]/(fds - 0)) - orot[1, 1]
        orott[2, 1] <- orot[2, 1] - (cx[2]/(fds))
    }
    else if (isTRUE(flgcx == FALSE) == TRUE) {
        ifelse(isTRUE(valued == TRUE) == TRUE && isTRUE(mean(lw) > 
            2L) == TRUE, lc <- (max(lw) * 0.1), lc <- 0)
        orott[1, 1] <- ((cx[1] + lc)/(fds - 0)) - orot[1, 1]
        orott[2, 1] <- orot[2, 1] - ((cx[1] + lc)/fds)
    }
    lst <- array(0L, dim = c(2, 2, r))
    dat <- data.frame(matrix(nrow = 0L, ncol = 2L))
    for (j in 1:r) {
        lst[, 1, j] <- orott[, 1]
        lst[, 2, j] <- d[j]
        dat[(nrow(dat) + 1L):(nrow(dat) + 2L), ] <- lst[, , j]
    }
    rm(j)
    rrot <- xyrt(dat, as.numeric(angp))
    if (isTRUE(pars[1, 1] != 0L | pars[1, 2] != 0L) == TRUE) {
        rrot[, 1] <- rrot[, 1] - xo
        rrot[, 2] <- rrot[, 2] - yo
    }
    for (i in 1:r) {
        if (isTRUE(valued == TRUE) == TRUE && isTRUE(directed == 
            TRUE) == TRUE) {
            if (isTRUE(b[i] %in% multiplex::men(b)[1]) == TRUE) {
                graphics::arrows(rrot[which(seq_len(nrow(dat))%%2L == 
                  1L)[i], 1], rrot[which(seq_len(nrow(dat))%%2L == 
                  1L)[i], 2], rrot[which(seq_len(nrow(dat))%%2L == 
                  0L)[i], 1], rrot[which(seq_len(nrow(dat))%%2L == 
                  0L)[i], 2], code = 2, length = 0, angle = 0, 
                  lty = vlt[i], lw = lw[i], col = grDevices::adjustcolor(ecol[i], 
                    alpha = alfa))
            }
            else if (isTRUE(b[i] %in% multiplex::men(b)[1]) == 
                FALSE) {
                graphics::arrows(rrot[which(seq_len(nrow(dat))%%2L == 
                  1L)[i], 1], rrot[which(seq_len(nrow(dat))%%2L == 
                  1L)[i], 2], rrot[which(seq_len(nrow(dat))%%2L == 
                  0L)[i], 1], rrot[which(seq_len(nrow(dat))%%2L == 
                  0L)[i], 2], code = 1, length = 0, angle = 0, 
                  lty = vlt[i], lw = lw[i], col = grDevices::adjustcolor(ecol[i], 
                    alpha = alfa))
            }
        }
        else {
            graphics::segments(rrot[which(seq_len(nrow(dat))%%2L == 
                1L)[i], 1], rrot[which(seq_len(nrow(dat))%%2L == 
                1L)[i], 2], rrot[which(seq_len(nrow(dat))%%2L == 
                0L)[i], 1], rrot[which(seq_len(nrow(dat))%%2L == 
                0L)[i], 2], lty = vlt[i], lw = lw[i], col = grDevices::adjustcolor(ecol[i], 
                alpha = alfa))
        }
        if (isTRUE(directed == TRUE) == TRUE) {
            if (isTRUE(n < 15) == TRUE) {
                Hd <- data.frame(x = c((-0.6 - (0.01538462 * 
                  n)), (-0.35 - (0.01538462 * n)), (-0.6 - (0.01538462 * 
                  n)), (0.4 - (0.01538462 * n))), y = c(-0.5, 
                  0, 0.5, 0)) * (hds)
            }
            else {
                Hd <- data.frame(x = c(-0.7, -0.45, -0.7, 0.3), 
                  y = c(-0.5, 0, 0.5, 0)) * (hds)
            }
            if (isTRUE(as.numeric(lw[i]) < 7L) == TRUE) {
                ifelse(isTRUE(as.numeric(lw[i]) <= 1L) == TRUE, 
                  Hd <- Hd * (as.numeric(lw[i]))/((as.numeric(lw[i]) * 
                    8.571) + 30), Hd <- Hd * (as.numeric(lw[i]))/((as.numeric(lw[i]) * 
                    8.571) + 40))
            }
            else if (isTRUE(as.numeric(lw[i]) >= 15L) == TRUE) {
                Hd <- Hd * (as.numeric(lw[i]))/(as.numeric(lw[i]) + 
                  120L)
            }
            else {
                Hd <- Hd * (as.numeric(lw[i]))/((as.numeric(lw[i]) * 
                  8.571) + 40)
            }
            if (isTRUE(i %in% flgcr) == TRUE) {
                prx1 <- rrot[which(seq_len(nrow(dat))%%2L == 
                  1L)[i], 1]
                pry1 <- rrot[which(seq_len(nrow(dat))%%2L == 
                  1L)[i], 2]
                hd1 <- xyrt((Hd), (as.numeric(angp) - 180L))
                hd1[, 1] <- hd1[, 1] + prx1
                hd1[, 2] <- hd1[, 2] + pry1
                graphics::polygon((hd1), col = grDevices::adjustcolor(ecol[i], 
                  alpha = alfa), border = NA)
                prx2 <- rrot[which(seq_len(nrow(dat))%%2L == 
                  0L)[i], 1]
                pry2 <- rrot[which(seq_len(nrow(dat))%%2L == 
                  0L)[i], 2]
                hd2 <- xyrt(Hd, (as.numeric(angp) - 0L))
                hd2[, 1] <- hd2[, 1] + prx2
                hd2[, 2] <- hd2[, 2] + pry2
                graphics::polygon((hd2), col = grDevices::adjustcolor(ecol[i], 
                  alpha = alfa), border = NA)
            }
            else {
                if (isTRUE(b[i] %in% multiplex::men(b)[1]) == 
                  FALSE) {
                  prx <- rrot[which(seq_len(nrow(dat))%%2L == 
                    1L)[i], 1]
                  pry <- rrot[which(seq_len(nrow(dat))%%2L == 
                    1L)[i], 2]
                  hd <- xyrt((Hd), (as.numeric(angp) - 180L))
                }
                else if (isTRUE(b[i] %in% multiplex::men(b)[1]) == 
                  TRUE) {
                  prx <- rrot[which(seq_len(nrow(dat))%%2L == 
                    0L)[i], 1]
                  pry <- rrot[which(seq_len(nrow(dat))%%2L == 
                    0L)[i], 2]
                  hd <- xyrt(Hd, (as.numeric(angp) - 0L))
                }
                hd[, 1] <- hd[, 1] + prx
                hd[, 2] <- hd[, 2] + pry
                graphics::polygon((hd), col = grDevices::adjustcolor(ecol[i], 
                  alpha = alfa), border = NA)
            }
        }
        else {
            invisible(NA)
        }
    }
    rm(i)
    x <- NULL
    rm(x)
    graphics::par(new = FALSE)
}

