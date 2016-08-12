xyrtb <-
function (pares, ang) 
{
    if (isTRUE(pares[1, 1] != 0 | pares[2, 1] != 0) == TRUE) {
        temp <- pares
        tpares <- c(0, 0) - pares[1, ]
        for (i in 1:nrow(temp)) {
            temp[i, ] <- pares[i, ] + tpares
        }
        rm(i)
        txrot <- temp[, 1] * cos(ang * (pi/180L)) - temp[, 2] * 
            sin(ang * (pi/180L))
        tyrot <- temp[, 2] * cos(ang * (pi/180L)) + temp[, 1] * 
            sin(ang * (pi/180L))
        xrot <- txrot
        yrot <- tyrot
        tm <- as.data.frame(cbind(xrot, yrot))
        tm[1, ] <- tm[1, ] - (c(0, 0) - pares[1, ])
        tm[2, ] <- tm[2, ] - (c(0, 0) - pares[1, ])
    }
    else {
        xrot <- pares[, 1] * cos(ang * (pi/180L)) - pares[, 2] * 
            sin(ang * (pi/180L))
        yrot <- pares[, 2] * cos(ang * (pi/180L)) + pares[, 1] * 
            sin(ang * (pi/180L))
        tm <- cbind(xrot, yrot)
    }
    attr(tm, "dimnames") <- NULL
    return(invisible(tm))
}
