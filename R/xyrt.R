xyrt <-
function (pares, ang) 
{
    xr <- pares[, 1] * cos(ang * (pi/180L)) - pares[, 2] * sin(ang * 
        (pi/180L))
    yr <- pares[, 2] * cos(ang * (pi/180L)) + pares[, 1] * sin(ang * 
        (pi/180L))
    return(invisible(cbind(xr, yr)))
}
