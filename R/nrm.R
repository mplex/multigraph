nrm <-
function (x, digits = 3) 
{
    if (is.array(x) == TRUE) {
        xnorm <- (x[, 1] - min(x[, 1]))/(max(x[, 1]) - min(x[, 
            1]))
        rat <- (max(x[, 1]) - min(x[, 1]))/(max(x[, 2]) - min(x[, 
            2]))
        ynorm <- ((x[, 2] - min(x[, 2]))/(max(x[, 2]) - min(x[, 
            2]))) * (rat)
        ifelse(isTRUE(rat > 0) == FALSE, ynorm <- ((x[, 2] - 
            min(x[, 2]))/(max(x[, 2]) - min(x[, 2]))) * (1L/rat), 
            NA)
        return(round(data.frame(X = xnorm, Y = ynorm), digits))
    }
    else if (is.vector(x) == TRUE) {
        return(round((x - min(x))/(max(x) - min(x)), digits))
    }
}
