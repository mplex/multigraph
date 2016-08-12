area <-
function (x, solo = TRUE) 
{
    if (isTRUE(solo == TRUE) == TRUE) {
        return(sqrt(((max(x[, 1]) - min(x[, 1])) * (max(x[, 2]) - 
            min(x[, 2])))/nrow(x)))
    }
    else if (isTRUE(solo == TRUE) == FALSE) {
        cds <- x
        rat <- (max(cds[, 1]) - min(cds[, 1]))/(max(cds[, 2]) - 
            min(cds[, 2]))
        cds[, 1] <- (cds[, 1] - min(cds[, 1]))/(max(cds[, 1]) - 
            min(cds[, 1]))
        ifelse(isTRUE(rat > 0) == TRUE, cds[, 2] <- ((cds[, 2] - 
            min(cds[, 2]))/(max(cds[, 2]) - min(cds[, 2]))) * 
            (1L/rat), cds[, 2] <- ((cds[, 2] - min(cds[, 2]))/(max(cds[, 
            2]) - min(cds[, 2]))) * (rat))
        nds <- data.frame(X = as.numeric(as.vector(cds[, 1])), 
            Y = as.numeric(as.vector(cds[, 2])))
        nds <- (2L/max(nds)) * nds
        return(sqrt(((max(nds[, 1]) - min(nds[, 1])) * (max(nds[, 
            2]) - min(nds[, 2])))/nrow(nds)))
    }
}
