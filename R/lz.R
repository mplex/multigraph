lz <-
function (x, delta, w) 
{
    ifelse(missing(delta) == TRUE, delta <- matrix(1, nrow(x), 
        nrow(x)), NA)
    L <- matrix(0, nrow = attr(w, "Size"), ncol = attr(w, "Size"))
    for (i in 1:nrow(x)) {
        D <- 0
        for (j in 1:nrow(x)) {
            if (isTRUE(i != j) == TRUE) {
                nrmz <- norm(x[i, ] - x[j, ], type = "2")
                delt <- as.matrix(w)[i, j] * as.matrix(delta)[i, 
                  j]
                L[i, j] <- ((-delt)/nrmz)
                D <- D - ((-delt)/nrmz)
            }
            else {
                nrmz <- 0
            }
        }
        rm(j)
        L[i, i] <- D
    }
    rm(i)
    L[which(L == Inf)] <- 0
    L
}
