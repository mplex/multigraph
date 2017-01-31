sts <-
function (nds, delta = NULL, mwd = NULL) 
{
    ifelse(is.null(delta) == TRUE, delta <- matrix(1, nrow(nds), 
        nrow(nds)), NA)
    strss <- 0
    mnds <- as.matrix(nds)
    if (is.null(mwd)) {
        mwd <- delta^(-2)
        mwd[which(mwd == Inf)] <- 0
    }
    for (j in 1:nrow(nds)) {
        for (i in 1:(j - 1)) {
            if (isTRUE(i > 0 && j > 1) == TRUE) {
                strss <- strss + as.matrix(mwd)[i, j] %*% (norm(mnds[i, 
                  ] - mnds[j, ], type = "2") - as.matrix(delta)[i, 
                  j])^2
            }
            else {
                NA
            }
        }
        rm(i)
    }
    rm(j)
    as.vector(strss)
}
