rng <-
function (r) 
{
    if (r == 1L) 
        return(0)
    if (r > 1L) {
        x <- vector()
        x <- append(x, (-1))
        for (i in 1:(r - 1)) x <- append(x, ((-1) + (2L/(r - 
            1L)) * i))
        return(x * (r/50L))
    }
    else stop("no negative values")
}
