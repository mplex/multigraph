popl <-
function (n, seed = seed) 
{
    alpha <- pi * (3L - sqrt(5L))
    set.seed(seed)
    phase <- stats::runif(1) * 2L * pi
    ptsx <- vector()
    ptsy <- vector()
    for (k in 1:n) {
        theta <- k * alpha + phase
        rr <- sqrt(k/n)
        ptsx <- append(ptsx, (rr * cos(theta)))
        ptsy <- append(ptsy, (rr * sin(theta)))
    }
    rm(k)
    return(cbind(ptsx, ptsy) - min(cbind(ptsx, ptsy)))
}
