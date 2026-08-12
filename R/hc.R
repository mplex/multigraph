## (GPL-3)  Antonio Rivero Ostoic, multiplex@post.com
hc <-
function (x, y, r, nsteps = 900, flip, ...) 
{
    rs <- seq(-0.75, pi * 1.25, len = nsteps)
    if (isTRUE(flip == TRUE) == TRUE) {
        hcr <- data.frame(xc = x - r * cos(rs), yc = y - r * 
            sin(rs))
        graphics::lines(hcr[, 1], (hcr[, 2] - r), ...)
    }
    else {
        hcr <- data.frame(xc = x + r * cos(rs), yc = y + r * 
            sin(rs))
        graphics::lines(hcr[, 1], hcr[, 2], ...)
    }
    graphics::par(new = FALSE)
}

