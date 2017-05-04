hc <-
function (x, y, r, nsteps = 900, ...) 
{
    rs <- seq(-0.75, pi * 1.25, len = nsteps)
    xc <- x + r * cos(rs)
    yc <- y + r * sin(rs)
    hcr <- data.frame(x = xc, y = yc)
    graphics::lines(hcr[, 1], hcr[, 2], ...)
    graphics::par(new = FALSE)
}
