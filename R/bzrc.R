bzrc <-
function (pair, cex, elv = 0.25, lng, ...) 
{
    ifelse(missing(cex) == TRUE, cex <- 1L, NA)
    ifelse(missing(lng) == TRUE, lng <- 50, NA)
    ifelse(isTRUE(lng <= 2) == TRUE, lng <- 3L, NA)
    if (isTRUE(max(round(pair[, 1], 3)) < 0.5) == TRUE) {
        p <- rbind(pair[1, ], c(elv * -1, mean(c(pair[1, 2], 
            pair[2, 2]))), pair[2, ])
    }
    else {
        p <- rbind(pair[1, ], c(1L + elv, mean(c(pair[1, 2], 
            pair[2, 2]))), pair[2, ])
    }
    attr(p, "dimnames")[[2]] <- NULL
    d <- (cex/100)
    if (isTRUE(as.numeric(p[1, 2] - p[3, 2]) > d) == TRUE) {
        d <- d * (0.1/as.numeric(p[1, 2] - p[3, 2]))
    }
    else {
        NA
    }
    sqq <- seq(0, 1, length = lng) - d
    sq <- sqq[which(sqq > d)]
    deg <- nrow(p) - 1
    num_s <- (nrow(p) - 1)/deg
    if (isTRUE(num_s - floor(num_s) > 0) == TRUE) 
        stop("Number of rows in parameter matrix do not match input degree.")
    b <- matrix(0, nrow = length(sq), ncol = ncol(p))
    seg <- matrix(1, nrow = num_s, ncol = 2)
    for (i in 1:num_s) seg[i, ] <- c((i - 1) * deg + 1, (i - 
        1) * deg + 1 + deg)
    for (j in 1:length(sq)) {
        if (sq[j] == 0) {
            s <- 1
        }
        else {
            s <- ceiling(sq[j])
        }
        if (s > num_s) 
            s <- num_s
        p_sub <- matrix(p[seg[s, 1]:seg[s, 2], ], nrow = deg + 
            1, ncol = ncol(p))
        b[j, ] <- colSums(choose(deg, 0:deg) * ((1 - (sq[j] - 
            s + 1))^(deg - 0:deg)) * (sq[j] - s + 1)^(0:deg) * 
            p_sub)
    }
    rm(j)
    graphics::lines(b[, 1], b[, 2], ...)
    graphics::par(new = FALSE)
}
