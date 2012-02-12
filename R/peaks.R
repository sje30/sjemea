##  https://stat.ethz.ch/pipermail/r-help/2005-November/083376.html
## Martin Machler, Nov 2005.

peaks <- function(series, span = 3, do.pad = TRUE) {
    if((span <- as.integer(span)) %% 2 != 1) stop("'span' must be odd")
    s1 <- 1:1 + (s <- span %/% 2)
    if(span == 1) return(rep.int(TRUE, length(series)))
    z <- embed(series, span)
    v <- apply(z[,s1] > z[, -s1, drop=FALSE], 1, all)
    if(do.pad) {
        pad <- rep.int(FALSE, s)
        c(pad, v, pad)
    } else v
}

peaksign <- function(series, span = 3, do.pad = TRUE)
{
    ## Purpose: return (-1 / 0 / 1) if series[i] is ( trough / "normal" / peak )
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 25 Nov 2005

    if((span <- as.integer(span)) %% 2 != 1 || span == 1)
        stop("'span' must be odd and >= 3")
    s1 <- 1:1 + (s <- span %/% 2)
    z <- embed(series, span)
    d <- z[,s1] - z[, -s1, drop=FALSE]
    ans <- rep.int(0:0, nrow(d))
    ans[apply(d > 0, 1, all)] <- as.integer(1)
    ans[apply(d < 0, 1, all)] <- as.integer(-1)
    if(do.pad) {
        pad <- rep.int(0:0, s)
        c(pad, ans, pad)
    } else ans
}


check.pks <- function(y, span = 3) {
  stopifnot(identical(peaks( y, span), peaksign(y, span) ==  1),
            identical(peaks(-y, span), peaksign(y, span) == -1))
}

for(y in list(1:10, rep(1,10), c(11,2,2,3,4,4,6,6,6))) {
    for(sp in c(3,5,7))
        check.pks(y, span = sp)
    stopifnot(peaksign(y) == 0)
}

y <- c(1,4,1,1,6,1,5,1,1) ; (ii <- which(peaks(y))); y[ii]
##- [1] 2 5 7
##- [1] 4 6 5
check.pks(y)

set.seed(7)
y <- rpois(100, lambda = 7)
check.pks(y)
py <- peaks(y)
plot(y, type="o", cex = 1/4, main = "y and peaks(y,3)")
points(seq(y)[py], y[py], col = 2, cex = 1.5)

p7 <- peaks(y,7)
points(seq(y)[p7], y[p7], col = 3, cex = 2)
mtext("peaks(y,7)", col=3)

set.seed(2)
x <- round(rnorm(500), 2)
y <- cumsum(x)
check.pks(y)

plot(y, type="o", cex = 1/4)
p15 <- peaks(y,15)
points(seq(y)[p15], y[p15], col = 3, cex = 2)
mtext("peaks(y,15)", col=3)
