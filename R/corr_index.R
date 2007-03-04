## Compute the correlation index.
## Code is to be cleaned up a bit.
## Sun 04 Mar 2007

corr.index <- function(s, distance.breaks, dt=0.05) {
  ## Make a correlation index object.
  dists = make.distances(s$layout$pos)
  dists.bins = bin.distances(dists, distance.breaks)

  spikes = s$spikes
  if (length(spikes) > 1) {
    corr.indexes = make.corr.indexes(spikes, dt)
    corr.id = cbind(my.upper(dists), my.upper(corr.indexes))
    corr.id.means = corr.get.means(corr.id)
  } else {
    corr.indexes = NA
    corr.id = NA
    corr.id.means = NA
  }

  res = list(dists=dists, dists.bins = dists.bins,
    corr.indexes = corr.indexes,
    dt = dt,
    corr.id = corr.id,
    corr.id.means = corr.id.means)

  res
}


plot.corr.index <- function(s, log='', identify=FALSE,
                            main=NULL,
                            ...) {
  ## Plot the correlation indices as a function of distance.
  ## If identify is TRUE, we can locate cell pairs on the plot using
  ## left mouse button.
  
  dists = s$corr$dists[which(upper.tri(s$corr$dists))]
  corrs = s$corr$corr.indexes[which(upper.tri(s$corr$corr.indexes))]

  if (is.null(main)) {
    main = paste(basename(s$file), "dt:", s$corr$dt)
  }
  
  xlabel = expression(paste("distance (", mu, "m)"))

  plot.default(dists, corrs, xlab=xlabel, ##log=log,
               ylab="correlation index", bty="n",
               main=main, col='red',
               ...)

  if (identify) {
    labels1 <- outer(seq(1, s$NCells), seq(1,s$NCells), FUN="paste")
    labs <- labels1[which(upper.tri(labels1))]
    identify(dists, corrs, labels=labs)
  }

  plotCI(s$corr$corr.id.means[,1], s$corr$corr.id.means[,2],
         s$corr$corr.id.means[,3],
         xlab="distance", ylab="correlation index", 
         pch=19, add=TRUE)
  corr.do.fit(s$corr$corr.id,plot=TRUE)
}

plot.corr.index.fit <- function(s, ...) {
  ### Show the correlation indexes and then the fit.
  plot.corr.index(s, identify=FALSE,col="red", log="")
  plotCI(s$corr.id.means[,1], s$corr.id.means[,2], s$corr.id.means[,3],
         xlab="distance", ylab="correlation index", 
         pch=19, add=TRUE)
  corr.do.fit(s$corr.id,plot=TRUE)
}


make.corr.indexes <- function(spikes, dt) {
  ## Return the correlation index values for each pair of spikes.
  ## The matrix returned is upper triangular.
  ## SPIKES should be a list of length N, N is the number of cells.
  ## "dt" is the maximum time for seeing whether two spikes are coincident.
  ## This is defined in the 1991 Meister paper.
  n <- length(spikes)
  if (n == 1) {
    ## If only one spike train, cannot compute the cross-corr indexes.
    0;
  } else {
    Tmax <- max(unlist(spikes))           #time of last spike.
    Tmin <- min(unlist(spikes))           #time of first spike.
    corrs <- array(0, dim=c(n,n))
    for ( a in 1:(n-1)) {
      n1 <- length(spikes[[a]])
      for (b in (a+1):n) {
        n2 <- length(spikes[[b]])
        corrs[a,b] <-  as.double(count.nab(spikes[[a]], spikes[[b]],dt) *
                                 (Tmax-Tmin)) /
                                   (as.double(n1) * n2 * (2*dt))
      }
    }
    if (any(is.na(corrs))) {
      stop("corrs has some NA values -- possible integer overflow in n1*n2, or zero spikes in one of the trains?")
    }
    corrs
  }
}



corr.index.means <- function(x) {
  ## Compute the mean,sd correlation index at each given distance.
  dists <- x$dists[which(upper.tri(x$dists))]
  corrs <- x$corr.indexes[which(upper.tri(x$corr.indexes))]

  dists.uniq <- unique(dists)
  num.dists <- length(dists.uniq)       #num of  different distances.

  ##print(dists.uniq)
  ## create 4-D array to store results.  Each row stores the
  ## distance, mean corr, sd, and num of values at that distance.

  res <- array(0,  dim=c(num.dists,4))
  colnames(res) <- c("dist","mean corr", "sd", "n")
  
  i <- 1

  for (d in dists.uniq) {
    ## find all correlations for pairs within 0.01um of given distance.
    cs <- corrs[ which(abs(dists-d)<0.01)]
    corrs.mean <- mean(cs)
    corrs.sd   <- sd(cs)
    res[i,] <- c(d, corrs.mean, corrs.sd, length(cs))
    i <- 1+i
  }

  res
}

corr.get.means <- function(id) {
  ## Compute the mean,sd of the correlation index at each distance.
  ## id is the array of [n,2] values.  Each row is [d,i].
  ## Returns a matrix.
  
  corr.get.means.helper <- function(x) {
    ## Helper function to create  mean and sd of one set of distances.
    indexes <- which(id[,1] == x)
    c(x, mean(id[indexes,2]), sd(id[indexes,2]), length(indexes))
    ##c(x, median(id[indexes,2]), mad(id[indexes,2]), length(indexes))
  }
  
  d.uniq <- sort(unique(id[,1]))
  means <- t(sapply(d.uniq, corr.get.means.helper))
  colnames(means) <- c("dist", "mean", "sd", "n")
  means
}

corr.do.fit <- function(id, plot=TRUE, show.ci=FALSE, ...) {
  ## Do the fit to the exponential and optionally plot it.  Any
  ## correlation index of zero is removed, since we cannot take the
  ## log of zero.  Hopefully there won't be too many of these.
  ## If SHOW.CI is true, do the fit with 95% confidence intervals.

  y.zero <- which(id[,2]==0)
  if (length(y.zero)>0) {
    id <- id[-y.zero,]
    warning(paste("removing", length(y.zero),"zero entries"))
  }
  x <- id[,1]
  y.log <- log(id[,2])
  fit <- lm(y.log ~ x)
  if (show.ci) {
    expt.new <- data.frame(x = seq(0, 850, 10))  #range of x for predictions.
    expt.clim <- predict(fit, expt.new, interval="confidence")
  }
  if (plot)  {
    if (show.ci) {
      ## Confidence intervals will show mean, so don't need
      ## to do both matlines and curve.
      matlines(expt.new$x, exp(expt.clim), lty=c(1,2,2),
               col="black")
    } else {
      curve(exp(fit$coeff[1]) * exp(x * fit$coeff[2]), add = TRUE,
            from=0, ...)
    }
  }
  fit
}

corr.check.fit <- function() {
  ## Simple test routine to see that the exponential fits are okay.
  a <- 40; b <- .01
  x <- seq(from=1, to=500, by=20)
  y <- a*exp(-b*x) + (2*rnorm(length(x)))
  plot(x,y, log="y")
  fit <- corr.do.fit( cbind(x,y), col=p9.col)
  
  ## should be similar to (a,b)
  print(exp(fit$coefficients))
}

my.upper <- function (x,diag=FALSE) {
  ## Return the upper triangular elements of a matrix on a
  ## column-by-column basis.
  ## e.g. my.upper(matrix(1:9, nrow=3), diag=TRUE)
  ## returns >>1 4 5 7 8 9<<
  if (is.matrix(x)) {
   x[ which(upper.tri(x,diag))]
  } else {
    stop(paste(deparse(substitute(x)),"is not a matrix"))
  }
}
