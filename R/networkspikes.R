## networkspikes.R --- identify and analsyse network spikes
## Sun 28 Jan 2007
## Taking ideas from Eytan & Marom (2006).


##ns.T = 0.003                             #bin time for network spikes
##ns.N = 10                               #number of active electrodes.

## 2007-07-27: Code merged in from second version, temp in
## ~/proj/sangermea/test_ns.R Old versions of code at bottom of this
## file can be deleted in a few months.
compute.ns <- function(s, ns.T, ns.N, sur, plot=FALSE) {
  ## Main entrance function to compute network spikes.
  ## Typical values:
  ## ns.T: 3
  ## ns.N: 10
  ## sur: 100

  counts <- spikes.to.count2(s$spikes, time.interval=ns.T)
  p <- find.peaks(counts, ns.N)
  ns <- list(counts=counts, ns.N=ns.N, ns.T=ns.T)
  class(ns) <- "ns"
  m <- mean.ns(ns, p, plot=plot, nrow=4, ncol=4, ask=F, sur=sur)
  if (is.null(m)) {
    ## No network spikes found.
    ns$brief <- c(n=0, peak.m=NA, peak.sd=NA, durn.m=NA, durn.sd=NA)
  } else {
    ns$mean <- m$ns.mean; ns$measures <- m$measures
    peak.val <- ns$measures[,"peak.val"]
    durn <- ns$measures[,"durn"]
    ns$brief <- c(n=nrow(ns$measures),
                  peak.m=mean(peak.val), peak.sd=sd(peak.val),
                  durn.m=mean(durn), durn.sd=sd(durn))

  }
  
  ns
}

spikes.to.count2 <- function(spikes,
                            time.interval=1, #time bin of 1sec.
                            beg=floor(min(unlist(spikes))),
                            end=ceiling(max(unlist(spikes)))
                            )
{
  ## Convert the spikes for each cell into a firing rate (in Hz)
  ## We count the number of spikes within time bins of duration
  ## time.interval (measured in seconds).
  ##
  ## Currently cannot specify BEG or END as less than the
  ## range of spike times else you get an error from hist().  The
  ## default anyway is to do all the spikes within a data file.
  ##
  ## C version, which should replace spikes.to.count
  ## Returns a time series object.
  
  ## time.breaks <- seq(from=beg, to=end, by=time.interval)
  nbins <- ceiling( (end-beg) / time.interval)
  
  z <- .C("ns_count_activity",
          as.double(unlist(s$spikes)),
          as.integer(s$nspikes),
          as.integer(s$NCells),
          as.double(beg), as.double(end), as.double(time.interval),
          as.integer(nbins),
          counts = integer(nbins),
          PACKAGE="sjemea")

  ## Return counts as a time series.
  res <- ts(data=z$counts, start=beg, deltat=time.interval)

  res
}

plot.ns <- function(ns, ...) {
  ## Plot function for "ns" class.
  plot(ns$counts, ...)
  abline(h=ns$ns.N, col='red')

  ##peak.times <- times[ ns$peaks[,1]]
  peak.times <- ns$measures[,"time"]
  peak.val   <- ns$measures[,"peak.val"]
  points(peak.times, peak.val, col='blue', pch=19)

}

summary.ns <- function(ns) {
  ## Summary function for "ns" class.
  n <- ns$brief["n"]
  cat(sprintf("%d network spikes\n", n))
  peak.m <- ns$brief["peak.m"]
  peak.sd <- ns$brief["peak.sd"]


  durn.m <- ns$brief["durn.m"]
  durn.sd <- ns$brief["durn.sd"]
  cat(sprintf("recruitment %.2f +/- %.2f\n", peak.m, peak.sd))
  cat(sprintf("FWHM %.3f +/- %.3f (s)\n", durn.m, durn.sd))
}

mean.ns <- function(ns, p, sur=100,
                    plot=TRUE, nrow=8, ncol=8, ask=FALSE) {
  ## Compute the mean network spikes, and optionally show the
  ## individual network spikes.

  ## This code does not check to worry if there is a spike right at either
  ## end of the recording.  naughty!

  if (is.null(p)) {
    if (is.null(ns$measures)) {
      cat("*** No network spikes found\n")
      return (NULL)
    } else {
      ## use info previously stored in measures.
      p <- ns$measures
    }
  }

  
  if (plot) {
    old.par <- par(mfrow=c(nrow,ncol), mar=c(2.5,1,1,1),ask=ask)
  }
  ave = rep(0, (2*sur)+1)
  npts = length(ns$counts)
  times <- time(ns$counts)
  measures = matrix(NA, nrow=nrow(p), ncol=4)
  colnames(measures) = c("time", "index", "peak.val", "durn")
  n.ns = 0                              #Number of valid network spikes found
  for (i in 1:nrow(p)) {
    peak.i = p[i,"index"]; lo = (peak.i-sur); hi = peak.i+sur

    ## Check that enough data can be found:
    if ( (lo >0) && ( hi < npts) ) {
      n.ns = n.ns + 1

      dat = ns$counts[lo:hi]
      peak.val = dat[sur+1]
      measures[n.ns, "time"] = times[peak.i]
      measures[n.ns, "index"] = peak.i
      measures[n.ns, "peak.val"] = peak.val
      

      if (plot) {
        plot(dat, xaxt='n', yaxt='n', ylim=c(0,60),
             bty='n', type='l',xlab='', ylab='')
        ##abline(v=sur+1)
        max.time <- ns$ns.T * sur
        axis(1, at=c(0,1,2)*sur,
             ##labels=c('-300 ms', '0 ms', '+300 ms'))
             labels=c(-max.time, 0, max.time))

      }
      
      hm = find.halfmax(dat, peak.n=sur+1, frac=0.5, plot=plot)
      measures[n.ns, "durn"] = hm$durn* ns$ns.T
      if (plot) {
        text <- sprintf("%d durn %.3f",
                        round(peak.val), measures[n.ns, "durn"])
        legend("topleft", text, bty='n')
      }

      ##dat2 = dat;
      ##dat2[1:(hm$xl-1)] = 0;
      ##dat2[(hm$xr+1):((2*sur)+1)] = 0;
      
      ##k = kurtosis(dat2)
      ##measures[n.ns, 1] = k
      ave = ave + dat


    }
  }

  if (n.ns < nrow(p)) {
    ## Some peaks could not be averaged, since they were at either
    ## beg/end of the recording.
    ## So, in this case, truncate the matrix of results to correct
    ## number of rows.
    measures = measures[1:n.ns,,drop=FALSE]
  }
  
  ## now show the average
  if (n.ns > 0) {
    ave = ave/n.ns
    if (plot) {
      plot(ave, xaxt='n', yaxt='n', bty='n', type='l',xlab='', ylab='')
      legend("topleft", paste("m", round(max(ave))), bty='n')
      find.halfmax(ave)
    }


    ##stripchart(measures[,1], ylab='K', method='jitter', vert=T, pch=19,
    ##main=paste('kurtosis', round(mean(measures[,1]),3)))
    if (plot) {
      stripchart(measures[,"durn"], ylab='durn (s)', method='jitter',
                 vert=TRUE, pch=19,
                 main=paste('FWHM durn', round(mean(measures[,"durn"]),3)))
    }

    if (plot) {
      par(old.par)
    }

  }

  
  ns.mean = ts(ave, start=(-sur*ns$ns.T), deltat=ns$ns.T)

  list(measures=measures, ns.mean=ns.mean)
}

show.ns <- function(ns, ...) {
  ## Show the individual network spikes after they have been computed.
  ##
  ## This is useful if you don't show the individual network spikes
  ## when they are first iterated over to calculate the mean.
 
  res <- mean.ns(ns, p=NULL, ...)
  NULL                                  #ignore result
}



find.peaks <- function(trace, ns.N) {

  max.peaks = 20000

  npts = length(trace)
  
  peaks = matrix(NA, nrow=max.peaks, ncol=2)
  colnames(peaks) <- c("index", "peak.val")
  n = 0

  inside = FALSE;

  for (i in 1:npts) {

    cur = trace[i]    

    if(inside) {
      ## currently in a peak.
      if (cur == 0) {
        ## no longer inside a peak, save results if peak was tall enough.
        inside=FALSE;

        if (peak > ns.N) {
          n=n+1
          if (n > max.peaks) {
            ## oh oh, need more room.
            browser()
          } else {
            peaks[n,] = c(peak.t, peak)
          }
        }
        
      } else {
        ## still in a peak
        if (cur > peak) {
          peak = cur; peak.t = i;
        }
      }
    } else {
      ## currently outside a peak.
      if (cur > 0) {
        inside = TRUE; peak = cur; peak.t = i
      }
    }
  }

  ## tidy up result at end.

  if (n > 0) {
    peaks = peaks[1:n,,drop=FALSE]
  } else {
    ## no peaks found.
    peaks = NULL
  }
}




find.halfmax.cute <- function(y) {
  ## Given a peak within DAT, find the FWHM.
  ## This is a cute method, but not robust enough -- since it assumes
  ## that the peak is unimodel -- which may not be the case.
  

  x = 1:length(y)
  p.t = 101                             #HACK!!!
  half.max = y[p.t]/2                   #HACK!!!
  f <- approxfun(x, y)
  f2 <- function(x) { f(x) - half.max }
  l <- uniroot(f2, lower=1, upper=p.t)
  r <- uniroot(f2, lower=p.t, upper=length(y))

  segments(l$root, f(l$root), r$root, f(r$root), col='blue')

}



find.halfmax <- function(y, peak.n=NULL, plot=TRUE, frac=0.5) {

  ## Given a peak somwhere with Y, find the FWHM.
  ##
  ## If PEAK.N is not null, it will be location of the peak -- this is helpful
  ## when there are multiple peaks within one window, and we want to find
  ## the FWHM of the non-major peak.
  ## By default, frac = 0.5, to find the half max.  Change this to some other
  ## value, e.g. 10% to find 10% onset and offset.
  ## 
  ##
  ## This may fail for a few reasons, e.g. not finding half-max values within
  ## the range, running out of data...
  ## all of which should be counted for!

  n = length(y)

  if (is.null(peak.n))
    peak.n = which.max(y)
  
  peak.val = y[peak.n]

  half.max = peak.val * frac
  
  ## Break the data into three segments:

  ## llllllllllllllllllPrrrrrrrrrrrrrrrrr
  ## P is the peak; examine curve to the left (lll) and to the right (rrr) to
  ## find when the peak has decayed to half max.

  left.y = y[1:(peak.n-1)]
  right.y = y[(peak.n+1):n]

  underhalf.l = which(left.y < half.max)
  xl1 = underhalf.l[length(underhalf.l)]   #get last point under halfmax.
  xl2 = xl1+1

  yl1 = y[xl1]; yl2 = y[xl2];
  dy = half.max - yl1

  ## below, (xl2 - xl1) should equal 1.
  dx = (dy  *(xl2-xl1))/ (yl2-yl1)

  ##xl.half = xl1 + dx
  xl.half = xl1 + dx

  ## Now work on right of curve.  find first point at which y falls below
  ## half max value.
  underhalf.r = which(right.y < half.max)
  xr2 = underhalf.r[1] + peak.n
  xr1 = xr2 - 1

  yr1 = y[xr1]; yr2 = y[xr2]
  dy = half.max - yr2

  dx = dy * (xr1-xr2)/(yr1-yr2)
  stopifnot(dx<0)
  xr.half = xr2 + dx

  
  if(plot) {
    ##abline(v=xl.half, col='green'); abline(v=xr.half, col='green'); #temp
    abline(h=peak.val * frac, col='red')
    segments(xl.half, half.max, xr.half, half.max, col='blue')
  }

  list(xl=xl.half, xr=xr.half, durn=xr.half-xl.half)
}

## now interpolate -- hard way
## a <- approx(x, y, n=length(y)*30)
## lines(a)

## amax.x = which.max(a$y)
## points(a$x[amax.x], a$y[amax.x], col='blue', pch=19)

## ## find right side down.
## half.max = max(y)/2
## rx <- which(a$y[-(1:amax.x)]< half.max)[1] + amax.x

## ## find left side down.
## lx <- which(a$y[1:amax.x]< half.max)
## lx <- lx[length(lx)]
## segments(a$x[lx], a$y[lx],  a$x[rx], a$y[rx], col='blue')

## The "R" way of interpolating -- nice!


check.ns.plot <- function(counts, p, xlim, ns.N) {

  plot(counts$times, counts$sum, type='l', xlim=xlim,
       xlab="time (s)", ylab='num active channels')
  points(counts$times[p[,1]], p[,2], pch=19, col='blue')
  abline(h=ns.N, col='red')               #threshold line.
}

ns.bin.peak <- function(p, nbins=12, wid=5) {
  ## Bin values in P into a set of NBINS bins, of size WID.
  ## Bins are right-closed (except for first bin, closed at both ends).
  ## Labels are added onto the bins.
  ##
  ## x <- c(0, 4,5, 20, 54,55, 60)
  ## ns.bin.peak(x, wid=10, nbins=7 )
  ##

  if ( is.null(p) ) {
    ## no valid values, so no need to make the histogram.
    ## This happens when there are no network spikes.
    p <- 0; invalid <- TRUE
  } else {
    invalid <- FALSE
  }
  
  b <- seq(from=0, by=wid, length=nbins+1)
  max.allowed <- max(b)
  if ( any( above <- which(p > max.allowed)) ) {
    stop("some values above max.allowed")
  }
  h <- hist(p, plot=F, breaks=b)
  c <- h$counts

  if (invalid) {
    ## no valid counts, so set all counts to zero.
    c <- c*0
  }
  
  l <- hist.make.labels(0, max.allowed, nbins)
  names(c) <- l
  c
}




######################################################################
## End of functions
######################################################################

kurtosis <- function (x, na.rm = FALSE) {
  ## Copied from e1071 package.
  if (na.rm) 
    x <- x[!is.na(x)]
  sum((x - mean(x))^4)/(length(x) * var(x)^2) - 3
}


spikes.to.count.old <- function(spikes,
                            time.interval=1, #time bin of 1sec.
                            beg=floor(min(unlist(spikes))),
                            end=ceiling(max(unlist(spikes)))
                            )
{
  ## First version: do not use!
  ##  The C version below is much faster:
  ##
  ## > unix.time(counts2 <- spikes.to.count(s$spikes, time.interval=ns.T))
  ##  user  system elapsed 
  ## 0.027   0.022   0.049 
  ## > unix.time(counts <- spikes.to.count.old(s$spikes, time.interval=ns.T))
  ## user  system elapsed 
  ## 11.826   7.780  19.613 
  ##
  ## Convert the spikes for each cell into a firing rate (in Hz)
  ## We count the number of spikes within time bins of duration
  ## time.interval (measured in seconds).
  ##
  ## Currently cannot specify BEG or END as less than the
  ## range of spike times else you get an error from hist().  The
  ## default anyway is to do all the spikes within a data file.
  ##
  ## This has adapted from make.spikes.to.frate

  
  spikes.to.counts <- function(spikes, breaks, time.interval) {
    ## helper function.  Convert one spike train into a count of how many
    ## spikes are within each bin.
    h <- hist(spikes, breaks=breaks,plot=FALSE,
              ## right=F makes closer agreement (but not total) with C version.
              right=F)
    res = h$counts

    ## We may want to check presence/absence of spike within a bin.
    multi.spikes = which(res>1)
    if (any(multi.spikes)) {
      res[multi.spikes] = 1
    }
    
    res
  }
  
  time.breaks <- seq(from=beg, to=end, by=time.interval)
  if (time.breaks[length(time.breaks)] < end) {
    ## extra time bin needs adding.
    ## e.g seq(1,6, by = 3) == 1 4, so we need to add 7 ourselves.
    time.breaks <- c(time.breaks,
                     time.breaks[length(time.breaks)]+time.interval)
  }

  very.bad = FALSE
  if (very.bad) {
    counts1 <- lapply(spikes, spikes.to.counts, breaks=time.breaks,
                      time.interval=time.interval)
    
    ## counts1 is a list; we want to convert it into an array.
    counts <- array(unlist(counts1),
                    dim=c(length(time.breaks)-1, length(counts1)))
    ## Do the average computation here.
    ## av.rate == average rate across the array.
    sum.rate <- apply(counts, 1, sum)
  } else {
    ncells <- s$NCells
    for (i in 1:ncells) {
      c = spikes.to.counts(s$spikes[[i]], time.breaks, time.interval)
      if (i == 1) {
        sum.rate = c
      } else {
        sum.rate = sum.rate + c
      }
    }

  }

  
  ## We can remove the last "time.break" since it does not correspond
  ## to the start of a time frame.
  res <- list(
              ##counts=counts,
              times=time.breaks[-length(time.breaks)],
              sum=sum.rate,
              time.interval=time.interval)
  res
}

show.ns.old <- function(p, counts, nrow=8, ncol=8, ask=FALSE, plot=TRUE) {
  ## Show the network spikes.

  ## This code does not check to worry if there is a spike right at either
  ## end of the recording.  naughty!

  if (is.null(p)) {
    cat("*** No network spikes found\n")
    return (NULL)
  }
  if (plot) {
    old.par <- par(mfrow=c(nrow,ncol), mar=c(2.5,1,1,1),ask=ask)
  }
  
  sur = 100                              #100 normally
  ave = rep(0, (2*sur)+1)
  npts = length(counts$sum)
  measures = matrix(NA, nrow=nrow(p), ncol=3)
  colnames(measures) = c("index", "peak.val", "durn")
  n.ns = 0                              #Number of valid network spikes found
  for (i in 1:nrow(p)) {
    peak.i = p[i,1]; lo = (peak.i-sur); hi = peak.i+sur

    ## Check that enough data can be found:
    if ( (lo >0) && ( hi < npts) ) {
      n.ns = n.ns + 1

      dat = counts$sum[lo:hi]
      peak.val = dat[sur+1]

      measures[n.ns, 1] = peak.i
      measures[n.ns, 2] = peak.val

      if (plot) {
        plot(dat, xaxt='n', yaxt='n', ylim=c(0,60),
             bty='n', type='l',xlab='', ylab='')
        ##abline(v=sur+1)
        axis(1, at=c(0,1,2)*sur,
             labels=c('-300 ms', '0 ms', '+300 ms'))
        legend("topleft", paste(round(peak.val)), bty='n')
      }
      

      hm = find.halfmax(dat, peak.n=sur+1, frac=0.5, plot=plot)
      measures[n.ns, 3] = hm$durn

      ##dat2 = dat;
      ##dat2[1:(hm$xl-1)] = 0;
      ##dat2[(hm$xr+1):((2*sur)+1)] = 0;
      
      ##k = kurtosis(dat2)
      ##measures[n.ns, 1] = k
      ave = ave + dat


    }
  }

  if (n.ns < nrow(p)) {
    ## Some peaks could not be averaged, since they were at either
    ## beg/end of the recording.
    ## So, in this case, truncate the matrix of results to correct
    ## number of rows.
    measures = measures[1:n.ns,,drop=FALSE]
  }
  
  ## now show the average
  if (n.ns > 0) {
    ave = ave/n.ns
    if (plot) {
      plot(ave, xaxt='n', yaxt='n', bty='n', type='l',xlab='', ylab='')
      legend("topleft", paste("m", round(max(ave))), bty='n')
      find.halfmax(ave)
    }

    if (plot) {
      par(old.par)
    }

    ##stripchart(measures[,1], ylab='K', method='jitter', vert=T, pch=19,
    ##main=paste('kurtosis', round(mean(measures[,1]),3)))
    if (plot) {
      stripchart(measures[,2], ylab='durn (bins)', method='jitter',
                 vert=TRUE, pch=19,
                 main=paste('FWHM durn', round(mean(measures[,2]),3)))
    }
  }


  list(measures=measures, ns.mean=ave)
}
