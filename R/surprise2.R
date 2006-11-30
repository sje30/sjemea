s.min = 5


plot.par <- function(allb, ylab, index, max=-1,title='') {
  ## plot result over channels, e.g. burst duration, or si.
  plot.channels <- min(length(allb), 20)
  durns <- list()
  for (i in 1:plot.channels) {
    ##print(i)
    b <- allb[[i]]
    if (length(b)==1) {
      res <- NA
    } else {
      if (length(b) ==4) {
        ## have only one burst!  TO CHECK...! (e.g. spike train 18)
        res <- b[index]
      } else {
        res <- b[,index]
      }
    }

    ## TODO -- strip out any INF values that creep into SI.
    infs <- which(res==Inf)
    ##print(infs)
    ##print(res)
    if (length(infs)>0)
      res <- res[-infs]
    
    durns[[i]] <- res
  }

  if (max>0) {
    durns <- sapply(durns, pmin, max)
    ##browser()
  }

  mins <- min(sapply(durns, min), na.rm=TRUE)
  maxs <- max(sapply(durns, max), na.rm=TRUE)
  stripchart(durns, method="jitter", pch=20, vert=T,main=title,
             ylim=c(mins,maxs),
             xlab='channel', ylab=ylab)

  ## return value.
  ##durns

}

spikes.to.bursts.surprise <- function(s) {

  ncells <- s$NCells

  ncells <- 20                           #temp
  allb <- list()
  for (train in 1:ncells) {
    cat(sprintf("** analyse train %d\n", train))
    spikes <- s$spikes[[train]]

    bursts <- find.bursts(spikes)

    ## SJE: remove any NAs from bursts that occur -- check later why
    ## they occur.
    if (length(bursts)>1) {             #i.e. not NA returned.
      ## we found some bursts.
      bad.rows <- which (is.na(bursts[,1]))
      ##browser()
      if (length(bad.rows)>0)
        bursts <- bursts[-bad.rows,]
    }
    allb[[train]] <- bursts
  }

  allb
}

find.bursts <- function(spikes,debug=FALSE) {

  ## e.g.
  ## find.bursts(s$spikes[[5]])
  ## init.
  nspikes = length(spikes)
  mean.isi = nspikes/ (spikes[nspikes] - spikes[1])
  threshold = mean.isi/2

  n = 1

  ## Create a temp array for the storage of the bursts.  Assume that it
  ## will not be longer than Nspikes/3.
  max.bursts <- floor(nspikes/3)
  bursts <- matrix(NA, nrow=max.bursts, ncol=4)
  burst <- 1

  ## note, no need to check to end of spike train!
  while ( n < nspikes-2) {              #end condition to check!!!
    if (debug)
      print(n)
    if( ((spikes[n+1] - spikes[n]  ) < threshold) &&
       ((spikes[n+2] - spikes[n+1]) < threshold)) {
      res <- find.burst(n, spikes, nspikes, mean.isi, threshold,debug)

      if (is.na(res[1])) {
        ## no burst found, just move on one spike.
        n <- n + 1
      } else {
        bursts[burst,] <- res
        burst <- burst + 1
        if (burst > max.bursts) {
          print("too many bursts")
          browser()
        }
        n <- n + res[2]                 #move to end of burst.
      }
    } else {
      ## no triple-spike.
      n = n + 1
    }
  }

  ## At end of spike train, now truncate bursts to right length.
  if (burst > 1) {
    res <- bursts[1:burst,]
    if(debug)
      print(dim(res))
    colnames(res) <- c("first", "len", "SI", "durn")
  } else {
    res <- NA
  }

  res
  
}

find.burst <- function(n, spikes, nspikes, mean.isi, threshold,debug) {
  ## Find a burst starting at spike N.
  ## TODO -- need to worry about running out of spikes at the end of the
  ## spike train!

  if (debug) 
    cat(sprintf("** find.burst %d\n", n))
  i=3  ## First three spikes are in burst.
  s = surprise(n, i, spikes, nspikes, mean.isi)
  if (s > s.min) {
    if (debug)
      cat(sprintf("starting phase 1 n %d s %.4f\n", n, s))
    ## Phase 1 - add spikes to the train.
    phase1 = TRUE
    while( phase1) {
      s.new = surprise(n, i+1, spikes, nspikes, mean.isi)
      ## TODO -- useful debug info.
      ## cat(sprintf("phase 1: n %d i %d s.new %.4f\n", n, i, s.new))
      if (s.new > s) {
        s = s.new
        i = i+1
      } else {
        phase1 = FALSE
      }
    }
    ## start deleting spikes from the start of the burst.
    phase2 = TRUE
    while(phase2) {
      s.new = surprise(n, i, spikes, nspikes, mean.isi)
      if (debug)
        cat(sprintf("phase 2: n %d i %d s.new %.4f\n", n, i, s.new))        
      if (s.new > s) {
        n = n+1; i = i-1
        s = s.new
        if (i==2) {
          ## perhaps set end of phase2 here in this case?
          print("i ==2 in phase 2")
          browser()
        }
      } else {
        phase2 = FALSE
      }
    }
    durn = spikes[n+i] - spikes[n]
    res <- c(n=n, i=i, s=s, durn=durn)
  } else {
    res <- rep(NA, 4)
  }
  ##browser()
  res
  
}

surprise <- function(n, i, spikes, nspikes, mean.isi) {
  ## Calculate surprise.

  if (n+i> nspikes) {
    ## Handle case when we are at the end of a spike train...
    ## Needs to be done in several situations. e.g. find.bursts
    s = 0
  }
  else {
    dur <- spikes[n+i] - spikes[n]
    lambda <- dur / mean.isi
    p <- ppois(i-2, lambda, lower.tail=FALSE)
    s = -log(p)
  }

  s
}


plot.spikes <- function(xlim=NULL, show.bursts=TRUE) {
  min.t <- spikes[1]
  max.t <- spikes[nspikes]

  if (is.null(xlim)) {
    xlim <- c(min.t, max.t)
  } 
  plot(NA, xlim=xlim, xlab='time', ylim=c(0,1))
  segments(spikes, rep(0.2, nspikes), spikes, rep(0.8, nspikes))

  
  if (show.bursts) {
    burst.x1 <- spikes[b[,1]]
    burst.x2 <- spikes[b[,1]+ b[,2]]
    burst.y <- rep(0.5, length=length(burst.x1))
    segments(burst.x1, burst.y, burst.x2, burst.y, col='red', lwd=3)
  }
}
       
bursts.to.active <- function(bursts, tmin, tmax, dt) {

  nbins = floor((tmax-tmin)/dt)+1

  active = vector(length=nbins)         #default all FALSE.

  nbursts = nrow(bursts)

  for (b in 1:nbursts) {
    burst.start =  spikes[ bursts[b,1]]
    burst.stop =  spikes[ bursts[b,1] + bursts[b,2]]

    cat(sprintf("burst %d from %f to %f\n", b, burst.start, burst.stop))

    start.bin = floor( (burst.start - tmin)/dt) + 1
    stop.bin =  floor( (burst.stop  - tmin)/dt) + 1
    bins = start.bin:stop.bin
    for (bin in bins)
      active[bin] = TRUE
  }
  names(active) <- seq(from=tmin, by=dt, length=nbins)
  active
}
