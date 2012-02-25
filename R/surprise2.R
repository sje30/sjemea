## Poisson surprise method for burst analysis -- L'egendy and Salcman (1985)
## Author: Stephen J Eglen
## Copyright: GPL

s.min = 5                               #threshold on suprise index.

burst.isi.threshold = FALSE             #do we want to use threshold on ISI?

##burst.isi.max = 0.1 #ISI within burst must be smaller than this.

burst.isi.max = NULL                    #set non-null to be the threshold between spikes.

## Threshold for burstiness; how many burst/minute are there to count as a
## bursty unit.
bursty.threshold = 1

######################################################################

burst.info <- c("beg", "len", "SI", "durn", "mean.isis")
burst.info.len = length(burst.info)



spikes.to.bursts <- function(s, method="si") {
  ## Entry function for burst analysis.
  ## Possible methods:
  ## "mi" - maxinterval
  ## "si" - surprise index.
  ## "logisi" - log of the ISI histogram
  
  ncells <- s$NCells

  if (method == "logisi") {
    isi.low <- logisi.compute(s)$Locmin
    logisi.par$isi.low <- isi.low
  }
  
  ##ncells <- 10                           #temp
  allb <- list()
  for (train in 1:ncells) {
    ## cat(sprintf("** analyse train %d\n", train))
    spikes <- s$spikes[[train]]

    bursts = switch(method,
      "mi" = mi.find.bursts(spikes),
      "si" = si.find.bursts(spikes),
      "logisi" = logisi.find.burst(spikes),
      "rs" = rs.find.bursts(spikes),
      stop(method, " : no such method for burst analysis")
    )

    allb[[train]] <- bursts
  }

  allb
}


spikes.to.bursts.surprise <- function(s) {
  ## Wrapper function for surprise method.
  spikes.to.bursts(s, method="si")
}



si.find.bursts <- function(spikes,debug=FALSE) {
  ## For one spike train, find the burst using SI method.
  ## e.g.
  ## find.bursts(s$spikes[[5]])
  ## init.
  nspikes = length(spikes)

  ### OOOOPS!!!!  mean.isi was actually the mean firing rate!
  mean.isi = mean(diff(spikes))
  
  ##mean.isi = nspikes/ (spikes[nspikes] - spikes[1])
  threshold = mean.isi/2

  n = 1

  ## Create a temp array for the storage of the bursts.  Assume that it
  ## will not be longer than Nspikes/3.
  max.bursts <- floor(nspikes/3)
  bursts <- matrix(NA, nrow=max.bursts, ncol=burst.info.len)
  burst <- 0

  ## note, no need to check to end of spike train!
  while ( n < nspikes-2) {              #end condition to check!!!
    if (debug)
      print(n)
    if( ((spikes[n+1] - spikes[n]  ) < threshold) &&
       ((spikes[n+2] - spikes[n+1]) < threshold)) {
      ##res <- find.burst(n, spikes, nspikes, mean.isi, threshold,debug)
      res <- si.find.burst2(n, spikes, nspikes, mean.isi, burst.isi.max,debug)

      if (is.na(res[1])) {
        ## no burst found, just move on one spike.
        n <- n + 1
      } else {
        ## found a burst.
        burst <- burst + 1
        if (burst > max.bursts) {
          print("too many bursts")
          browser()
        }
        bursts[burst,] <- res
        ## WRONG!!! This assumes no phase 2!!!
        ## n <- n + res[2]                 #move to end of burst.
        ## move on to spike after n.
        ##n <- res["n"] + res["i"]
        ## TODO: why does names(res) change?
        ## Could try "0 + res[1] + res[2]"???
        n <- res[1] + res[2]
        names(n) <- NULL                #TODO HHHHH???
        ##browser()
      }
    } else {
      ## no triple-spike.
      n = n + 1
    }
  }

  ## At end of spike train, now truncate bursts to right length.

  if (burst > 0) {
    res <- bursts[1:burst,,drop=FALSE]
    colnames(res) <- burst.info
  } else {
    res <- NA
  }

  res
  
}


si.find.burst2 <- function(n, spikes, nspikes, mean.isi, threshold=NULL,
                        debug=FALSE) {
  ## Find a burst starting at spike N.
  ## Include a better phase 1.


  ## Determine ISI threshold.
  if (is.null(threshold)) 
    isi.thresh = 2 * mean.isi
  else
    isi.thresh = threshold
  
  if (debug) 
    cat(sprintf("** find.burst %d\n", n))
  
  i=3  ## First three spikes are in burst.
  s = surprise(n, i, spikes, nspikes, mean.isi)

  ## Phase 1 - add spikes to the train.
  phase1 = TRUE
  ##browser()

  ## in Phase1, check that we still have spikes to add to the train.
  while( phase1 ) {

    ##printf("phase 1 s %f\n", s);
    
    i.cur = i;

    ## CHECK controls how many spikes we can look ahead until SI is maximised.
    ## This is normally 10, but will be less at the end of the train.
    check = min(10, nspikes-(i+n-1))

    looking = TRUE; okay = FALSE;
    while (looking) {

      if (check==0) {
        ## no more spikes left to check.
        looking=FALSE;
        break;
      }
      check=check-1; i=i+1
      s.new = surprise(n, i, spikes, nspikes, mean.isi)
      if (debug) 
        printf("s.new %f s %f n %d i %d check %d\n", s.new, s, n, i, check)

      if (s.new > s) {
        okay=TRUE; looking=FALSE;
      } else {
        ## See if we should keep adding spikes?
        if ( (spikes[i] - spikes[i-1]) > isi.thresh ) {
          looking = FALSE;
        }
          
      }
    }
    ## No longer checking, see if we found an improvement.
    if (okay) {
      if (s > s.new) {
        ## This should not happen.
        printf("before s %f s.new %f\n", s, s.new)
        browser()
      }
      s = s.new
    } else {
      ## Could not add more spikes onto the end of the train.
      phase1 = FALSE
      i = i.cur
    }
  }


  ## start deleting spikes from the start of the burst.
  phase2 = TRUE
  while(phase2) {
    if (i==3) {
      ## minimum length of a burst must be 3.
      phase2=FALSE
    } else {
      s.new = surprise(n+1, i-1, spikes, nspikes, mean.isi)
      if (debug)
        cat(sprintf("phase 2: n %d i %d s.new %.4f\n", n, i, s.new))        
      if (s.new > s) {
        if (debug) 
          print("in phase 2 acceptance\n")
        n = n+1; i = i-1
        s = s.new
      } else {
        ## removing front spike did not improve SI.
        phase2 = FALSE
      }
    }
  }
  

  ## End of burst detection; accumulate result.
  if ( s > s.min) {


    ## compute the ISIs, and then the mean ISI.
    
    ## Fencepost issue: I is the number of spikes in the burst, so if
    ## the first spike is N, the last spike is at N+I-1, not N+I.
    isis = diff(spikes[n+(0:(i-1))])
    mean.isis = mean(isis)
    
    durn = spikes[n+i-1] - spikes[n]
    res <- c(n=n, i=i, s=s, durn=durn, mean.isis=mean.isis)

    if (debug) 
      print(res)

  } else {
    ## burst did not have high enough SI.
    res <- rep(NA, burst.info.len)
  }
  ##browser()
  res
  
}

surprise <- function(n, i, spikes, nspikes, mean.isi) {
  ## Calculate surprise index for spike train.

  ##stopifnot(n+i <= nspikes)
  dur <- spikes[n+i-1] - spikes[n]
  lambda <- dur / mean.isi
  p <- ppois(i-2, lambda, lower.tail=FALSE)
  s = -log(p)

  s
}


       
######################################################################
## General methods, not just for Surprise Index Method.

calc.burst.summary <- function(s) {
  ## Compute the summary burst information.  Use a separate 
  ## call (write.csv() for example) to write the burst information to file.
  
  allb <- s$allb
  
  ## Create a table of output results.

  channels <- s$channels
  spikes <- as.vector(s$nspikes)

  duration <- s$rec.time[2]  - s$rec.time[1]

  mean.freq <- round(spikes/duration, 3)

  nbursts <- sapply(allb, num.bursts)

  bursts.per.sec <- round(nbursts/duration,3)
  bursts.per.min <- bursts.per.sec * 60


  bursty = ifelse(bursts.per.min >= bursty.threshold, 1, 0)

  durations <- burstinfo(allb, "durn")
  mean.dur <- round(sapply(durations, mean), 3)
  sd.dur <- round(sapply(durations, sd), 3)

  
  
  ##mean.isis <- burstinfo(allb, "mean.isis")
  ##mean.mean.isis <- round(sapply(mean.isis, mean), 3)
  ##sd.mean.isis <- round(sapply(mean.isis, sd), 3)
  ISIs = calc.all.isi(s, allb)
  mean.ISIs = sapply(ISIs, mean)
  sd.ISIs = unlist( sapply(ISIs, sd, na.rm=TRUE))

  
  ns <- burstinfo(allb, "len")
  mean.spikes <- round(sapply(ns, mean), 3)
  sd.spikes   <- round(sapply(ns, sd), 3)
  total.spikes.in.burst <- sapply(ns, sum)
  per.spikes.in.burst <- round(100 *(total.spikes.in.burst / spikes), 3)

  si <- burstinfo(allb, "SI")
  mean.si <- round(sapply(si, mean), 3)


  IBIs <- calc.all.ibi(s, allb)
  mean.IBIs <- sapply(IBIs, mean)
  sd.IBIs <- sapply(IBIs, sd, na.rm=TRUE)
  cv.IBIs <- round(sd.IBIs/ mean.IBIs, 3)
  ## round afterwards...
  mean.IBIs <- round(mean.IBIs, 3); sd.IBIs <- round(sd.IBIs, 3)
  
  df <- data.frame(channels=channels, spikes=spikes, mean.freq=mean.freq,
                   nbursts=nbursts,
                   bursts.per.sec=bursts.per.sec,
                   bursts.per.min=bursts.per.min,
                   bursty = bursty,
                   mean.dur=mean.dur,
                   sd.dur=sd.dur,
                   mean.spikes=mean.spikes,
                   sd.spikes=sd.spikes,
                   per.spikes.in.burst=per.spikes.in.burst,
                   per.spikes.out.burst=round(100.0-per.spikes.in.burst,3),
                   mean.si=mean.si,
                   mean2.isis=mean.ISIs,
                   sd.mean.isis=sd.ISIs,
                   mean.IBIs=mean.IBIs,
                   sd.IBIs=sd.IBIs,
                   cv.IBIs=cv.IBIs
                   )
  ##write.csv(df, file=outfile)

  df

}

mean.burst.summary = function(allb.sum) {
  ## Summarise the burst information.  This does not handle per.spikes.in.burst
  subset = allb.sum[which(allb.sum$bursty==1),]
  
  fields = c("spikes", "mean.dur", "cv.IBIs", "bursts.per.min", "per.spikes.in.burst")
  res = rep(0, length(fields)*2)
  names(res) = paste(rep(fields, each=2), c("m", "sd"), sep=".")
  n = 1
  for (field in fields) {
    dat = subset[[field]]
    if (length(dat) > 0 ) {
      mean = mean(dat, na.rm=TRUE); sd = sd(dat, na.rm=TRUE);
    } else {
      mean = sd = NA;
    }
    res[n] = mean; res[n+1] = sd
    n = n +2
  }

  res

}




burstinfo <- function(allb, index) {
  ## Extra some part of the Burst information, for each channel.
  ## index will be the name of one of the columns of burst info.
  ## This is a HELPER function for calc.burst.summary
  sapply(allb, function(b) {
    if (length(b)>1) {
      b[,index]
    } else {
      0
    }
  }, simplify=FALSE)
}
  



calc.ibi <- function(spikes, b) {
  ## Compute the interburst intervals (IBI) for one spike train.
  ## Only valid if more than one burst.

  nburst = num.bursts(b)
  if ( nburst == 0) {
    res = NA                            #no bursts
  } else {
    if (nburst == 1) {
      res = NA                          #cannot compute  IBI w/only 1 burst.
    } else {
      ## find end spike in each burst.
      end = b[,"beg"] + b[,"len"] - 1

      ## for NBURST bursts, there will be NBURST-1 IBIs.
      start.spikes = b[2:nburst,"beg"]
      end.spikes   = end[1:(nburst-1)]
      ## NEX uses a strange definition of IBI -- it counts the time between
      ## the first spike of burst N and the first spike of burst N+1 as the
      ## IBI.  If we want to use that definition, use the following line:
      ##end.spikes   = b[1:(nburst-1),"beg"]
      res = spikes[start.spikes] - spikes[end.spikes]
    }
  }
  res
}

calc.all.ibi <- function (s, allb) {
  ## Compute IBI for all spike trains.
  nchannels <- s$NCells
  IBIs = list()
  for (i in 1:nchannels) {
    IBIs[[i]]  = calc.ibi(s$spikes[[i]], allb[[i]])
  }

  IBIs
}


calc.all.isi <- function (s, allb) {
  ## Compute ISI within bursts for all spike trains.

  calc.isi = function(spikes, b) {
    ## for one spike train, get all ISIs within bursts in that train.
    if (num.bursts(b)==0) {
      return ( NA )
    }

    ## as.vector is needed below in case each burst is of the same
    ## length (in which case an array is returned by apply).  In
    ## normal cases, bursts are of different lengths, so "apply"
    ## returns a list.
    
    isis = as.vector(
      unlist(apply(b, 1,
      function(x) {
        diff(spikes[ x[1]:x[2]])
      } )))
  }
  
  nchannels <- s$NCells
  ISIs = list()
  for (i in 1:nchannels) {
    ISIs[[i]]  = calc.isi(s$spikes[[i]], allb[[i]])
  }

  ISIs
}


num.bursts <- function(b) {
  ## Return the number of bursts found for a spike train.
  if(is.na(b[1]))
    0
  else
    nrow(b)
}

## Plotting code.


plot.burst.info <- function(allb, index, ylab=NULL, max=-1,title='') {
  ## Plot result of burst analysis in a stripchart, one col per channel.
  ## Feature plotted is given by index, e.g. "durn", "len".

  ##plot.channels <- min(length(allb), 70)
  plot.channels <- length(allb)         #plot all channels.
  
  values <- list()
  for (i in 1:plot.channels) {
    b <- allb[[i]]
    if (num.bursts(b)==0) {
      res <- NULL
    } else {
        res <- b[,index]
    }

    ## TODO -- strip out any INF values that creep into SI.
    infs <- which(res==Inf)
    ##print(infs)
    ##print(res)
    if (length(infs)>0)
      res <- res[-infs]
    
    values[[i]] <- res
  }

  if (max>0) {
    values <- sapply(values, pmin, max)
  }
  mins <- min(sapply(values, min), na.rm=TRUE)
  maxs <- max(sapply(values, max), na.rm=TRUE)

  if(is.null(ylab))
    ylab=index

  stripchart(values, method="jitter", pch=20, vert=TRUE,main=title,
             ylim=c(mins,maxs),
             xlab='channel', ylab=ylab)

  ## return value.
  ##values

}

bursts.to.active <- function(bursts, tmin, tmax, dt) {
  ## ??? Not sure if this is used right now.
  spikes <- NULL  ## TODO: fix spikes before it can be used.
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


######################################################################
## Old code below.

find.burst <- function(n, spikes, nspikes, mean.isi, threshold,debug) {
  ## Find a burst starting at spike N.

  if (debug) 
    cat(sprintf("** find.burst %d\n", n))
  i=3  ## First three spikes are in burst.
  s = surprise(n, i, spikes, nspikes, mean.isi)
  if (s > s.min) {
    if (debug)
      cat(sprintf("starting phase 1 n %d s %.4f\n", n, s))
    ## Phase 1 - add spikes to the train.
    phase1 = TRUE
    i.asn = n+i-1 #current end index of spike train.

    ## in Phase1, check that we still have spikes to add to the train.
    while( phase1 && (i.asn < nspikes) ) {

      if (!burst.isi.threshold ||
          (( spikes[i.asn+1] - spikes[i.asn]) < burst.isi.max))
        s.new = surprise(n, i+1, spikes, nspikes, mean.isi)
      else
        s.new = 0                       #handles case when ISI > threshold.
      ## TODO -- useful debug info.
      ## cat(sprintf("phase 1: n %d i %d s.new %.4f\n", n, i, s.new))
      if (s.new > s) {
        s = s.new
        i = i+1; i.asn = i.asn+1;
      } else {
        phase1 = FALSE
      }
    }
    ## start deleting spikes from the start of the burst.
    phase2 = TRUE
    while(phase2) {
      s.new = surprise(n+1, i-1, spikes, nspikes, mean.isi)
      if(is.na(s.new))
        browser()
      if (debug)
        cat(sprintf("phase 2: n %d i %d s.new %.4f\n", n, i, s.new))        
      if (s.new > s) {
        if (debug) 
          print("in phase 2 acceptance\n")
        n = n+1; i = i-1
        s = s.new
        if (i==2) {
          ## perhaps set end of phase2 here in this case?
          ## this will happen!!!! TODO
          print("i ==2 in phase 2")
          phase2=FALSE
          ##browser()
        }
      } else {
        phase2 = FALSE
      }
    }

    ## End of burst detection; accumulate result.

    ## compute the ISIs, and then the mean ISI.

    ## Fencepost issue: I is the number of spikes in the burst, so if
    ## the first spike is N, the last spike is at N+I-1, not N+I.
    isis = diff(spikes[n+(0:(i-1))])
    mean.isis = mean(isis)
    
    durn = spikes[n+i-1] - spikes[n]
    res <- c(n=n, i=i, s=s, durn=durn, mean.isis=mean.isis)
  } else {
    res <- rep(NA, burst.info.len)
  }
  ##browser()
  res
  
}

plot.spikes <- function(xlim=NULL, show.bursts=TRUE) {
  spikes <- b <- NULL  ## TODO: This requires SPIKES and B to be defined somwehwere...

  nspikes = length(spikes)
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
