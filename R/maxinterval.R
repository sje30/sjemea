## Implement maxinterval method for burst detection.
## Fri 23 Feb 2007

mi.find.bursts <- function(spikes,debug=FALSE) {

  ## For one spike train, find the burst using max interval method.
  ## e.g.
  ## find.bursts(s$spikes[[5]])
  ## init.
  ## params currently in MI.PAR

  par = mi.par
  beg.isi = mi.par$beg.isi
  end.isi = mi.par$end.isi
  min.ibi = mi.par$min.ibi

  min.durn = mi.par$min.durn
  min.spikes = mi.par$min.spikes
  
  nspikes = length(spikes)

  n = 1

  ## Create a temp array for the storage of the bursts.  Assume that it
  ## will not be longer than Nspikes/2.
  max.bursts <- floor(nspikes/2)
  bursts <- matrix(NA, nrow=max.bursts, ncol=3)
  colnames(bursts) = c("beg", "end", "IBI")
  burst <- 0

  ## Phase 1 -- detection.
  ## Find ISIs closer than beg.isi, and end with end.isi.

  n = 2
  nspikes = length(spikes)
  in.burst = FALSE

  ## LAST.END is the time of the last spike in the previous burst.
  ## This is used to calculate the IBI.
  ## For the first burst, this is no previous IBI

  last.end = NA;                        #for first burst, there is no IBI.
  while ( n < nspikes) {
    
    next.isi = spikes[n] - spikes[n-1]
    if (in.burst) {
      if (next.isi > end.isi) {
        ## end of burst
        end = n-1; in.burst = FALSE

        
        ibi =  spikes[beg] - last.end; last.end = spikes[end]
        res = c(beg, end, ibi)
        ##printf("found burst %d %d %.3f \n", beg, end, ibi)
        burst = burst + 1
        if (burst > max.bursts) {
          print("too many bursts!!!")
          browser()
        }
        bursts[burst,] <- res
      }
    } else {
      ## not yet in burst.
      if (next.isi < beg.isi) {
        beg = n-1; in.burst = TRUE
      }
    }
    n = n+1
  }

  ## At the end of the burst, check if we were in a burst.
  if (in.burst) {
    end = nspikes; in.burst = FALSE
    ibi =  spikes[beg] - last.end; last.end = spikes[end]
    res = c(beg, end, ibi)
    ##printf("found burst %d %d %.3f \n", beg, end, ibi)
    burst = burst + 1
    if (burst > max.bursts) {
      print("too many bursts!!!")
      browser()
    }
    bursts[burst,] <- res
  }

  

  if (burst > 0 ) {
    ## truncate to right length.
    bursts = bursts[1:burst,,drop=F]
  }
  if (debug) {
    print("End of phase1\n")
    print(bursts)
  }
  
  
  ## Phase 2 -- merging.
  if (burst > 0 ) {
    ibis = bursts[,3]
    ##print(ibis)
    merge.bursts = which(ibis < min.ibi)
  
    if (any(merge.bursts)) {
      ## Remove bursts one by one.
      ## This might be inefficient.
      offset = 0
      ##browser()
      for (burst in merge.bursts) {

        ##print(bursts)
        burst = burst-offset; offset=offset+1
        ## update "end" information.
        bursts[burst-1, "end"] = bursts[burst, "end"]
        bursts = bursts[-burst,,drop=F] #delete this row.

        ##bursts[burst,"end"] = NA #signal that is corrupt, but removed later.

      }
    }
  }
  if (debug) {
    print("End of phase 2\n")
    print(bursts)
  }


  ## phase 3 -- remove small bursts: those less than min duration, or
  ## min number of spikes.

  len = bursts[,"end"] - bursts[,"beg"] + 1
  durn = spikes[bursts[,"end"]] - spikes[bursts[,"beg"]]
  bursts = cbind(bursts, len, durn)

  rejects = which ( (durn < min.durn) | ( len < min.spikes) )
  
  if (any(rejects)) {
    bursts = bursts[-rejects,,drop=FALSE]
  }


  ## Compute mean ISIS
  len = bursts[,"end"] - bursts[,"beg"] + 1
  durn = spikes[bursts[,"end"]] - spikes[bursts[,"beg"]]
  mean.isis = durn/(len-1)                  #TODO -- CHECK.

  SI = rep(1, length(mean.isis ))
  bursts = cbind(bursts, mean.isis, SI)
  ## If any bursts were rejected, we now need to recalc IBI.
    
  ## End
  bursts
  
}

