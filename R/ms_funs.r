## Functions for reading multisite data (both Markus/Rachel's and Jay's)
## Mon 10 Sep 2001


## todo -- are the cell positions inverted?

## * Markus's functions.

mm.WrapTime <- 16** 4 * 128        #Clock wraparound in tics, ca 419s.
mm.sample.rate <- 20000.0               #20KHz sample rate.
mm.burst.sep <- 2


## Size of the various data types.
longsize <- 4; floatsize <- 4; intsize <- 2; unsize <- 2


mm.readpos <- function(posfile) {
  x <- read.table(posfile,sep="\t", skip=5, header=F)
  res <- cbind(x$V2, x$V3)
  if (dim(x)[2] == 4) rownames(res) <- as.character(x$V4)
  class(res) <- "mm.pos"
  res
}


mm.readpos.compare <- function(NCells, boxes, posfile) {
  ## Generate the multisite positions.
  ## Read in the position file if it was given to compare with my
  ## assignment of channels to electrode positions.
  guess.pos <- array(0, dim=c(NCells,2))
  ref.pos <- mm.readpos("ms_sje_pos.text")
  for (i in 1:NCells) {
    matches <- which(boxes[,1] == i)
    ##cat(paste("matches for cell",i,":",matches, "\n"))
    if (length(matches) == 0) {
      stop(paste("no matches found for cell",i))
    }
    channel <- boxes[matches[1],2]
    guess.pos[i,] <- ref.pos[channel,]
  }
  if (is.character(posfile)) {

    if (!file.exists(posfile))
      stop(paste("Position file",posfile,"does not exist"))
    
    ## now read in the spike position file.  Always return this
    ## if it is available.
    pos <- mm.readpos(posfile)
    if (NCells != dim(pos)[1]) {
      stop(paste("NCells against size in this.pos differs",
                 NCells, dim(pos)[1], "\n"))
    }

    ## compare the computed positions with those points read in from
    ## Markus' program.
    diffs <- pos - guess.pos
    dists <- apply(diffs, 1, function(x) { sqrt(sum(x**2))})
    if (any(dists)) {
      warning(paste("some cell positions wrong\n",
                    paste(which(dists >0),
                          signif(dists[which(dists>0)],4), "\n",
                          collapse=" "),"\n"))
    }
  } else {
    ## no data file, so just take guess.
    pos <- guess.pos
  }

  class(pos) <- "mm.pos"
  pos
}

plot.mm.pos <- function(x, use.rownames=F) {
  range <- c(-300,300)
  plot(x[,1], x[,2], xlim=range, ylim=range, xlab="", ylab="", type="n")
  if (use.rownames)
    text(x[,1], x[,2], rownames(x))
  else
    text(x[,1], x[,2])
}




read.ms.mm.data <- function(cellname, posfile=NULL) {
  ## Read in the multisite data and return a list with all the relevant
  ## data.  Determine the format of the file then call the appropriate
  ## routine.


  if(is.null(posfile) ) {
    posfile <- paste(cellname, ".pos", sep='')
    if (!file.exists(posfile))
      posfile <- NULL
    else
      cat(paste("guess posfile:",posfile, "\n"))
  }
    
  fp <- file(cellname , 'rb')
  Format <- readBin(fp, integer(), 1, longsize, endian="big")
  close(fp)

  if (Format == 2) 
    read.ms.mm.data.format2(cellname, posfile)
  else 
    read.ms.mm.data.format1(cellname, posfile)

}

read.ms.mm.data.format2 <- function(cellname, posfile=NULL) {
  ## Read in the multisite data and return a list with all the relevant
  ## data.

  ## Get the total size of the file so it can be compared with value
  ## of seek() once all the data have been read in.
  filesize <- file.info(cellname)$size
  
  fp <- file(cellname , 'rb')
  seek(fp,0)
  Format <- readBin(fp, integer(), 1, longsize, endian="big")
  t <- readBin(fp, integer(), 4, longsize, endian="big")
  FileIndex <- t[1]; BoxIndex <- t[2]; RecIndex <- t[3]; StatIndex <- t[4]

  ## Now read the NFiles...
  seek(fp, 64)
  t <- readBin(fp, integer(), 4, intsize, endian="big")
  NFiles <- t[1]; NBoxes <- t[2]; NRecords <- t[3]; NCells <- t[4]

  if (NFiles>1)
    warning(paste("NFiles larger than 1 - check ok? - e.g. endTimes",
                  NFiles, "\n"))
  
  t <- readBin(fp, integer(), 2, longsize, endian="big")
  NEvents <- t[1]; NSpikes <- t[2]
  cat(paste("NEvents", NEvents, "NSpikes", NSpikes, "\n"))

  ## Read in the fileinfo.
  if (seek(fp) != FileIndex)
    stop("error - current file position different from expected FileIndex")

  seek(fp, FileIndex)
  for (r in 1:NFiles) {
    vrn <- readBin(fp, integer(), 1, intsize, endian="big")
    pfilename <- readChar(fp, 64)
    pdirname <- readChar(fp, 64)
    flcrdat  <- readBin(fp, integer(), 1, longsize, endian="big")
    t <- readBin(fp, integer(), 3, intsize, endian="big")
    LowRecord <- t[1]; nrec <- t[2]; LastRecord <- t[3];
    cat(paste("File", r, "name", pfilename, "dir", pdirname,
              "nrec", nrec, "LastRecord", LastRecord, "\n"))
  }

  ## Read in the Box Info ####################
  if (seek(fp) != BoxIndex)
    stop("error - current file position different from BoxIndex")

  ## Make an array to store the boxes.
  boxes <- array(0, dim= c(NBoxes, 7))
  ## todo -- determine how these boxes relate to position of neurons.
  seek(fp, BoxIndex)
  for (r in 1:NBoxes) {
    t <- readBin(fp, integer(), 7, intsize, endian="big")
    group <- t[1]; channel <- t[2]; plott <- t[3]
    ##cat(paste("Box", r, "Group", group, "Chan", channel,
    ##"Plot", plott, "bounds", t[4], t[5], t[6], t[7], "\n"))
    boxes[r,] <- t
  }


  ## now parse the RecIndex... ####################

  if (seek(fp) != RecIndex)
    stop(paste("seek position different from expected RecIndex",
               seek(fp), RecIndex))

  ## RecIndex points to an array of length NRecords, each which points
  ## to the start of the rth record.

  seek(fp,RecIndex)
  RecordIndexes <- readBin(fp, integer(), NRecords, longsize, endian="big")

  ## Parse each record...

  startclock <- integer(NRecords)
  endclock   <- integer(NRecords)
  starttimes <- integer(NRecords)       #to be calculated...
  endtimes   <- integer(NRecords)
  nevents    <- integer(NRecords)
  nspikes    <- integer(NRecords)

  spikecount <- 0
  eventcount <- 0
  laststop   <- 0

  ## Make an empty list of size NCells.  Each element will be a list.
  allspikes <- list()
  for (i in 1:NCells) {
    allspikes[i] <- list()
  }
  EndTime <- 0;                        #todo: start of each file?

  for (r in 1:NRecords) {

    if ((laststop >0) && (seek(fp) != laststop)) {
      stop(paste("Error: RecordIndex and position of last byte differ",
                 "Record", r, "start", start, "laststop", laststop))
    }
    
    seek(fp, RecordIndexes[r])
    startclock[r] <- readBin(fp, integer(), 1, unsize, signed=F, endian="big")
    endclock[r]   <- readBin(fp, integer(), 1, unsize, signed=F, endian="big")
    nevents[r]    <- readBin(fp, integer(), 1, longsize, endian="big")
    nspikes[r]    <- readBin(fp, integer(), 1, longsize, endian="big")


    ShiftTime <- (EndTime %/% mm.WrapTime) * mm.WrapTime
    while ( ( (startclock[r] * 128) + ShiftTime) < EndTime)
      ShiftTime <- ShiftTime + mm.WrapTime

    starttimes[r] <- (startclock[r] * 128) + ShiftTime
    
    cat(paste(r, "clock", startclock[r], endclock[r], "#events", nevents[r],
              "#spikes", nspikes[r], "\n"))
    spikecount <- spikecount + nspikes[r]
    eventcount <- eventcount + nevents[r]
    ## Read in the number of spikes from each cell in record r.
    spikespercell <- readBin(fp, integer(), NCells, longsize, endian="big")

    ## Read in the time of event.
    eventsinrecord <- readBin(fp, integer(), nevents[r],longsize, endian="big")

    ## width of events
    we <- readBin(fp, integer(), nevents[r], intsize, endian="big")

    ## peak of events
    pe <- readBin(fp, integer(), nevents[r], intsize, endian="big")

    ## time of each spike from each cell in record r
    TLast <- -1
    for (cell in 1:NCells) {
      nspikescell <- spikespercell[cell]
      spiketimes <- readBin(fp, integer(), nspikescell, longsize, endian="big")


      if (nspikescell >0 ) {
        spiketimes <- spiketimes + ShiftTime
        lastspiketime <- spiketimes[nspikescell]
        if ( TLast < lastspiketime)
          TLast <- lastspiketime
      }
      
      if (r == 1)
        allspikes[cell] <- list(spiketimes)
      else
        allspikes[cell] <- list(c(allspikes[[cell]],spiketimes))
    }


    ## Now check the end time.
    ##cat(paste("before: EndTime", EndTime, "TLast", TLast, "\n"))
    if (EndTime < TLast)
      EndTime <- TLast

    endtimes[r] <- ( (endclock[r]+1)* 128) +
      ((EndTime %/% mm.WrapTime) * mm.WrapTime)

    if (endtimes[r] < EndTime) {
      endtimes[r] <- endtimes[r] + mm.WrapTime
      cat(paste("wraptime added", mm.WrapTime, "\n"))
    }

    EndTime <- endtimes[r]
    
    ## this is the end of the loop for this record.
    laststop <- seek(fp)                     # used for counting purposes.



  }
  if (spikecount != NSpikes)
    stop(paste ("spikecount differs from expected value",
                spikecount, Nspikes))

  if (eventcount != NEvents)
    stop(paste ("eventcount differs from expected value",
                eventcount, Nevents))

  ## Check the C values
  if (seek(fp) != StatIndex)
    stop(paste ("StatIndex problem", stop, StatIndex))

  C <- readBin(fp, integer(), NCells, intsize, endian="big")
  SpikesInCell <- readBin(fp, integer(), NCells, longsize, endian="big")

  if (( sum(SpikesInCell) != NSpikes))
    stop("Error in the total number of spikes in cell")

  ## Can also check SpikesInCell with the sum of spikes
  count.allspikes <- sapply(allspikes, length)
  if (sum(abs(count.allspikes - SpikesInCell)) > 0)
    stop("Counts of spikes differs...")

  Pe <- readBin(fp, numeric(), NCells, floatsize, endian="big")
  Wi <- readBin(fp, numeric(), NCells, floatsize, endian="big")
  PP <- readBin(fp, numeric(), NCells, floatsize, endian="big")
  WP <- readBin(fp, numeric(), NCells, floatsize, endian="big")
  WW <- readBin(fp, numeric(), NCells, floatsize, endian="big")
  CrossF <- readBin(fp, numeric(), NCells*63, floatsize, endian="big")
  CrossR <- readBin(fp, numeric(), NCells*63, floatsize, endian="big")

  dim(CrossF) <- c(NCells,63)
  dim(CrossR) <- c(NCells,63)
  if ( seek(fp) != filesize)
    stop(paste("difference at end of file", seek(fp), filesize))

  ## End of processing this file.
  close(fp)

  pos <- mm.readpos.compare(NCells, boxes, posfile)
  dists <- make.distances(pos)
  dists.bins <- apply(dists, 2, function(x) { ceiling((x-35)/70) + 1})

  ## Convert spike times into seconds.
  allspikes <- lapply(allspikes, function(x) { x / mm.sample.rate})

  ## check that the spikes are monotonic.
  check.spikes.monotonic(allspikes)
  bursts <- lapply(allspikes, function(x) spikes.to.bursts(x, mm.burst.sep))
  
  res <- list (NFiles=NFiles, NBoxes=NBoxes, NRecords = NRecords,
               NSpikes=NSpikes, NEvents=NEvents,
               startclock=startclock, endclock=endclock,
               nevents=nevents, nspikes=nspikes,
               starttimes=starttimes,
               endtimes=endtimes,
               NCells=NCells, boxes=boxes, C=C,
               spikes=allspikes,
               bursts=bursts,
               CrossF=CrossF, CrossR=CrossR, Pe=Pe,
               file=cellname,
               pos=pos, dists=dists, dists.bins=dists.bins)
  class(res) <- "mm.s"
  res
}

read.ms.mm.data.format1 <- function(cellname, posfile=NULL) {
  ## Read in the multisite data and return a list with all the relevant
  ## data.  This assumes the data is in format 1.

  ## Get the total size of the file so it can be compared with value
  ## of seek() once all the data have been read in.
  filesize <- file.info(cellname)$size

  fp <- file(cellname , 'rb')
  seek(fp,0)
  t <- readBin(fp, integer(), 4, intsize, endian="big")
  
  NFiles <- t[1]; NBoxes <- t[2]; NRecords <- t[3]; NCells <- t[4]

  t <- readBin(fp, integer(), 2, longsize, endian="big")
  NSpikes <- t[1]; NEvents <- t[2];
  
  cat(paste("NFiles", NFiles, "NBoxes", NBoxes, "NRecords", NRecords,
            "NCells", NCells, "\n"))
  cat(paste("NEvents", NEvents, "NSpikes", NSpikes, "\n"))


  ## Read in the file information
  for (r in 1:NFiles) {
    vrn <- readBin(fp, integer(), 1, intsize, endian="big")
    pfilename <- readChar(fp, 64)
    pdirname <- readChar(fp, 64)
    flcrdat  <- readBin(fp, integer(), 1, longsize, endian="big")
    t <- readBin(fp, integer(), 3, intsize, endian="big")
    LowRecord <- t[1]; nrec <- t[2]; LastRecord <- t[3];
    cat(paste("File", r, "name", pfilename, "dir", pdirname,
              "nrec", nrec, "LastRecord", LastRecord, "\n"))
  }

  ## Read in the Box Info ####################
  boxes <- array(0, dim= c(NBoxes, 7))
  for (r in 1:NBoxes) {
    t <- readBin(fp, integer(), 7, intsize, endian="big")
    group <- t[1]; channel <- t[2]; plott <- t[3]
    cat(paste("Box", r, "Group", group, "Chan", channel,
    "Plot", plott, "bounds", t[4], t[5], t[6], t[7], "\n"))
    boxes[r,] <- t
  }

  startclock <- readBin(fp, integer(), NRecords, unsize,signed=F, endian="big")
  endclock   <- readBin(fp, integer(), NRecords, unsize,signed=F, endian="big")

  cat("start times\n")
  print(startclock); print(endclock)
  SpikesInRecords <- readBin(fp, integer(), NRecords, longsize, endian="big")
  SpikesInCell    <- readBin(fp, integer(),   NCells, longsize, endian="big")


  ## This seems to duplicate the information in the boxes.
  C              <- readBin(fp, integer(),   NCells,  intsize, endian="big")
  print(C)
  ##cat(paste("after reading C, file pos is", seek(fp), "\n"))


  ##   Pe <- readBin(fp, numeric(), NCells, doublesize, endian="big")
  ##   print(Pe[1:20])
  ##   Wi <- readBin(fp, numeric(), NCells, doublesize, endian="big")
  ##   PP <- readBin(fp, numeric(), NCells, doublesize, endian="big")
  ##   WP <- readBin(fp, numeric(), NCells, doublesize, endian="big")
  ##   WW <- readBin(fp, numeric(), NCells, doublesize, endian="big")

  ##   CrossF <- readBin(fp, numeric(), NCells*63, doublesize, endian="big")
  ##   CrossR <- readBin(fp, numeric(), NCells*63, doublesize, endian="big")
  ##   dim(CrossF) <- c(NCells,63)
  ##   dim(CrossR) <- c(NCells,63)

  ## double is 10 bytes according to my calculations.
  ## Markus acknowledges that the double format in ThinkC (MAC) is curious
  ## so  for now, I'm just reading in blocks of 10 bytes.  131 is
  ## derived from (63 + 63) + 5 
  tempstuff <- readBin(fp, integer(), NCells*131*10/2, 2, endian="big")
  
  ## Number of spikes from each cell in each record
  N <- readBin(fp, integer(), NCells*NRecords, longsize, endian="big")
  dim(N) <- c(NCells,NRecords)
  ##cat("N\n");print(N)

  ## The N array is useful for knowing which spikes belong to which cell
  ## and which record.
  if (any(apply(N, 2, sum) != SpikesInRecords))
    stop("SpikesInRecords and Col sum of N differ")
  
  ## Time of each spike from each cell in each record
  my.num.spikes <- sum(N)
  cat(paste("my.num.spikes", my.num.spikes, "\n"))
  T <- readBin(fp, integer(), my.num.spikes, longsize, endian="big")

  breaks <- as.vector(N)
  high <- cumsum(breaks)
  low <- c(1, high[1:(length(high)-1)]+1)

  spikes <- apply(rbind(low,high), 2, function (i) {T[i[1]:i[2]]})


  starttimes <- integer(NRecords)       #to be calculated...
  endtimes   <- integer(NRecords)
  nevents    <- integer(NRecords)
  nspikes    <- integer(NRecords)

  ## Make an empty list of size NCells.  Each element will be a list.
  allspikes <- list()
  for (i in 1:NCells) allspikes[i] <- list()
  EndTime <- 0;                        #todo: start of each file?

  for (r in 1:NRecords) {

    ShiftTime <- (EndTime %/% mm.WrapTime) * mm.WrapTime
    while ( ( (startclock[r] * 128) + ShiftTime) < EndTime)
      ShiftTime <- ShiftTime + mm.WrapTime

    starttimes[r] <- (startclock[r] * 128) + ShiftTime
    
    cat(paste(r, "clock", startclock[r], endclock[r], "\n"))

    ## time of each spike from each cell in record r
    TLast <- -1
    for (cell in 1:NCells) {
      ## spikes for this record and cell.
      spiketimes <- spikes[[((r-1)*NCells)+cell]]
      nspikescell <- length(spiketimes)

      if (nspikescell >0 ) {
        spiketimes <- spiketimes + ShiftTime
        lastspiketime <- spiketimes[nspikescell]
        if ( TLast < lastspiketime)
          TLast <- lastspiketime
      }
      if (r == 1)
        allspikes[cell] <- list(spiketimes)
      else
        allspikes[cell] <- list(c(allspikes[[cell]],spiketimes))
    }

    ## Now check the end time.
    ##cat(paste("before: EndTime", EndTime, "TLast", TLast, "\n"))
    if (EndTime < TLast)
      EndTime <- TLast

    endtimes[r] <- ( (endclock[r]+1)* 128) +
      ((EndTime %/% mm.WrapTime) * mm.WrapTime)

    if (endtimes[r] < EndTime) {
      endtimes[r] <- endtimes[r] + mm.WrapTime
      cat(paste("wraptime added", mm.WrapTime, "\n"))
    }

    EndTime <- endtimes[r]
    
  }                                   #next record.
  
  ## Number of events in each record
  EventsInRecord <- readBin(fp, integer(), NRecords, longsize, endian="big")


  my.num.events <- sum(EventsInRecord)
  if (my.num.events > 0) {
    warning(paste("we have some events... oh oh!", my.num.events, "\n"))
    ## todo: need to read in TE, WE, PE if we have any events.
  }
  TE <- readBin(fp, integer(), my.num.events, longsize, endian="big")
  WE <- readBin(fp, integer(), my.num.events, longsize, endian="big")
  PE <- readBin(fp, integer(), my.num.events, longsize, endian="big")

  ## should now be at the end of the file, so can check file length.

  if ( seek(fp) != filesize)
    stop(paste("difference at end of file", seek(fp), filesize))

  ## End of processing this file.
  close(fp)

  pos <- mm.readpos.compare(NCells, boxes, posfile)
  dists <- make.distances(pos)
  dists.bins <- apply(dists, 2, function(x) { ceiling((x-35)/70) + 1})
  hist.r <- table(dists.bins)
  prob.r <- hist.r / sum(hist.r)
  
  allspikes <- lapply(allspikes, function(x) { x / mm.sample.rate})
  
  ## check that the spikes are monotonic.
  check.spikes.monotonic(allspikes)
  bursts <- lapply(allspikes, function(x) spikes.to.bursts(x, mm.burst.sep))
  
  res <- list (NFiles=NFiles, NBoxes=NBoxes, NRecords = NRecords,
               NSpikes=NSpikes, NEvents=NEvents,
               startclock=startclock, endclock=endclock,
               nevents=NEvents, nspikes=NSpikes,
               starttimes=starttimes,
               endtimes=endtimes,
               NCells=NCells, boxes=boxes, C=C,
               spikes=allspikes,
               bursts=bursts,
               T=T,
               ##CrossF=CrossF, CrossR=CrossR, Pe=Pe,
               file=cellname,
               N=N,
               pos=pos, dists=dists, dists.bins=dists.bins,
               prob.r=prob.r
               )
  class(res) <- "mm.s"
  res
}


plot.mm.s <- function(x, whichcells=1:x$NCells, mintime=0,
                     maxtime=max(unlist(x$spikes), na.rm=T),
                      show.bursts=F) {

  ## When evaluating maxtime, some cells may have no spikes; their
  ## lists get converted to NA in unlist() so those NA values need
  ## removing.

  N <- length(whichcells)
  ticpercell <- 1/N; deltay <- ticpercell * 0.9;
  yminadd <- ticpercell

  if (show.bursts)
    spikes <- x$bursts
  else
    spikes <- x$spikes
  
  plot( c(mintime, maxtime), c(0,1), xlab='', type='n', main=x$file)

  ymin <- 0
  for (cell in whichcells) {

    ts <- spikes[[cell]]
    n <- length(ts)
    xs <- double(n*3)
    ymax <- ymin + deltay

    ## NA items allow for breaks in the lines:
    ##
    ## x1   y
    ## x1   y+dy
    ## NA   NA
    ## x2   y
    ## x2   y+dy
    ## NA   NA
    
    xs[ seq(from=1, by=3, length=n) ] <- ts
    xs[ seq(from=2, by=3, length=n) ] <- ts
    xs[ seq(from=3, by=3, length=n) ] <- NA

    ys <- double(n*3)
    ys[ seq(from=1, by=3, length=n) ] <- ymin
    ys[ seq(from=2, by=3, length=n) ] <- ymax
    ys[ seq(from=3, by=3, length=n) ] <- NA
    lines(xs, ys)
    ymin <- ymin + yminadd
  }

  allys <- seq(from=0, by=yminadd, length=N)
  allxs <- 0
  text(allxs, allys, whichcells)
  ##dev.copy2eps(file=paste(cellname, ".ps", sep=''))
}



summary.mms <- function(x) {
  cat(paste("Spike data:", x$file, "\n"))
  cat(paste("NCells", x$NCells, "NRecords", x$NRecords, "\n"))
}

######################################################################
## Jay's functions.
######################################################################

## reload the overlap function if needed.
if (is.loaded(symbol.C("count_overlap"))) dyn.unload("jay/count_overlap.so")
dyn.load("jay/count_overlap.so");

jay.read.spikes <- function(filename, scale=100) {
  ## Read in Jay's data set.
  fp <- file(filename, open="r")
  max.channels <- 64
  channels <- character(max.channels)
  ## first read in the channel names
  num.channels <- 0
  read.channel.names <- 1
  while(read.channel.names) {
    x<-scan(fp, "", n=1, sep='\t', quiet=T)
    ## If first letter of item is not "c" then assume we have now
    ## reached the timestamps.
    if (substr(x,1,1) != 'c') {
      read.channel.names <- 0
      rest <- scan(fp, sep='\t', quiet=T); close(fp)
      ## last element of `rest' is redundant, but we need to keep
      ## x - this is the first element.
      times <- c(as.double(x), rest[1:length(rest)-1])
      ntimes <- length(times)
      dim(times) <- c(num.channels, ntimes/num.channels)
      channels <- channels[1:num.channels] #truncate to right size.
      
    } else {
      ## still reading the channel names.
      num.channels <- num.channels + 1
      if (num.channels > max.channels) {
        stop(paste("num.channels has exceeded max.channels"))
      } else {
        channels[num.channels] <- x
      }
    }
  }

  spikes <- apply(times, 1, filter.for.na)  

  ## Count the number of spikes per channel, and label them.
  nspikes <- sapply(spikes, length)
  names(nspikes) <- channels

  ## Parse the channel names to get the cell positions.
  cols <- as.integer(substring(channels, 4,4)) * scale
  rows <- as.integer(substring(channels, 5,5)) * scale
  pos <- cbind(cols, rows)

  ## check that the spikes are monotonic.
  check.spikes.monotonic(spikes)
  
  res <- list( channels=channels,
              spikes=spikes, nspikes=nspikes, NCells=num.channels,
              file=filename,
              pos=pos)

  res

}

filter.for.na <- function(x) {
  ## Truncate each row of X so that trailing NA entries are removed.
  x.na <- which(is.na(x))
  if (any(x.na))
    x[1:x.na[1]-1]
  else
    x
}


make.distances <- function(posns)
{

  ## POSNS should be a (N,2) array.  Returns a NxN upper triangular
  ## array of the distances between all pairs of cells.
  
  n <- dim(posns)[1]
  dists <- array(0, dim=c(n,n))
  for ( a in 1:n-1)
    for (b in (a+1):n) {
      delta <- posns[a,] - posns[b,]
      dists[a,b] <- sqrt( sum(delta**2))
    }

  dists
}

make.corr.indexes <- function(spikes)
{
  ## Return the correlation index values for each pair of spikes.
  ## The matrix returned is upper triangular.
  ## SPIKES should be a list of length N, N is the number of cells.
  n <- length(spikes)
  dt <- 0.05
  Tmax <- max(unlist(spikes))           #time of last spike.
  corrs <- array(0, dim=c(n,n))
  for ( a in 1:(n-1)) {
    n1 <- length(spikes[[a]])
    for (b in (a+1):n) {
      n2 <- length(spikes[[b]])
      corrs[a,b] <-  (count.nab(spikes[[a]], spikes[[b]]) * Tmax) /
        (n1 * n2 * 0.1)
    }
  }

  corrs
}

make.p.t.cond.r <- function(spikes)
{
  ## Return the correlation index values for each pair of spikes.
  ## The matrix returned is upper triangular.
  ## SPIKES should be a list of length N, N is the number of cells.
  n <- length(spikes)
  n.distances <- max(distance.bins)

  spikepairs <- integer(n.distances)
  nhists <- integer(n.distances)
  tmax <- 2.0;                          #2 second maximum
  n.timebins <- 200

  ## Make a list to store the histograms for each distance bin.
  allhists <- list();
  for (i in 1:n.distances)
    allhists[[i]] <- integer(n.timebins)
                                    
  for (a in 1:(n-1)) {
    n.a <- length(spikes[[a]])
    for (b in (a+1):n) {
      n.b <- length(spikes[[b]])
      bin <- distance.bins[a,b]
      hist <- hist.ab(spikes[[a]], spikes[[b]], tmax, n.timebins)
      allhists[[bin]] <- allhists[[bin]] + hist
      nhists[bin] <- nhists[bin] + 1
      spikepairs[bin] <- spikepairs[bin] + (n.a * n.b)
    }
  }
  list (nhists = nhists,
        allhists = allhists,
        spikepairs = spikepairs)
}



count.nab <- function(ta, tb, tmax=0.05) {
  ## C routine to count the overlap N_ab
  z <- .C("count_overlap",
          as.double(ta),
          as.integer(length(ta)),
          as.double(tb),
          as.integer(length(tb)),
          as.double(tmax),
          res = integer(1))
  z$res
}

hist.ab <- function(ta, tb, tmax, nbins) {
  ## C routine to count the overlap N_ab
  z <- .C("bin_overlap",
          as.double(ta),
          as.integer(length(ta)),
          as.double(tb),
          as.integer(length(tb)),
          as.double(tmax),
          res = integer(nbins),
          as.integer(nbins))
  z$res
}


check.spikes.monotonic <- function(spikes) {
  ## Check to see that all spike times are monotonically increasing.
  ## check.spikes.monotonic( list(c(1,3,5), c(1,5,4)))
  results <- sapply( spikes, function(x) { any(diff(x) <0)})
  if (any(results)) {
    stop(paste("some spike trains are non-monotonic:", which(results),"\n"))
  }
}
  

spikes.to.bursts <- function(spikes, burst.sep=2) {
  ## Convert spikes to bursts.
  ## burst.sep is the threshold time between spikes for finding bursts.
  ## spikes.to.bursts(c(1,2,3, 7,8, 11,12,13,14, 19,20, 23,24))
  f <- which( diff(spikes) > burst.sep) +1
  spikes[c(1,f)]
}


######################################################################
# movie-related functions.




spikes.to.rates <- function(spikes) {
  h <- hist(spikes, breaks=movie.breaks,plot=F)
  h$counts
}



op.picture <- function(pos, rates, iteration) {

  ps.scale <- 0.5 ### 1.0                      #overall scale factor for plot.
  
  ps.min.x <- 40; ps.min.y <- 40
  ps.wid <-  560 * ps.scale; ps.ht <- 560 * ps.scale;
  ps.max.x <- ps.min.x + ps.wid
  ps.max.y <- ps.min.y + ps.ht
  ps.centre.x <- 0.5 * (ps.min.x + ps.max.x)
  ps.centre.y <- 0.5 * (ps.min.y + ps.max.y)
  

  ps.header <- paste("%!PS-Adobe-3.0 EPSF-3.0\n",
                     "%%Title: Main canvas\n",
                     "%%BoundingBox: ", ps.min.x, " ", ps.min.y, " ",
                     ps.max.x, " ", ps.max.y, "\n",
                     "%%CreationDate: ", date(), "\n",
                     "%%EndComments\n\n",
                     ##"%% /d { 3 1 roll   moveto 10.0 div drawbox} def\n\n",
                     "/d { 3 mul 20 min 0 360 arc fill } def\n\n",
                     "%%EndProlog\n%%Page: 1 1\n",
                     ps.centre.x, " ", ps.centre.y, " translate\n", sep='')

  ps.trailer <- "showpage\n%%Trailer"

  this.rates <- rates[iteration,]
  ncells <- length(this.rates)

  fname <- paste("frame", formatC(iteration, width=4,flag="0"), sep='')
  zz <- file(paste(fname,".ps",sep=''), "w")  # open an output file connection
  cat(ps.header, file = zz)
  for (i in 1:ncells) {
   p <- paste(pos[i,1], pos[i,2], this.rates[i], "d\n")
   cat(p, file = zz)
  }

  cat(ps.trailer, file = zz)
  close(zz)

  system(paste("mypstopnm -pbm ", paste(fname,".ps",sep='')))
  system(paste("ppmtogif ", paste(fname,".pbm",sep=''),">",
               paste(fname,".gif",sep='')))
         
  
  fname
}

  
