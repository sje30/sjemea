## Read the multisite data.
## Sat 25 Aug 2001

## showConnections(all=T)

## If 2nd argument to seek is omitted, it just returns the current position.

longsize <- 4; floatsize <- 4; intsize  <- 2; unsize   <- 2

read.ms.data <- function(cellname) {
  ## Read in the multisite data and return a list with all the relevant
  ## data.
  
  filesize <- file.info(cellname)$size
  fp <- file(cellname , 'rb')
  seek(fp,0)
  Format <- readBin(fp, integer(), 1, longsize, endian="big")
  t <- readBin(fp, integer(), 4, longsize, endian="big")
  FileIndex <- t[1]; BoxIndex <- t[2]; RecIndex <- t[3]; StatIndex <- t[4]

  if (Format != 2) {
    warning(paste("Format not equal to 2",Format))
  }

  ## Now read the NFiles...
  seek(fp, 64)
  t <- readBin(fp, integer(), 4, intsize, endian="big")
  NFiles <- t[1]; NBoxes <- t[2]; NRecords <- t[3]; NCells <- t[4]

  t <- readBin(fp, integer(), 2, longsize, endian="big")
  NEvents <- t[1]; NSpikes <- t[2]

  ## Read in the fileinfo.
  if (seek(fp) != FileIndex)
    warning("error - current file position different from expected FileIndex")

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
    warning("error - current file position different from BoxIndex")

  ## Make an array to store the boxes.
  boxes <- array(0, dim= c(NBoxes, 4))
  ## todo -- determine how these boxes relate to position of neurons.
  seek(fp, BoxIndex)
  for (r in 1:NBoxes) {
    t <- readBin(fp, integer(), 7, intsize, endian="big")
    group <- t[1]; channel <- t[2]; plott <- t[3]
    cat(paste("Box", r, "Group", group, "Chan", channel,
              "Plot", plott, "bounds", t[4], t[5], t[6], t[7], "\n"))
    boxes[r,] <- t[4:7]
  }


  ## now parse the RecIndex... ####################

  if (seek(fp) != RecIndex)
    warning(paste("seek position different from expected RecIndex",
                  seek(fp), RecIndex))

  ## RecIndex points to an array of length NRecords, each which points
  ## to the start of the rth record.

  seek(fp,RecIndex)
  RecordIndexes <- readBin(fp, integer(), NRecords, longsize, endian="big")

  ## Parse each record...

  spikecount <- 0
  eventcount <- 0
  laststop   <- 0

  ## Make an empty list of size NCells.  Each element will be a list.
  allspikes <- list()
  for (i in 1:NCells) {
    allspikes[i] <- list()
  }

  for (r in 1:NRecords) {
    if ((laststop >0) && (seek(fp) != laststop)) {
      warning("Error: RecordIndex and position of last byte read do not match")
      warning(paste("Record", r, "start", start, "laststop", laststop))
    }
    
    seek(fp, RecordIndexes[r])
    startclock <- readBin(fp, integer(), 1, unsize, endian="big")
    endclock   <- readBin(fp, integer(), 1, unsize, endian="big")
    nevents    <- readBin(fp, integer(), 1, longsize, endian="big")
    nspikes    <- readBin(fp, integer(), 1, longsize, endian="big")


    cat(paste(r, "clock", startclock, endclock, "#events", nevents,
              "#spikes", nspikes, "\n"))
    spikecount <- spikecount + nspikes
    eventcount <- eventcount + nevents

    ## Read in the number of spikes from each cell in record r.
    spikespercell <- readBin(fp, integer(), NCells, longsize, endian="big")

    ## Read in the time of event.
    eventsinrecord <- readBin(fp, integer(), nevents, longsize, endian="big")

    ## width of events
    we <- readBin(fp, integer(), nevents, intsize, endian="big")

    ## peak of events
    pe <- readBin(fp, integer(), nevents, intsize, endian="big")

    ## time of each spike from each cell in record r
    for (cell in 1:NCells) {
      nspikescell <- spikespercell[cell]
      spiketimes <- readBin(fp, integer(), nspikescell, longsize, endian="big")
      if (r == 1)
        allspikes[cell] <- list(spiketimes)
      else
        allspikes[cell] <- list(c(allspikes[[cell]],spiketimes))
    }
    ## this is the end of the loop for this record.
    laststop <- seek(fp)                     # used for counting purposes.
  }
  if (spikecount != NSpikes)
    warning(paste ("spikecount differs from expected value",
                   spikecount, Nspikes))

  if (eventcount != NEvents)
    warning(paste ("eventcount differs from expected value",
                   eventcount, Nevents))

  ## Chceck the C values
  if (seek(fp) != StatIndex)
    warning(paste ("StatIndex problem", stop, StatIndex))

  C <- readBin(fp, integer(), NCells, intsize, endian="big")
  SpikesInCell <- readBin(fp, integer(), NCells, longsize, endian="big")

  if (( sum(SpikesInCell) != NSpikes))
    warning("Error in the total number of spikes in cell")

  ## Can also check SpikesInCell with the sum of spikes
  count.allspikes <- sapply(allspikes, length)
  if (sum(abs(count.allspikes - SpikesInCell)) > 0)
    warning("Counts of spikes differs...")

  Pe <- readBin(fp, numeric(), NCells, floatsize, endian="big")
  Wi <- readBin(fp, numeric(), NCells, floatsize, endian="big")
  PP <- readBin(fp, numeric(), NCells, floatsize, endian="big")
  WP <- readBin(fp, numeric(), NCells, floatsize, endian="big")
  WW <- readBin(fp, numeric(), NCells, floatsize, endian="big")
  CrossF <- readBin(fp, numeric(), NCells*63, floatsize, endian="big")
  CrossR <- readBin(fp, numeric(), NCells*63, floatsize, endian="big")

  if ( seek(fp) != filesize)
    warning(paste("difference at end of file", seek(fp), filesize))

  ## End of processing this file.
  close(fp)

  res <- list (NFiles=NFiles, NBoxes=NBoxes, NRecords = NRecords,
               NCells=NCells, boxes=boxes, C=C,
               allspikes=allspikes,
               file=cellname)

  class(res) <- "mms"
  res
}



cellname <- 'c08(b08,9)'
cellname <- 'P0,8.19.9l/c07,17'
cellname <- 'P1,_8.20.91/c06'

s <- read.ms.data(cellname)

plot.mms(s)


plot.mms <- function(x, whichcells=1:x$NCells) {
  maxtime <- 8463658                    #todo -- what is this number really?
  maxtime <-13000000
  mintime <- 0

  deltay  <- 0.05
  yminadd <- deltay+0.02
  N <- length(whichcells)
  ticpercell <- 1/N; deltay <- ticpercell * 0.9;
  yminadd <- ticpercell
  ##deltay  <- 0.05; yminadd <- deltay+0.02

  plot( c(mintime, maxtime), c(0,1), xlab='', type='n', main=x$file)

  ymin <- 0
  for (cell in whichcells) {
    ts <- x$allspikes[[cell]]
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






