## IIT.R code

iit.read.spikes <- function(filename, ids=NULL,
                              time.interval=1, beg=NULL, end=NULL) {
  ## Read in data from IIT.
  ## FILENAME: spike times, stored in matlab sparse arrays.
  ## IDS: an optional vector of cell numbers that should be analysed
  ## -- the other channels are read in but then ignored.

  if (any(grep('mat$', filename))) {
    if (!require(R.matlab))
      stop('The R.matlab package is needed for this routine.  Please install.')
    
    z <- readMat(filename)

    ## 2010-05-14: some channels are included, but empty, so let's remove
    ## them.
    empty.channels <- which(sapply(z, length)==0)
    if (any(empty.channels))
      z <- z[-empty.channels]
    
    frame.rates = 7702  ## was 7800.
    ## each element - E- is a sparse matrix with one column, so find out
    ## where the non-zero elements are.
    ##spikes <- lapply(z, function(e) { which(e[,1]>0)/frame.rates})
    ## this is a hack right now, not sure why we need to add 1 to the sparse
    ## index values, but it works!
    spikes <- lapply(z, function(e) { (e@i+1)/frame.rates})
  } else {
    ## Assume the file is a regular csv and just read it in.
    spikes <- read.csv(filename)
    spikes <- lapply(spikes, jay.filter.for.na)
  }

  ## Now remove spikes outside of range required.
  spikes.range <- range(unlist(spikes))
  if (!is.null(end)) {
    spikes <- lapply(spikes, jay.filter.for.max, max=end)
  } else {
    end <- spikes.range[2]
  }

  if (!is.null(beg)) {
    spikes <- lapply(spikes, jay.filter.for.min, min=beg)
  } else {
    beg <- spikes.range[1]
  }

  ## Remove any channels that have zero spikes (which can happen if
  ## the beg, end range is too narrow, or if a datafile is empty,
  ## which happens sometime for the Feller data.
  spikes <- remove.empty.channels(spikes)


  if (!is.null(ids)) {
    ## Filter out some channel names, either inclusively or exclusively.
    spikes <- filter.channel.names(spikes, ids)
  }
  
  rec.time <- c(beg, end)

  channels <- names(spikes)

  
  ## Count the number of spikes per channel, and label them.
  nspikes <- sapply(spikes, length)

  ## meanfiring rate is the number of spikes divided by the (time of
  ## last spike - time of first spike).  
  meanfiringrate <- nspikes/ ( end - beg)

  rates <- make.spikes.to.frate(spikes, time.interval=time.interval,
                                beg=beg, end=end)

  ## Parse the channel names to get the cell positions.
  layout <- make.iit.layout(channels)

  ## TODO; worry about multiple units recorded at the same location?
  unit.offsets <- NULL                  #default value.
  
  ## check that the spikes are monotonic.
  check.spikes.monotonic(spikes)


  rates <- make.spikes.to.frate(spikes, time.interval=time.interval,
                                beg=beg, end=end)

  
  res <- list( channels=channels,
              spikes=spikes, nspikes=nspikes, NCells=length(spikes),
              meanfiringrate=meanfiringrate,
              file=filename,
              layout=layout,
              rates=rates,
              unit.offsets=unit.offsets,
              rec.time=rec.time
              )
  class(res) <- "mm.s"

  max.dist <- sqrt(2*(64*42)^2)
  breaks = seq(from=0, to=max.dist, by=50)
  res$corr = corr.index(res, breaks)

  res

}


make.iit.layout <- function(positions) {
  ## make the layout for SANGER MEA (cf. make.sanger1.layout)
  ## This is a hexagonal grid.
  spacing <- 42;  #20um diam electrodes
  xlim <- ylim <- c(0, 64*spacing)



  ## parse the channel names into rows and columns.
  regexp <- "Ch([0-9]+)\\.([0-9]+)"
  r <- as.integer(gsub(regexp, "\\1", positions))
  c <- as.integer(gsub(regexp, "\\2", positions))

  rows = (r-1)*spacing
  cols = (c-1)*spacing
  electrode.num <- ((r-1)*64) + c ## start electrodes from number 1
  
  pos <- cbind(x=rows, y=cols, electrode.num=electrode.num)
  
  rownames(pos) <- positions
  
  layout <- list(xlim=xlim, ylim=ylim, spacing=spacing,
                 pos=pos)

  class(layout) <- "mealayout"

  layout

}


