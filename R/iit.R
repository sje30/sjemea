## IIT.R code


######################################################################
## APS version of funcitons, using hdf5 or .Rdata format
######################################################################
aps.read.spikes <- function(filename, ids=NULL,
                             min.rate=-Inf, max.rate=Inf,
                             time.interval=1, beg=NULL, end=NULL) {
  ## Read in data from APS array
  ## FILENAME: spike times, stored in either HDF5 or RData.
  ## IDS: an optional vector of cell numbers that should be analysed
  ## -- the other channels are read in but then ignored.

  h5file <- any(grep('\\.h5$', filename))
  
  if (h5file) {
    require(rhdf5)
    ## read in the hdf5 file.
    obj <- h5read(path.expand(filename), '/')
    obj$spiketimes <- chop(obj$spikes, obj$sCount)
  } else {
    ## Assume it is an Rdata file.
    load(filename)
    obj <- x$x[,,1]
  }
  
  ## Find out which channels are valid; keep only these ones.
  validChannels <- as.vector(obj$validChannels)

  
  spikes2 <- obj$spiketimes[validChannels]
  spikes <- sapply(spikes2, as.vector)

  epos <- obj$epos
  epos <- epos[validChannels,]
  epos2 <- paste("Ch", epos[,1], ".", epos[,2], sep='')
  names(spikes) <- epos2

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


  
  ## Count the number of spikes per channel, and label them.
  nspikes <- sapply(spikes, length)

  ## meanfiring rate is the number of spikes divided by the (time of
  ## last spike - time of first spike).  
  meanfiringrate <- nspikes/ ( end - beg)

  ok.trains <- (meanfiringrate >= min.rate) & (meanfiringrate <= max.rate)
  poor.firing <- which(!ok.trains)
  if (any(poor.firing)) {
    ## Remove the spike trains that are firing too low or too high.
    printf("Removing %d spike trains with low/high firing rates\n",
           length(poor.firing))
    spikes <- spikes[ok.trains]
    nspikes <- nspikes[-poor.firing]
    meanfiringrate <- meanfiringrate[-poor.firing]
  }

  

  rates <- make.spikes.to.frate(spikes, time.interval=time.interval,
                                beg=beg, end=end)


  ## TODO; worry about multiple units recorded at the same location?


  ## Now compute the layout, now that we now which cells are left.
  channels <- names(spikes)
  layout <- make.iit.layout(channels)

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


make.aps.layout <- function(positions) {
  ## make the layout for APS MEA
  ## This is a regular square grid.
  spacing <- 42;  #20um diam electrodes
  xlim <- ylim <- c(0, 64*spacing)


  r <- positions[,1]
  c <- positions[,2]

  rows = (r-1)*spacing
  cols = (c-1)*spacing
  electrode.num <- ((r-1)*64) + c ## start electrodes from number 1
  
  pos <- cbind(x=rows, y=cols, electrode.num=electrode.num)

  rownames(pos) <- paste(positions[,1], positions[,2], sep='_')
  
  layout <- list(xlim=xlim, ylim=ylim, spacing=spacing,
                 pos=pos)
  
  class(layout) <- "mealayout"

  layout

}



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


