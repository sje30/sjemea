sun.read.spikes <- function(filename, ids=NULL,
                            time.interval=1,
                            beg=NULL, end=NULL,
                            min.rate=0) {
  ## Read in Sun's data sets.  (e.g. Sun et al, PNAS 2008: beta 2).

  ## This is adaptred from jay.read.spikes, as the file formats are
  ## very similar; the main difference being the way channels are
  ## named, and that electrodes are spaced 200um apart.
  
  ## IDS is an optional vector of cell numbers that should be analysed
  ## -- the other channels are read in but then ignored.

  ## Read in all the data at once, and then separate into channel
  ## names and spike times.
  

  data <- scan(filename, what=character(0), sep='\t')
  
  ## find out which data is a channel name, and which is 
  headers <- grep('^[a-zA-Z]', data)
  nchannels <- length(headers)
  ## safety check: headers must increment by 1
  stopifnot( nchannels==1 ||
            isTRUE(all.equal(diff(headers), rep(1, nchannels-1))))
  
  ## can have multiple units on one location, given by 'a', 'b' etc.
  
  channels <- data[headers]

  ## Now get the spike times.
  data <- data[-headers]                  #take off the channel names
  
  npts <- length(data)
  data <- data[-npts]             #get rid of trailing CR at end of file.
  npts <- npts -1
  
  data <- as.double(data)
  ## reshape data into matrix - one col is one electrode
  
  dim(data) <- c(nchannels, npts/nchannels)
  data <- t(data)
  colnames(data) <- channels
  

  spikes <- apply(data, 2, jay.filter.for.na)

  if (!is.null(end))
    spikes <- lapply(spikes, jay.filter.for.max, max=end)

  if (!is.null(beg))
    spikes <- lapply(spikes, jay.filter.for.min, min=beg)


  
  if (!is.null(ids) ) {
    if (any(ids>length(spikes)))
      stop(paste("some ids not in this data set:",
                 paste(ids[ids>length(spikes)],collapse=" ")))
    
    spikes <- spikes[ids];
    channels <- channels[ids];
  }

  spikes.range <- range(unlist(spikes))
  if (is.null(beg))  beg <-  spikes.range[1]
  if (is.null(end))  end <-  spikes.range[2]
  rec.time <- c(beg, end)
  if (min.rate > 0 ) {
    
    ## Check for inactive channels.
    nspikes <- sapply(spikes,length)
    durn <- diff(rec.time)
    rates <- nspikes/durn
    inactive <- which(rates < min.rate)
    if (any(inactive)) {
      paste("Removing spikes with low firing rates: ",
            paste(inactive, collapse=' '))
      spikes   =   spikes[-inactive]
      channels = channels[-inactive]
    }
    
    
  }


  
  ## Count the number of spikes per channel, and label them.
  nspikes <- sapply(spikes, length)
  names(nspikes) <- channels

  ## meanfiring rate is the number of spikes divided by the (time of
  ## last spike - time of first spike).  
  meanfiringrate <- nspikes/ ( end - beg)

  ## Parse the channel names to get the cell positions.

  layout <- make.sun.layout( channels)
  
  ## temporary test: shuffle electrode positions.
  ## pos <- pos[sample(1:num.channels),]

  unit.offsets <- NULL                  #default value.
  ## check that the spikes are monotonic.
  check.spikes.monotonic(spikes)

  rates <- make.spikes.to.frate(spikes, time.interval=time.interval,
                                beg=beg, end=end)
  
  ## Ignore any shifting of cells that were assigned to the same
  ## electrode.
  
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

  distance.breaks = c(0, 150, 250, 350, 450, 550, 650, 1000, 
    2000)  
  res$corr = corr.index(res, distance.breaks)

  res


}
sun2.read.spikes <- function(filename, ids=NULL,
                            time.interval=1,
                            beg=NULL, end=NULL,
                            min.rate=0) {
  ## Read in Sun's data sets.  (e.g. Sun et al, PNAS 2008: beta 2).

  ## This is adaptred from jay.read.spikes, as the file formats are
  ## very similar; the main difference being the way channels are
  ## named, and that electrodes are spaced 200um apart.
  
  ## IDS is an optional vector of cell numbers that should be analysed
  ## -- the other channels are read in but then ignored.

  ## Read in all the data at once, and then separate into channel
  ## names and spike times.
  

  data <- scan(filename, what=character(0), sep='\t')
  
  ## find out which data is a channel name, and which is 
  headers <- grep('^[a-zA-Z]', data)
  nchannels <- length(headers)
  ## safety check: headers must increment by 1
  stopifnot( nchannels==1 ||
            isTRUE(all.equal(diff(headers), rep(1, nchannels-1))))
  
  ## can have multiple units on one location, given by 'a', 'b' etc.
  
  channels <- data[headers]

  ## Now get the spike times.
  data <- data[-headers]                  #take off the channel names
  
  npts <- length(data)
  data <- data[-npts]             #get rid of trailing CR at end of file.
  npts <- npts -1
  
  data <- as.double(data)
  ## reshape data into matrix - one col is one electrode
  
  dim(data) <- c(nchannels, npts/nchannels)
  data <- t(data)
  colnames(data) <- channels
  

  spikes <- apply(data, 2, jay.filter.for.na)

  if (!is.null(end))
    spikes <- lapply(spikes, jay.filter.for.max, max=end)

  if (!is.null(beg))
    spikes <- lapply(spikes, jay.filter.for.min, min=beg)


  
  if (!is.null(ids) ) {
    if (any(ids>length(spikes)))
      stop(paste("some ids not in this data set:",
                 paste(ids[ids>length(spikes)],collapse=" ")))
    
    spikes <- spikes[ids];
    channels <- channels[ids];
  }

  spikes.range <- range(unlist(spikes))
  if (is.null(beg))  beg <-  spikes.range[1]
  if (is.null(end))  end <-  spikes.range[2]
  rec.time <- c(beg, end)
  if (min.rate > 0 ) {
    
    ## Check for inactive channels.
    nspikes <- sapply(spikes,length)
    durn <- diff(rec.time)
    rates <- nspikes/durn
    inactive <- which(rates < min.rate)
    if (any(inactive)) {
      paste("Removing spikes with low firing rates: ",
            paste(inactive, collapse=' '))
      spikes   =   spikes[-inactive]
      channels = channels[-inactive]
    }
    
    
  }


  
  ## Count the number of spikes per channel, and label them.
  nspikes <- sapply(spikes, length)
  names(nspikes) <- channels

  ## meanfiring rate is the number of spikes divided by the (time of
  ## last spike - time of first spike).  
  meanfiringrate <- nspikes/ ( sapply(spikes, max) - sapply(spikes, min))

  ## Parse the channel names to get the cell positions.

  layout <- make.sun.layout( channels)
  
  ## temporary test: shuffle electrode positions.
  ## pos <- pos[sample(1:num.channels),]
  
  ## check that the spikes are monotonic.
  check.spikes.monotonic(spikes)

  rates <- make.spikes.to.frate(spikes, time.interval=time.interval,
                                beg=beg, end=end)
  
  ## Ignore any shifting of cells that were assigned to the same
  ## electrode.
  
  res <- list( channels=channels,
              spikes=spikes, nspikes=nspikes, NCells=length(spikes),
              meanfiringrate=meanfiringrate,
              file=filename,
              layout=layout,
              rates=rates,
              rec.time=rec.time
              )
  class(res) <- "mm.s"

  distance.breaks = c(0, 150, 250, 350, 450, 550, 650, 1000, 
    2000)  
  res$corr = corr.index(res, distance.breaks)

  res


}


make.sun.layout <- function(positions) {
  ## make the layout for Sun's MEA
  ## This is like Jay's, but with 200 um spacing. 

  xlim <- ylim <- c(50, 1700)
  spacing <- 200

  electrode <- as.integer(substring(positions, 6,7))
  cols <- as.integer(substring(positions, 6,6)) * spacing
  rows <- as.integer(substring(positions, 7,7)) * spacing

  pos <- cbind(x=cols, y=rows, electrode.num=electrode)
  
  rownames(pos) <- positions
  layout <- list(xlim=xlim, ylim=ylim, spacing=spacing,
                 pos=pos)

  class(layout) <- "mealayout"

  layout

}

make.sun.layout.old <- function(positions) {
  ## make the layout for Sun's MEA
  ## This is like Jay's, but with 200 um spacing. 

  xlim <- ylim <- c(50, 1700)
  spacing <- 200

  cols <- as.integer(substring(positions, 6,6)) * spacing
  rows <- as.integer(substring(positions, 7,7)) * spacing
  pos <- cbind(cols, rows)
  
  rownames(pos) <- positions
  layout <- list(xlim=xlim, ylim=ylim, spacing=spacing,
                 pos=pos)

  class(layout) <- "mealayout"

  layout

}

