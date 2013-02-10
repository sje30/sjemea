## HDFreader functions, and helpers.
## 2013-01-04

h5.read.spikes <- function(h5file, ids=NULL,
                               time.interval=1, beg=NULL, end=NULL, corr.breaks) {
  # Read in a HDF5 file.

  chop <- function(v, counts) {
    ## chop(9:1, c(3,2,4))
    stopifnot(sum(counts)==length(v))
    end <- cumsum(counts)
    beg <- c(1, 1+end[-length(end)])
    begend <- cbind(beg, end)
    apply(begend, 1, function(x) v[x[1]:x[2]])
  }


  ## data$spikes is a 1d array, rather than a vector, so convert it to 1d vector below.
  data <- h5read(path.expand(h5file),
                 name='/')       #read in all of data at once.


  ## unroll the spikes data structure.
  spikes <- chop(as.vector(data$spikes), data$sCount)  
  names(spikes) <- data$names
  

  arrayinfo <- get.array.info(data)
  layout <- arrayinfo$layout
  if(missing(corr.breaks)) {
    corr.breaks <-   arrayinfo$corr.breaks
  }
  
  s <- construct.s(spikes, ids, time.interval, beg, end, corr.breaks, layout, filename=h5file)
}



construct.s <- function(spikes, ids, time.interval, beg, end,
                        corr.breaks, layout, filename) {

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

  ## No more spike trains should be removed, so we can refilter the layout to ensure
  ## that only spikes that are left in the array are kept in the layout.
  channels <- names(spikes)
  keep <- match(channels, rownames(layout$pos))
  layout$pos <- layout$pos[keep,]
  
  rec.time <- c(beg, end)

  ##channels <- names(spikes)

  
  ## Count the number of spikes per channel, and label them.
  nspikes <- sapply(spikes, length)

  ## meanfiring rate is the number of spikes divided by the (time of
  ## last spike - time of first spike).  
  meanfiringrate <- nspikes/ ( end - beg)

  rates <- make.spikes.to.frate(spikes, time.interval=time.interval,
                                beg=beg, end=end)

  ## TODO; worry about multiple units recorded at the same location?
  unit.offsets <- NULL                  #default value.
  
  ## check that the spikes are monotonic.
  check.spikes.monotonic(spikes)


  
  res <- list( channels=names(spikes),
              spikes=spikes, nspikes=nspikes, NCells=length(spikes),
              meanfiringrate=meanfiringrate,
              file=filename,
              layout=layout,
              rates=rates,
              unit.offsets=unit.offsets,
              rec.time=rec.time
              )
  class(res) <- "mm.s"

  if(length(corr.breaks) == 1) {
    res$corr = NULL
  }  else {
    res$corr = corr.index(res, corr.breaks)
  }

  res

}
