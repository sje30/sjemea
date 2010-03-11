## Code for reading in data from Alan Litke's array.

litke.spiketimes <- function(file) {
  ## read in the data, and return spike times and array position.
  require(R.matlab)
  z = readMat(file)
  if(length(z)>1)
    stop('too many items in the matlab file')
  
  a = z[[1]]

  names = a[,1]
  posns = a[,2]
  spikes = apply(a, 1, litke.keep.spikes)
  res = list(spikes=spikes, names=names, posns=posns)

}

litke.keep.spikes <- function(row, sample.rate = 20000) {
  ## Return elements 3 to N, where N+1 items are zero.
  ## sample rate of 20Khz is used to convert spikes to time in seconds.
  spikes = row[-(1:2)]
  end = which(spikes == 0)
  spikes = spikes[-end]
  spikes / sample.rate
}

oneoff.make.litke1.layout <- function(file=NULL) {
  ## One-off function needed to create the litke1.txt layout file..
  ## make.litke1.layout('Electrode_Map.mat')
  ## How to hide from external world?  Put dot at start?  or do not export?
  if (is.null(file)) {
    file = system.file("data", "Electrode_Map.mat", package='sjemea')
  }

  require(R.matlab)
  map = readMat(file)

  write.table(map$electrodeMap)
  stopifnot(all.equal(layout[,1], 1:512))
  table = data.frame(electrode=layout[,1], x=layout[,2], y=layout[,3])
  write.table(table, file='/tmp/litke1layout.txt', sep='\t',
              quote=FALSE, row.names=FALSE)
  table
}

make.litke1.layout <- function(positions, names) {
  ## POSITIONS are the names of the electrodes that were recorded.
  ## make the layout for SANGER MEA (cf. make.sanger1.layout)
  ## This is a hexagonal grid.
  xlim <- c(-1000, 1000)
  ylim <- c(-500, 500)
  spacing <- 60
  data(litke1layout)
  
  columns <- match(positions, litke1layout$electrode)
  ## TODO what does match return if no matching electrode found?
  ##
  pos <- cbind(litke1layout$x[columns], litke1layout$y[columns])
  
  rownames(pos) <- names
  
  layout <- list(xlim=xlim, ylim=ylim, spacing=spacing,
                 pos=pos)

  class(layout) <- "mealayout"

  layout

}

show.litke.layout <- function() {
  ## Show the layout of the Litke array.
  data(litke1layout)
  plot(litke1layout$x, litke1layout$y, asp=1,
       type='n', xlab='', ylab='')
  text(litke1layout$x, litke1layout$y, litke1layout$electrode, cex=0.7)
  title('Litke1 layout')
}


litke.read.spikes <- function(filename, ids=NULL,
                              time.interval=1, beg=NULL, end=NULL) {
  ## Read in data from Alan Litke (+ Ben Stafford).
  ## FILENAME: matlab data file
  ## channel per file.
  ## IDS: an optional vector of cell names that should be analysed
  ## -- the other channels are read in but then ignored.


  data <- litke.spiketimes(filename)
  posns <- data$posns
  names <- data$names
  pn <- cbind(names=names, posns=posns)

  spikes <- data$spikes
  names(spikes) <- data$names
  
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
  ## the beg, end range is too narrow).
  spikes <- remove.empty.channels(spikes)

  if (!is.null(ids)) {
    ## Filter out some channel names, either inclusively or exclusively.
    spikes <- filter.channel.names(spikes, ids)
  }
  
  rec.time <- c(beg, end)

  channels <- names(spikes)
  ## now go back and find the positions corresponding to these names
  ## that survived.

  columns <- match(channels, pn[,"names"])
  positions <- pn[columns,"posns"]
  
  ## Count the number of spikes per channel, and label them.
  nspikes <- sapply(spikes, length)
  ## meanfiring rate is the number of spikes divided by the (time of
  ## last spike - time of first spike).  
  meanfiringrate <- nspikes/ ( end - beg)

  rates <- make.spikes.to.frate(spikes, time.interval=time.interval,
                                beg=beg, end=end)

  layout <- make.litke1.layout(positions, names(spikes))

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

  litke.breaks = seq(from=0, to=2050, by=50)
  res$corr = corr.index(res, litke.breaks)

  res

}

