## Code for reading in Data from Marla Feller's group.
## 2009-07-21

## The layout for EJ's array can be found in  ejc_ElectrodeNumbering.jpg
## This file is installed as part of the R package.
## system.file("examples/ejc_ElectrodeNumbering.jpg",package="sjemea")

make.ejc.layout <- function(positions) {
  ## make the layout for SANGER MEA (cf. make.sanger1.layout)
  ## This is a hexagonal grid.
  xlim <- ylim <- c(-300, 300)
  spacing <- 60

  cols.an <- toupper(substring(positions,1,2))
  columns <- match(cols.an, ejcmealayout$name)
  ## round columns to nearest integer - i.e. um.
  pos <- cbind(x=as.integer(ejcmealayout$x[columns]),
               y=as.integer(ejcmealayout$y[columns]),
               electrode.num=ejcmealayout$number[columns])
  
  rownames(pos) <- positions
  
  layout <- list(xlim=xlim, ylim=ylim, spacing=spacing,
                 pos=pos)

  class(layout) <- "mealayout"

  layout

}


make.ejcmealayout <- function(file=NULL,sep=60,
                            write.table=FALSE) {
  ## Create the data file that has the array layout for all electrodes.
  ## This will not normally need to be called, as the data file is stored
  ## in the package, and can be accessed using:
  ## data(ejcmealayout)
  ##
  ## FILE: location of the array layout file.
  ## SEP: electrode separation.
  ## write.table: set to TRUE if you want the file written to the /tmp
  ## directory.

  ## if file is not specified, get the file from the data directory.
  if (is.null(file)) {
    file <- system.file("data", "ejc_layout.txt", package = "sjemea")
  }
    
  table <- read.csv(file, sep='\t')

  ## convert the row, col, fields into x,y position, assuming a
  ## hexagonal grid layout; odd-numbered rows need shifting left.

  dx <-  2 * cos(pi/6) * sep
  x <- table$c * dx
  x <- x - ifelse( (table$r %%2) == 1, dx/2, 0)

  y <- table$r * sep * cos(pi/3)
  table <- cbind(table, x=round(x,2), y=round(y,2))
  table

  if (write.table) {
    write.table(table, file='/tmp/ejcmealayout.txt', sep='\t',
                quote=FALSE, row.names=FALSE)
  }

  table
}


ejcmealayout.check <- function() {
  ## Check the layout looks okay.
  data(ejcmealayout)
  with(ejcmealayout, {
       plot(x, y, pch=19, cex=.2, asp=1, main='ejcmealayout',
            xlab='', ylab='', bty='n')
       text(x, y, paste(number, name))
     })
  
}


feller.spiketimes <- function(dir) {
  ## Read all the spike data contained in DIR, and return as a list.

  filenames <- dir(path=dir)

  ## TODO:
  ## Assume last three characters contain the electrode position.  How
  ## valid is this?  This will be something like "c3a" or "c3b" with
  ## the first two chars representing the location.
  electrodes <- substring(filenames, first=nchar(filenames)-2)

  ## Check if there any duplicate electrodes; this would indicate that
  ## there maybe several experiments patched together into same file.
  if (any(duplicated(electrodes)))
    stop("You have duplicate electrodes in ", dir)
  

  ## perhaps best to order electrodes by their number, rather than 2
  ## char id, as the number varies fairly smoothly with position on the
  ## array.

  ## Convert each electrode name in the filename into upper case
  ## (dropping any trailing a/b/c for multiple cells on a channel) and
  ## find equivalent channel number.

  data(ejcmealayout)

  ## take first two letters of electrode name to get the AN value - A
  ## is an alphabetic char [A-G] and N is a number [1-8].
  
  electrode.an <- toupper(substring(electrodes, 1,2))
  channel.num <- match( electrode.an, ejcmealayout$name)
  channel.num <- ejcmealayout$number[channel.num]
  ## If we load in the data, sorted by channel number, nearby neurons
  ## should look most similar.

  ## When ranking by channel.nu, there will be ties!  Need to check that these
  ## are resolved okay. TODO
  rank <-  rank(channel.num)
  channels <- data.frame(file=filenames, electrode=electrodes,
                         num=channel.num, rank=rank)

  channels.ordered <- channels[order(rank), ]

  ## channels.ordered is a nice dataframe that shows you what is about
  ## to be read in.
  print(channels.ordered)
  
  ## set to 1 if you want to check that the data is read in okay. Else
  ## divide by 20,000 to get time in seconds.
  
  to.secs <-  1/20000
  
  spikes <- lapply(as.character(channels.ordered$file), function(f) {
    file <- paste(dir, f, sep='/')
    vals <- scan(file, quiet=TRUE) * to.secs
    vals
  })

  names(spikes) <- channels.ordered$electrode

  ## todo: could check that last 3 chars of filename match names(spikes)
  spikes
}

remove.empty.channels <- function(spikes) {
  ## Remove any spike trains that are empty, i.e. zero spikes in them.
  ## This can happen if the beg, end range is too narrow, or if a
  ## datafile is empty, which happens sometime for the Feller data.
  ## TODO: currentlly only used by the Feller reader, perhaps it could
  ## also be used for other routines too?
  
  nspikes <- sapply(spikes, length)
  empty <- which(nspikes==0)
  if ( any(empty) ) {
    spikes <- spikes[-empty]
  }
  spikes 
}


names.to.indexes <- function(names, elems, allow.na=FALSE) {
  ## Return the indexes of where each element of ELEMS is within NAMES.
  ## If the first element of ELEMS is '-', then return all indexes except
  ## those matching ELEMS.
  ## Example:
  ## names = c('a', 'b', 'c', 'd', 'e')
  ## names.to.indexes(names, c('d', 'b', 'a'))  ## 4 2 1
  ## names.to.indexes(names, c( '-', 'c', 'a')) ## 2 4 5

  ## to check if first element is "-", we have to use this more
  ## complex expression, as elems[1] == "-" is an error if the first element
  ## by chance is NA.
  if ( isTRUE(all.equal("-", elems[1])) ) {
    invert = TRUE
    elems = elems[-1]
  } else {
    invert = FALSE

  }

  indexes = match(elems, names)
  if (!allow.na) {
    if (any(is.na(indexes)))
      stop('some indexes not found.')
  }
  
  if (invert)
    indexes = setdiff(1:(length(names)), indexes)

  indexes
  
}

filter.channel.names <- function(spikes, ids) {
  ## Filter out some channel names.
  ## Keep only the channels mentioned in IDS.
  ## If the elements of IDS are numeric, they are assumed to be the
  ## indexes of the spike trains; otherwise, they are assumed to be the 
  ## names of cells.
  ## e.g.
  ## spikes2 <- filter.channel.names(spikes, c('-', 'g4a', 'a6a'))
  ## spikes2 <- filter.channel.names(spikes, c('g4a', 'a6a'))
  ## spikes2 <- filter.channel.names(spikes, c(5, 3, 1))
  ## first call throws away two channels; second call keeps just two channels.
  ## third just keeps the three trains mentioned.

  if (any(is.character(ids)))
    ids = names.to.indexes(names(spikes), ids)
  
  spikes[ids]
}
             
  


feller.read.spikes <- function(filename, ids=NULL,
                               time.interval=1, beg=NULL, end=NULL) {
  ## Read in data from Marla Feller.
  ## FILENAME: directory that contains the spike times, one
  ## channel per file.
  ## IDS: an optional vector of cell numbers that should be analysed
  ## -- the other channels are read in but then ignored.


  spikes <- feller.spiketimes(filename)

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
  layout <- make.ejc.layout(channels)

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

  feller.breaks = seq(from=0, to=500, by=50)
  res$corr = corr.index(res, feller.breaks)

  res

}


