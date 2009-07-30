## Code for reading in Data from Marla Feller's group.
## 2009-07-21

make.ejc.layout <- function(positions) {
  ## make the layout for SANGER MEA (cf. make.sanger1.layout)
  ## This is a hexagonal grid.
  xlim <- ylim <- c(-300, 300)
  spacing <- 60

  cols.an <- toupper(substring(positions,1,2))
  columns <- match(cols.an, ejcmealayout$name)
  pos <- cbind(ejcmealayout$x[columns], ejcmealayout$y[columns])
  
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

  browser()
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


