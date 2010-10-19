## ncl_mea.R --- Specifics of analysising MEA data from Newcastle.
## Author: Stephen J Eglen
## Copyright: GPL

mcd.data.to.array <- function(file, beg=NULL, end=NULL) {

  ## Read in the MCD data file.  Return the spike trains and the channel
  ## names that are used.
  ## Spikes outside of the time window [BEG, END] (if non-null) are ignored.
  ## BEG and END are given in seconds.
  ## Recall that blank lines will not be read into the data matrix.
  
  data <- read.table(file, as.is=TRUE, skip=2, fill=TRUE)
  

  ## Row number of the start of each spike train.
  channel.start <- which(data[,2] =='Spikes')
  n.channels <- length(channel.start)
  
  ## For the last channel, we add where the end of the file is.
  ## This is so that our check for no spikes will also work for the last
  ## channel in the data set.
  channel.start <- c(channel.start, nrow(data))

  ## Name of each channel.
  channels <- data[channel.start,4]
  channel.ok <- rep(0, n.channels)
  spikes <- list(); n <- 0
  
  for (channel  in 1:n.channels) {
    beg.spike <- channel.start[channel]+2 #first spike for this channel
    if (beg.spike < channel.start[channel+1]) {
      ## We have some data.
      end.spike <- channel.start[channel+1] - 1 #last line

      ## convert data from msec to seconds.
      this.train <- as.numeric(data[beg.spike:end.spike,1]) / 1000
      if ( (length(this.train)>0) && !is.null(beg)) {
        rejs <- this.train < beg
        if (any(rejs))
          this.train <- this.train[!rejs]
      }
      ##browser()
      if ( (length(this.train)>0) && !is.null(end)) {
        rejs <- this.train > end
        if (any(rejs))
          this.train <- this.train[!rejs]
      }

      if ((length(this.train)>0)) {
        n <- n + 1
        spikes[[n]] <- this.train
        channel.ok[n] <- channels[channel]
      }

    }

  }

  channel.ok <- channel.ok[1:n]         #truncate to right length.
  ## TODO: could filter out spike times at this point, once all data has
  ## been read in.  Then can use jay.filter.for.max/min code as elsewhere.
  res <- list(spikes=spikes, channels=channel.ok) 
}





ncl.read.spikes <- function(filename, ids=NULL,
                            time.interval=1, beg=NULL, end=NULL) {

  ## Read in Ncl data set.  IDS is an optional vector of cell numbers
  ## that should be analysed -- the other channels are read in but
  ## then ignored.


  ## Tue 25 Jul 2006 -- I'm not sure if this is the cleanest way to read in the
  ## data files, as it probably can be read in using just
  ## read.csv(filename, sep='\t')
  ## but since this works, fine.
  

  dat <- mcd.data.to.array(filename, beg, end)
  spikes <- dat$spikes
  channels <- dat$channels
  names(spikes) <- channels

  spikes.range <- range(unlist(spikes))
  if (is.null(beg))  beg <-  spikes.range[1]
  if (is.null(end))  end <-  spikes.range[2]
  rec.time <- c(beg, end)

  
  ## Count the number of spikes per channel, and label them.
  nspikes <- sapply(spikes, length)

  ## TODO
  ## names(nspikes) <- channels 

  ## meanfiring rate is the number of spikes divided by the (time of
  ## last spike - time of first spike).  
  ##meanfiringrate <- nspikes/ ( sapply(spikes, max) - sapply(spikes, min))
  ## todo -- find min, max time.

  meanfiringrate <- nspikes/ ( end - beg)

  ## Parse the channel names to get the cell positions.
  layout <- make.sanger1.layout(substring(channels, 1,2))
  
  ## check that the spikes are monotonic.
  check.spikes.monotonic(spikes)


  rates <- make.spikes.to.frate(spikes, time.interval=time.interval,
                                beg=beg, end=end)
  
  ## See if we need to shift any units.  this affects only the
  ## visualisation of the units in the movies.  We assume that "shifted"
  ## positions are stored in the file with same name as data file
  ## except that the .txt is replaced with .sps.  Then each line of this
  ## file contains three numbers:
  ## c dx dy
  ## where c is the cell number to move, and dx,dy is the amount (in um)
  ## by which to move the cells.  If you edit the file, this function
  ## must be called again for the new values to be read in.
  ## The shifted positions are used only by the movie functions and
  ## by the function plot.shifted.jay.pos(s) [this shows all units].


  ## Tue 19 Dec 2006: this assumes filename ends in .txt; do not worry
  ## about this for now.
  
##   shift.filename <- sub("\\.txt$", ".sps", filename)
   unit.offsets <- NULL                  #default value.
##   if (FALSE && file.exists(shift.filename)) { #TODO -- why running?
##     updates <- scan(shift.filename)
##     ## must be 3 data points per line
##     stopifnot(length(updates)%%3 == 0)
##     updates <- matrix(updates, ncol=3, byrow=TRUE)
##     units <- updates[,1]
##     if (any(units> length(spikes))) {
##       stop(paste("some units not in recording...",
##                  paste(units[units>=length(spikes)],collapse=",")))
##     }
##     unit.offsets <- layout$pos*0               #initialise all elements to zero.
##     unit.offsets[units,] <- updates[,2:3]
##   }
  
  
  
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

  ncl.breaks = c(0, 150, 250, 350, 450, 550, 650, 1000, 2000)
  res$corr = corr.index(res, ncl.breaks)

  res

}
