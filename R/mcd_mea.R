## mcd_mea.R --- import data from MC_DataTool
## Author: Stephen J Eglen, P Jarzebowski
## Copyright: GPL


##' Read txt data file with spike times exported from MCD file with MC_DataTool.
##' 
##' Returns the spike trains and the channel names in the imported recording.
##' Spikes outside of the time window [BEG, END] (if non-null) are ignored.
##' Blank lines, pre and post-spike voltages will be ignored.
##' 
##' @param file input file
##' @param beg time in seconds from which data should be read in, NULL if should
##' read from the beginning
##' @param end time in seconds until which data should be read in, NULL if 
##' should read until the end of the file
##' @param pre.spike.ms time in ms before a spike for which the input file
##' has extra lines with voltage readout, 0 if no pre-spike times in the file
##' @param post.spike.ms time in ms after a spike for which the input file
##' has extra lines with voltage readout, 0 if no post-spike times in the file
mcd.data.to.array <- function(file, beg=NULL, end=NULL, pre.spike.ms=0, 
                              post.spike.ms=0) {
  
  data <- read.table(file, as.is=TRUE, skip=2, fill=TRUE)
  

  ## Row number of the start of each spike train.
  channel.start <- which(data[,2] =='Spikes')
  n.channels <- length(channel.start)
  
  ## For the last channel, we add where the end of the file is.
  ## This is so that our check for no spikes will also work for the last
  ## channel in the data set.
  channel.start <- c(channel.start, nrow(data) + 1)

  ## Name of each channel.
  channels <- data[channel.start,4]
  channel.ok <- rep(0, n.channels)
  spikes <- list(); n <- 0
  
  for (channel in 1:n.channels) {
    beg.spike <- channel.start[channel]+2 #first spike for this channel
    if (beg.spike < channel.start[channel+1]) {
      ## We have some data.
      end.spike <- channel.start[channel+1] - 1 #last line

      this.times <- as.numeric(data[beg.spike:end.spike,1])
      this.train <- spike.train.from.cut.offs(this.times,
                                              pre.spike.ms, 
                                              post.spike.ms)
      
      ## convert data from msec to seconds
      this.train <- round(this.train / 1000, 3)
      ## make sure no duplicates after lowering the precision
      this.train <- unique(this.train)
      if ( (length(this.train)>0) && !is.null(beg)) {
        rejs <- this.train < beg
        if (any(rejs))
          this.train <- this.train[!rejs]
      }

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


##' Returns a spike train read from spike cut offs.
##' 
##' @param timestamps is a vector of recorded times which is preceeded by pre- and 
##' postspike entries
##' @param pre.spike.ms time in ms before a spike for which the times were
##' recorded
##' @param post.spike.ms time in ms after a spike for which the times were
##' recorded
spike.train.from.cut.offs <- function(timestamps, pre.spike.ms, post.spike.ms) {
  this.train <- rep(NA, length(timestamps))
  spike.count <- 0
  current.pre.spike <- -pre.spike.ms - post.spike.ms - 1
  time.step <- 1
  if (length(timestamps) > 1) {
    time.step <- timestamps[2] - timestamps[1]
  }
  for (i in 1:length(timestamps)) {
    if (current.pre.spike + pre.spike.ms + post.spike.ms < timestamps[i]) {
      current.pre.spike = timestamps[i]
    }
    if (current.pre.spike + pre.spike.ms <= timestamps[i] && 
        timestamps[i] < current.pre.spike + pre.spike.ms + time.step) {
      spike.count <- spike.count + 1
      this.train[[spike.count]] <- timestamps[i]
    }
  }
  
  this.train[1:spike.count]
}


##' Creates s object from the NCL dataset input file.
##'
##' @seealso \code{\link{mcd.read.spikes}}
ncl.read.spikes <- function(filename, ids=NULL, time.interval=1, beg=NULL, 
                            end=NULL) {
  mcd.read.spikes(filename, ids, time.interval, beg, end)
}


##' Creates s object from txt data file exported from MCD file with MC_DataTool.
##' 
##' Returns the spike trains and the channel names in the imported recording.
##' Spikes outside of the time window [BEG, END] (if non-null) are ignored.
##' Blank lines, pre and post-spike voltages will be ignored.
##' 
##' @param filename input file
##' @param ids IDS optional vector of cell numbershat should be analysed -- the 
##' other channels are read in but then ignored.
##' @param time.interval time bin used for calculation of the firing rate
##' @param beg time in seconds from which data should be read in, NULL if should
##' read from the beginning
##' @param end time in seconds until which data should be read in, NULL if 
##' should read until the end of the file
##' @param pre.spike.ms time in ms before a spike for which the input file
##' has extra lines with voltage readout. 
##' @param post.spike.ms time in ms after a spike for which the input file
##' has extra lines with voltage readout. 
mcd.read.spikes <- function(filename, ids=NULL, time.interval=1, beg=NULL, 
                            end=NULL, pre.spike.ms=0, post.spike.ms=0) {

  ## Tue 25 Jul 2006 -- I'm not sure if this is the cleanest way to read in the
  ## data files, as it probably can be read in using just
  ## read.csv(filename, sep='\t')
  ## but since this works, fine.
  

  dat <- mcd.data.to.array(filename, beg, end, pre.spike.ms, post.spike.ms)
  spikes <- dat$spikes
  channels <- as.character(dat$channels)
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
  layout$array <- "MCS_8x8_200um"
  
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
