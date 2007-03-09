## Code created for Sanger.
## Fri 19 Jan 2007


######################################################################
## Code for Sanger MEA analysis,
######################################################################

sanger.init <- function() {
  ## Run for initialisation of analysis of Sanger data.

  windows <- .Platform$OS.type !='unix'

  if (windows) {
    ## Set up on Windows machine.
    ## Scripts are stored in meadev/scripts.
    setwd("c:/meadev/scripts")
    assign("mea.data.dir", "c:/meadev/data/", env = .GlobalEnv)
    assign("mea.table.dir", "c:/meadev/tables/", env = .GlobalEnv)
    assign("mea.op.dir", "c:/meadev/op/", env = .GlobalEnv)
  } else {
    ## Set up for linux.
    assign("mea.data.dir",   "~/proj/sangermea/data/", env = .GlobalEnv)
    assign("mea.table.dir",  "~/proj/sangermea/tables/", env = .GlobalEnv)
    assign("mea.op.dir",     "~/proj/sangermea/op/", env = .GlobalEnv)
  }

  ## Create the cache of datafiles.
  assign("mea.data.files",  make.meafile.cache(mea.data.dir),
         env  = .GlobalEnv)
}

mea.op <- function(f) {
  ## Create a new output file.
  sprintf("%s%s", mea.op.dir, f)
}

make.sanger1.layout <- function(positions) {
  ## make the layout for SANGER MEA


  xlim <- ylim <- c(50, 1700)
  spacing <- 200

  cols <- as.integer(substring(positions, 1,1)) * spacing
  rows <- (9-as.integer(substring(positions, 2,2))) * spacing
  pos <- cbind(cols, rows)
  
  rownames(pos) <- positions
  
  layout <- list(xlim=xlim, ylim=ylim, spacing=spacing,
                 pos=pos)

  class(layout) <- "mealayout"

  layout

}

sanger.read.spikes <- function(filename, ids=NULL,
                               time.interval=1,
                               beg=NULL, end=NULL,
                               min.rate=0) {

  ## Read in Sanger data set.  
  ## IDS (IGNORE FOR NOW)is
  ## an optional vector of cell numbers that should be analysed -- the
  ## other channels are read in but then ignored.

  ## MIN.RATE (when >0) is the mininum firing rate (Hz) that a channel
  ## must have for it to be analysed.  e.g. a sensible threshold would
  ## be 1/60 to indicate it must on average make one spike/minute to
  ## be analysed.
  

  ## Tue 25 Jul 2006 -- I'm not sure if this is the cleanest way to read in the
  ## data files, as it probably can be read in using just
  ## read.csv(filename, sep='\t')
  ## but since this works, fine.
  
  fp <- file(filename, open="r")

  header <- readLines(fp, n=1)
  header.split <- strsplit(header, '\t')[[1]]
  n.cols <- length(header.split)
  stopifnot( header.split[n.cols] == 'Sweep_Stop')
  ## We have two extra args to worry about at the end.
  channels <- header.split[1:(n.cols-2)]
  num.channels <- length(channels)

  
  ## The first line of data is special in that it has the Sweep_Start
  ## and Sweep_stop times.
  line1 <- scan(fp, "", n=n.cols, sep='\t', quiet=TRUE)
  ## Throw away last two numbers on line 1, as they are the
  ## sweep_start and sweep_stop times.
  sweep.start <- as.numeric(line1[n.cols-1])
  sweep.stop  <- as.numeric(line1[n.cols])
  ##line1 <- line1[1:num.channels]
  

  ## Now read in remaining times.
  rest <- scan(fp, sep='\t', quiet=TRUE)
  close(fp)
  ## NA in the files indicates no more spikes for that channel.
  times <- c(as.double(line1), rest)
  stopifnot( (length(times) %% (num.channels+2)) == 0)

  ## Reformat to a matrix.
  ntimes <- length(times)
  dim(times) <- c(n.cols, ntimes/n.cols)

  ## Remove the last two rows, as they just the sweep_start and
  ## sweep_stop times.  drop=F is needed so that if there is one
  ## channel, the resulting times will remain an array.
  times <- times[1:num.channels,,drop=FALSE]

  ## Now convert each row of TIMES into a spike train, removing any
  ## NAs from the end of each row.  This used to be an apply() command
  ## but is messy when either you have just one spike train, or
  ## perhaps two spike trains with same number of spikes in each
  ## train.

  spikes <- list()
  for (i in 1:num.channels) {
    spikes[[i]] <- jay.filter.for.na(times[i,])
  }


  if (!is.null(end)) {
    spikes <- lapply(spikes, jay.filter.for.max, max=end)
  } else {
    end <- sweep.stop
  }

  if (!is.null(beg)) {
    spikes <- lapply(spikes, jay.filter.for.min, min=beg)
  } else {
    beg <- sweep.start
  }

  if (min.rate > 0 ) {
    
    ## Check for inactive channels -- those with a mean firing rate
    ## below some average rate.

    ## This catches the odd situation when a channel has no spikes on
    ## it -- this can happen when a duration (beg, end) is given where
    ## no spikes occur on that channel.
    
    
    nspikes <- sapply(spikes,length)
    durn <- sweep.stop - sweep.start
    rates <- nspikes/durn
    inactive <- which(rates < min.rate)
    if (any(inactive)) {
      cat(paste("Removing spikes with low firing rates: ",
                paste(inactive, collapse=' '), "\n"
                ))
      spikes = spikes[-inactive]
      channels = channels[-inactive]
    }
    
    
  }


  
  if (!is.null(ids) ) {
    if (any(ids>length(spikes)))
      stop(paste("some ids not in this data set:",
                 paste(ids[ids>length(spikes)],collapse=" ")))
    
    spikes <- spikes[ids];
    channels <- channels[ids];
  }

  ## Count the number of spikes per channel, and label them.
  nspikes <- sapply(spikes, length)
  names(nspikes) <- channels

  ## meanfiring rate is the number of spikes divided by the (time of
  ## last spike - time of first spike).  
  meanfiringrate <- nspikes/ ( end - beg)

  ## Parse the channel names to get the cell positions.
  layout <- make.sanger1.layout(substring(channels, 4, 5))

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
  
  shift.filename <- sub("\\.txt$", ".sps", filename)
  unit.offsets <- NULL                  #default value.
  if (FALSE && file.exists(shift.filename)) { #TODO -- why running?
    updates <- scan(shift.filename)
    ## must be 3 data points per line
    stopifnot(length(updates)%%3 == 0)
    updates <- matrix(updates, ncol=3, byrow=TRUE)
    units <- updates[,1]
    if (any(units> length(spikes))) {
      stop(paste("some units not in recording...",
                 paste(units[units>=length(spikes)],collapse=",")))
    }
    unit.offsets <- pos*0               #initialise all elements to zero.
    unit.offsets[units,] <- updates[,2:3]
  }

  ## Compute CV of ISI.
  mean.isi = sapply(spikes, function(s) { mean(isi(s))})
  
  cv.isi = sapply(spikes, cv.isi)
  
  res <- list( channels=channels,
              totalspikes=sum(nspikes),
              spikes=spikes, nspikes=nspikes, NCells=length(spikes),
              meanfiringrate=meanfiringrate,
              file=filename,
              ##pos=pos,
              layout=layout,
              rates=rates,
              unit.offsets=unit.offsets,
              rec.time=c(beg, end),
              mean.isi=mean.isi,
              cv.isi = sapply(spikes, cv.isi)
              )
  class(res) <- "mm.s"

  ## Compute the correlation index.
  distance.breaks = c(0, 150, 250, 350, 450, 550, 650, 1000, 2000)
  res$corr = corr.index(res, distance.breaks)

  res

}


make.meafile.cache <- function(dir) {
  ## Remake the file cache.
  ## Search recursively through DIR to find all filenames.
  files <- dir(dir, recursive=TRUE, full.names=TRUE)
  mea.key <- basename(files)

  res <- cbind(mea.key, files)
  res
}


meafile <- function(file) {
  ## Use the file cache to find where the file is stored.
  ## This saves us having to always use the fullpath to a file.
  row <- which(file==mea.data.files[,1])
  row.n <- length(row)

  ## Do some safety checks.
  if ( row.n ==1) {
    ## return the full filename.
    mea.data.files[row,2]
  } else {
    if (row.n==0) {
      stop(sprintf("Cannot find file %s in mea.data.files\n", file))
    } else {
      stop(sprintf("File %s has %d entries in mea.data.files\n",
                   file, row.n))
    }
  }
}


meatable <- function(file) {
  ## Find the condition table.
  ## Uses global MEA.TABLE.DIR.
  file <- paste(mea.table.dir, file,sep='')
  if (!file.exists(file))
    stop(file, " not found")
  file
}

read.cond.tab <- function(file) {
  dat <- read.csv(file, as.is=T)


  ## Handle some sanger specific stuff:

  ## "Age (DIV)" is quite cumbersome, so shorten it to "Age"
  long.age <- pmatch("Age..DIV.", names(dat))
  if (length(long.age)==1)
    names(dat)[long.age] = "Age"

  ## strip trailing empty lines.
  empty.lines <- which(dat[,1] == "")
  if ( any(empty.lines) )
    dat <- dat[-empty.lines,]

  ## Remove any rows that should be ignored.
  ignore <- which(dat$Ignore == 1)
  if (any (ignore))

    dat <- dat[-ignore,]
  dat
}

