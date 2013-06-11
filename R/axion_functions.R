## General functions useful for processing the Axion data.
## 2013-01-29

## This matrix stores information about the layout of the electrodes into
## wells.
## TODO: this should be redundant.
## axion.layouts <- matrix(c(12, 3, 4, 8, 8,
##                           48, 6, 8, 4, 4),
##                         nrow=2, byrow=TRUE)
## colnames(axion.layouts) <- c("well", "max.well.r", "max.well.c", "max.elec.r",
##                              "max.elec.c")

## This variable stores all the information related to the wells; typically this
## is accessed through the plateinfo(arrayname) function.
.plateinfo <- list("Axion 48 well"=list(
                     n.well=48,
                     wells=paste( rep(LETTERS[6:1], each=8), rep(1:8,6), sep=''),
                     n.well.r=6,
                     n.well.c=8,
                     layout=c(8,6),
                     n.elec.r=4,
                     n.elec.c=4),
                 "Axion 12 well"=list(
                     n.well=12,
                     wells=paste( rep(LETTERS[3:1], each=4), rep(1:4,3), sep=''),
                     n.well.r=3,
                     n.well.c=4,
                     layout=c(4,3),
                     n.elec.r=8,
                     n.elec.c=8))

                     
plateinfo <- function(arrayname) {
    ## Return useful information related to arrayname
    ## 
    ## plateinfo("Axion 12 well")
    res <- .plateinfo[[arrayname]]
    if (is.null(res)) {
        stop("arrayname not recognised:", arrayname)
    } else {
        res
    }
}
    
## * Scripts to convert data into HDF5.
axion.elec.name.to.xy <- function(name, plateinfo) {
  ## Convert electrode name to  (x,y) position.
  ## plateinfo stores all the information about the plates.
  ## and hence the well layout of the plate.

  max.well.row <-  plateinfo$n.well.r
  max.well.col <-  plateinfo$n.well.c
  max.elec.row <-  plateinfo$n.elec.r
  max.elec.col <-  plateinfo$n.elec.c
  
  
  well.r <- max.well.row - match(substring(name, 1,1), LETTERS)
  well.c <- as.integer(substring(name, 2,2)) - 1
  elec.r <- as.integer(substring(name, 5,5)) - 1
  elec.c <- as.integer(substring(name, 4,4)) - 1

  gap <- 1
  spacing <- 200                        #electrode spacing.
  well.wid <- (max.elec.col+gap)*spacing
  well.ht  <- (max.elec.row+gap)*spacing

  x <- (well.c*well.wid) + (elec.c*spacing)
  y <- (well.r*well.ht)  + (elec.r*spacing)

  cbind(x,y)

}

map2list <- function(file) {
  ## read in the data file, and then make the list of spike times.
  data <- scan(file, what=character(), sep='\t')
  max.channels <- 1997
  channels <- grep('^Elect_', data[1:max.channels])
  nchannels <- length(channels)
  stopifnot(nchannels != max.channels)  #i.e. max.channels not large enough.
  names <- substring(data[channels], 7, 11)
  well <- as.factor(substring(names, 1, 2))
  data3 <- data[-c(1:nchannels, length(data))]
  longest.spike <- length(data3)/nchannels #should be integer
  spikes.matrix <- matrix(data3, nrow=longest.spike, byrow=T)

  spikes <- apply(spikes.matrix, 2, function(x) sjemea::jay.filter.for.na(as.numeric(x)))
  names(spikes) <- names

  spikes
}

axion.spikesum <- function(spikes) {
  ## Generate a simple summary of the spikes list.
  len <- length(spikes)
  all.range <- sapply(spikes, range)
  nspikes <- sum(sapply(spikes, length))
  min <- min(all.range[1,])
  max <- max(all.range[2,])
  str <- sprintf("summary: %d electrodes %d spikes, min %.4f max %.4f",
                 len, nspikes, min, max)
  str
}

axion.spikesum2 <- function(spikes) {
  ## Generate a simple summary of the spikes list.
  ## This version returns a vector, rather than a string.  This is more
  ## useful for building a data frame of values.
  len <- length(spikes)
  all.range <- sapply(spikes, range)
  nspikes <- sum(sapply(spikes, length))
  min <- min(all.range[1,])
  max <- max(all.range[2,])
  str <- sprintf("summary: %d electrodes %d spikes, min %.4f max %.4f",
                 len, nspikes, min, max)
  ##str
  c(nelectrodes=len, nspikes=nspikes, time.min=min, time.max=max)
}

axion.spikestodf <- function(spikes) {
  ## Convert a list of spikes to a 2-column  (elec, time) data frame.
  names <- names(spikes)
  names(spikes) <- NULL
  nspikes <- sapply(spikes, length)
  data.frame(elec=rep(names, times=nspikes), time=unlist(spikes))
}


old.file2Rdata <- function(file) {
  ## This seems to take as long to read in as the raw file!!!
  data <- scan(file, what=character(), sep='\t')
  Rdatafile <- sprintf("%s.Rdata", file)
  save(data, file=Rdatafile)
  Rdatafile
}

axion.map.to.h5 <- function(key, make.venn=TRUE, debug=TRUE) {
  ## Return the name of the file.
  h5file <- sprintf("%s/epa%s.h5", h5.dir, key)
  wildcard <- sprintf('^%s.*mapTimestamps(.xz)?$', key)
  f <- list.files(path=data.dir, pattern=wildcard, full.names=TRUE)
  ## Allow for the files to be compressed, as R can transparently uncompress
  ## the files.  xz compression is pretty good here.  use the command "xz" to
  ## compress a file, or "unxz" to uncompress it.
  if (debug) {
    print(f)
  }
  spikes.sep <- lapply(f, map2list)       #can take up to a minute or so.

  short.filenames <- gsub(".mapTimestamps", "", basenamepy(f)$base)
  summary.table <- t(sapply(spikes.sep, axion.spikesum2))
  rownames(summary.table) <- short.filenames

  ## Put all the 2-dataframes together into one big dataframe
  ma <- do.call("rbind", lapply(spikes.sep, axion.spikestodf))
  
  ## Here we have our combined spike train.
  s2 <- split(ma$time, ma$elec)

  numelec <- length(s2)
  total.spikes <- sum(sapply(s2, length))
  ##time.range <- range(unlist(s2))         #quite slow.
  ## find total range
  time.ranges <- sapply(s2, range)
  time.min <- min(time.ranges[1,])
  time.max <- max(time.ranges[2,])
  cat(printf("Total number of spikes: %d\n", total.spikes))
  cat(printf("Unique number of electrodes: %d\n", numelec))
  cat(printf("Time range [%.3f %.3f] (seconds)\n", time.min, time.max))
  print(summary.table)

  ## Now save the HDF5 file.
  map.to.h5(s2, h5file)
  if (debug) {
    ## keep a record of the summary.table
    d <- as.data.frame(summary.table)
    d2 <- data.frame(file=rownames(summary.table), d, stringsAsFactors=FALSE)
    h5write(d2, path.expand(h5file), "summary.table")
  }
  if (make.venn && is.element(length(f), c(2,3))) {
    ## Try to make a Venn diagram of the resulting files, showing
    ## how many electrodes were in each recording.
    ## We make them only if there are two or three files in the recording.
    ## They will be stored in the current directory.
    
    require(VennDiagram)
    require(grid)
    if (length(f) == 2) {
      elec1 <- names(spikes.sep[[1]])
      elec2 <- names(spikes.sep[[2]])

      setdiff(elec1, elec2)
      setdiff(elec2, elec1)
      cat.pos <- c(-10, 10)
    }

    if (length(f) == 3) {
      elec1 <- names(spikes.sep[[1]])
      elec2 <- names(spikes.sep[[2]])
      elec3 <- names(spikes.sep[[3]])

      setdiff(elec1, elec2)
      setdiff(elec1, elec3)
      setdiff(elec2, elec1)
      setdiff(elec2, elec3)
      cat.pos <- c(-10, 10, 180)
    }

    x = lapply(spikes.sep, names)
    sapply(x, length)
    ##names <- paste('set', 1:(length(x)))
    names(x) <- short.filenames

    pdf(file=sprintf('%s_venn.pdf', key))
    p <- venn.diagram(x = x, euler.d = FALSE, scaled=FALSE, 
                      cat.pos = cat.pos, filename = NULL, reverse=FALSE)
    grid.newpage()
    grid.draw(p)
    dev.off()
  }

  ## Return the name of the hdf file created.
  h5file
}

map.to.h5 <- function(spikes, h5file) {
  ## Given a list of spikes, save the HDF5 file.

  h5file <- path.expand(h5file)
  if (file.exists(h5file))
    unlink(h5file)

  nspikes <- sapply(spikes, length)
  channels <- names(spikes)
  wells <- axion.guess.well.number(channels)
  array <- sprintf("Axion %d well", wells)
  plateinfo <- plateinfo(array)
  epos <- axion.elec.name.to.xy(channels, plateinfo)
  h5createFile(h5file)

  ## Let's compress the spike train by first creating chunks and
  ## compression options.  TODO, fix, it, not yet working for me.
  ## 2013-01-24
  sum.spikes <- sum(nspikes)
  ##h5createDataset(h5file, "/spikes", dims=sum.spikes, chunk=1e5, level=7)
  h5write(unlist(spikes), h5file, "/spikes")
  h5write(nspikes, h5file, "/sCount")
  h5write(epos, h5file, "/epos")
  h5write(channels, h5file, "/names")
  h5write(array, h5file, "/array")
  print(h5ls(h5file))

  
}

axion.guess.well.number <- function(channels) {
  ## Given the channel names, guess the number of wells on the plate.
  ## This works on the logic that certain electrode names will only be
  ## found on certain plates. e.g. the electrode name "D6_33" can only appear
  ## on a well with 48 arrays.
  ##
  ## axion.guess.well.number("D3_33")  ## should be 48.
  ## axion.guess.well.number("B3_53")  ## should be 12
  ## axion.guess.well.number("A2_11") ## this is ambiguous.
  
  well.r <- match(substring(channels, 1,1), LETTERS)
  well.c <- as.integer(substring(channels, 2,2))
  elec.r <- as.integer(substring(channels, 5,5))
  elec.c <- as.integer(substring(channels, 4,4))

  max.well.r <- max(well.r)
  max.well.c <- max(well.c)

  max.elec.r <- max(elec.r)
  max.elec.c <- max(elec.c)

  nplates <- length(.plateinfo)
  well <- 0
  for (i in 1:nplates) {
    plateinfo <- .plateinfo[[i]]
    if (max.well.r <= plateinfo$n.well.r &&
        max.well.c <= plateinfo$n.well.c &&
        max.elec.r <= plateinfo$n.elec.r &&
        max.elec.c <= plateinfo$n.elec.c) {
      well <- plateinfo$n.well
      break;
    }
  }
  if (well == 0) {
    stop("Cannot guess number of wells on plate.")
  }

  well

}


axion.electrodes.on.well <- function(well, electrodes) {
  ## Return names of electrodes that are on well WELL.
  matches <- grep(well, electrodes)
  electrodes[matches]
}


axion.elec2well <- function(elec) {
  ## Extract well name from ELECtrode name.
  substring(elec, 1, 2)
}


## Burst analysis summaries.
axion.filter.bursts <- function(allb, min.bursts=0, well) {
  ## Apply a couple of filter to list of burst info.
  nbursts <- sapply(allb, function(b) {
    if(is.na(b[1])) 0 else nrow(b)})
      
  if (!missing(well)) {
    well.idx <-  axion.elec2well(names(allb))==well
  } else {
    well.idx <- rep(TRUE, length(allb))
  }

  index <- which ( (nbursts >=min.bursts) & well.idx)
  allb[index]
}

axion.burst.summary.for.well <- function(well) {
  b1 <- axion.filter.bursts(allb, min.bursts=10, well=well)

  durns <- lapply(b1, function(b) b[,"durn"])
  stripchart(durns, vertical=TRUE, ylab='burst duration (s)', 
             method='jitter', pch=20)

  nbursts <- lapply(b1, function(b) nrow(b))
  stripchart(nbursts, vertical=TRUE, ylab='Number of bursts', 
             method='jitter', pch=20)

  ibi <- lapply(b1, function(b) b[,"IBI"])
  stripchart(ibi, vertical=TRUE, ylab='Interburst Interval (s)', 
             method='jitter', pch=20)

  mean.isi <- lapply(b1, function(b) b[,"mean.isis"])
  stripchart(mean.isi, vertical=TRUE, ylab='mean ISI within burst (s)', 
             method='jitter', pch=20)
}
