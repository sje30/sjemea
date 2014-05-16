## Code for reading in data from Alan Litke's array.

litke.spiketimes <- function(file) {
  ## read in the data, and return spike times and array position.
  require(R.matlab)
  z = readMatRda(file)
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
  if (any(end)) ## For longest spike train, there will be no spikes to remove
    spikes = spikes[-end]
  spikes / sample.rate
}

oneoff.make.litke1.layout <- function(file=NULL) {
  ## One-off function needed to create the litke1.txt layout file..
  ## make.litke1.layout('Electrode_Map.mat')
  ## How to hide from external world?  Put dot at start?  or do not export?
  if (is.null(file)) {
    file = system.file("extdata", "Electrode_Map.mat", package='sjemea')
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
  name <- "litke_hex_60um"
  array.info <- get.array.info(list(array=name))
  ##xlim <- c(-1000, 1000)                #duplicated in ./arrays.R
  ##ylim <- c(-500, 500)
  ##spacing <- 60
  xlim <- array.info$layout$xlim
  ylim <- array.info$layout$ylim
  spacing <- array.info$layout$spacing

  ##litke1layout <- NULL                  # silence the checker.
  data(litke1layout)
  
  columns <- match(positions, litke1layout$electrode)
  ## TODO what does match return if no matching electrode found?
  ##
  pos <- cbind(litke1layout$x[columns], litke1layout$y[columns],
               positions)
  
  rownames(pos) <- names
  colnames(pos) <- c("x", "y", "electrode.num")
  
  layout <- list(xlim=xlim, ylim=ylim, spacing=spacing,
                 pos=pos, array=name)

  class(layout) <- "mealayout"

  layout

}

show.litke.layout <- function() {
  ## Show the layout of the Litke array.
  litke1layout <- NULL; rm(litke1layout) #silence R CMD CHECK
  data(litke1layout)
  plot(litke1layout$x, litke1layout$y, asp=1,
       type='n', xlab='', ylab='')
  text(litke1layout$x, litke1layout$y, litke1layout$electrode, cex=0.7)
  title('Litke1 layout')
}





##' Read in the .mat file containing MEA data from Alan Litke.
##' 
##' Read in the matlab file and process as a regular MEA.
##' 
##' 
##' @aliases litke.read.spikes iit.read.spikes aps.read.spikes
##' @param filename Name of the text file to be read in.
##' @param ids Optional vector of cell numbers that can be analysed, rather
##' than analysing all electrodes in the recording.  Warning: Not implemented
##' in all readers.
##' @param time.interval Bin width (in seconds) for estimating firing rate.
##' Defaults to 1 second.
##' @param beg Optional start time.
##' @param end Optional end time.
##' @param corr.method Which correlation method to use?  Only 'ci' supported
##' for the Likte reader.
##' @return Return the data structure 's'.
##' @section METHOD: The data provided by Alan Litke is stored as a matlab
##' file.  There should be only one item in the file, which is the matrix
##' containing the spikes.
##' 
##' The matrix contains one row per electrode.  The first item in the row is
##' the cell number; the second item is the electrode that recorded the cell.
##' The remaining numbers are the spike times (which must be divided by 20000
##' to get time in seconds).  All spike trains are zero-padded to be the same
##' length.
##' 
##' As the matlab format is used, the R.matlab and Rcompression packages are
##' required.  These can be installed from R using:
##' 
##' install.packages(c('R.utils', 'R.matlab')) install.packages('Rcompression',
##' repos = "http://www.omegahat.org/R")
##' 
##' (The IIT data also use a matlab format, and so the above installations are
##' also needed for IIT.)
##' 
##' For the APS function, the data are stored in either HDF5 format or in an
##' RData file, converted from the .matlab equivalent (the matlab files take
##' far too long to read in, as they have long cell array structures).
##' 
##' The location of each electrode is given in the data file `litke1layout',
##' see examples below.
##' @seealso \code{\link{jay.read.spikes}}
##' @references This reader used for the Stafford et al. (2009) data in Neuron.
##' @keywords math
##' @examples
##' 
##' ## Show the layout of the array.
##' show.litke.layout()
##' 
##' ## The location of the file storing the electrode positions
##' system.file('data', 'litke1layout.txt', package='sjemea')
##' 
##' @export litke.read.spikes
litke.read.spikes <- function(filename, ids=NULL,
                              time.interval=1, beg=NULL, end=NULL,
                              corr.method="ci") {
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

  if (corr.method == "ci") {
    litke.breaks = c(0, 0.001, seq(from=50, to=2050, by=50)) #duplicate ./arrays.R
    
    res$corr = corr.index(res, litke.breaks)
  }

  res

}


######################################################################
## Convert matlab files to equivalent RDa files hidden under a "rda" subfolder.
## This allows us to read in matlab files quickly once they have been converted.
######################################################################


mat.to.rda <- function(f, check=TRUE, verbose=TRUE) {
    mat2rda = readMat(f)
    f2 = mat.rda.name(f)
    save(mat2rda, file=f2)
    if (check) {
        z = new.env()
        load(f2, z)
        stopifnot(all.equal(mat2rda, get("mat2rda", z)))
    }
    f2
}

mat.rda.name <- function(f) {
  ## Given the matlab file f, return equivalent rda filename
  parts = basenamepy(f)
  sprintf("%s/rda/%s.rda", parts$dir, parts$base)
}

mat.to.rda.dir <- function(dir) {
  ## Process whole directory of mat files.
  ## Convert each one to a rda file.

  ## Create the rda directory if it doesn't exist.
  rda.dir = sprintf("%s/rda", dir)
  if (!file.exists(rda.dir))
    dir.create(rda.dir)

  files = list.files(dir, pattern = "mat$", full.names=TRUE)
  lapply(files, mat.to.rda)
  
}


readMatRda <- function(file, verbose=TRUE) {
  ## wrapper around readMat; check first if the .rda file has been made.

  f2 = mat.rda.name(file)
  if (file.exists(f2)) {
    if (verbose)
      cat(sprintf("Reading %s as rda file\n", f2))
    
    env = new.env()
    load(f2, env)
    z = get("mat2rda", env)
  } else {
    z = readMat(file)
  }

  z
}
