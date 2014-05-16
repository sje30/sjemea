## R script to read in a nex file.
## 2014-03-14


##' Read in a Neuroxplorer data file.
##'
##' This code has been adapted from the matlab version of the script
##' provided by http://www.neuroexplorer.com/code.html
##' 
##' @param f File name
##' @param array.name by default this is the 8 by 8 array with 100um separation.
##' @param channel.regexp Regular expression for finding the channel numbers.
##' @return An s object.
##' @author Stephen Eglen
nex.read.spikes <- function(f, 
                            array.name='MCS_8x8_100um',
                            channel.regexp="_([0-9])([0-9])") {

  nexdata = readnex(f)

  names <- sapply(nexdata$neurons, function(x) x$name)

  spacing <- get.array.info(list(array=array.name))$layout$spacing
  row.cols <- parse.channel.names(names, regexp=channel.regexp,
                                  spacing=spacing)
  ## Check that all electrode names are valid.
  valid.name <- !is.na(row.cols[,1])

  ## Show what electrode names are ignored.
  if (any(!valid.name)) {
    cat("The following electrode names are removed\n")
    print(names[!valid.name])
  }
  ## Keep the information only for valid electrodes.
  names <- names[valid.name]
  row.cols <- row.cols[valid.name,]
  neurons <- nexdata$neurons[valid.name]
  spikes <- sapply(neurons, function(x) x$timestamps)
  nspikes <- sapply(spikes, length)

  epos <- row.cols
  data <- list()
  data$epos <- epos
  data$names <- names
  data$array <- array.name
  arrayinfo <- get.array.info(data)
  layout <- arrayinfo$layout

  names(spikes) <- names
  s <- construct.s(spikes, ids=NULL, time.interval=1,
                   beg=nexdata$tbeg, end=nexdata$tend,
                   corr.breaks=arrayinfo$corr.breaks,
                   layout=layout, filename=f)

  s

}


parse.channel.names <- function(names, regexp="_([0-9])([0-9])", spacing=100) {
  ## Extract the two numbers from names like "Ch_48a" --> (4, 8).
  ## spacing is our multiplier to convert the channel numbers into distances.
  ## NAMES is of length N; return matrix of size Nx2 with row and col information.
  ## If no row,col could be parsed, return NA in that row.

  m <- regexec(regexp, names)
  matches <- regmatches(names, m)
  rows <- sapply(matches, function(o) if (length(o)==0) { NA } else {as.numeric(o[[2]])})
  cols <- sapply(matches, function(o) if (length(o)==0) { NA } else {as.numeric(o[[3]])})

  m <- spacing * cbind(rows,cols)
  rownames(m) <- names
  m
}

readnex <- function(file) {
  ## Read in the entire NEX file.
  

  fid <-  file(file, 'rb')
  on.exit(close(fid))
  
  
  nexFile <- list()

  nexFile$magic = readBin(fid, "int", 1, endian='little')
  stopifnot(nexFile$magic == 827868494)


  nexFile$version = readBin(fid, "int", 1, endian='little')
  nexFile$comment = readChar(fid, 256, useBytes=TRUE)
  
  nexFile$freq = readBin(fid, 'double', 1, endian='little');
  nexFile$tbeg = readBin(fid, 'int', 1, endian='little') / nexFile$freq;
  nexFile$tend = readBin(fid, 'int', 1, endian='little') / nexFile$freq;

  nvar = readBin(fid, 'int', 1, endian='little');
  nexFile$nvar = nvar
  cat(sprintf("Reading in %d variables from %s\n", nvar, file))
  # skip location of next header fields and padding
  seek(fid, 260, 'current')

  ## Create an empty list.  All of them will be longer than needed
  ## so can be truncated at end.
  neurons <- vector("list", nvar)
  events <- vector("list", nvar)
  intervals <- vector("list", nvar)
  waveforms <- vector("list", nvar)
  conts <- vector("list", nvar)
  
  neuronCount = 0;
  eventCount = 0;
  intervalCount = 0;
  waveCount = 0;
  popCount = 0;
  contCount = 0;
  markerCount = 0;


  for (variableIndex in 1:nvar) {
    # read variable header
    type = readBin(fid, 'int', 1, endian='little');
    varVersion = readBin(fid, 'int',1, endian='little');
    name = readChar(fid, 64, useBytes=TRUE)
    # remove first zero and all characters after the first zero
    # TODO name(end+1) = 0;
    #name = name(1:min(find(name==0))-1);
    offset = readBin(fid, 'int', 1, endian='little');
    n = readBin(fid, 'int', 1, endian='little');
    wireNumber = readBin(fid, 'int', 1, endian='little');
    unitNumber = readBin(fid, 'int', 1, endian='little');
    gain = readBin(fid, 'int', 1, endian='little');
    filter = readBin(fid, 'int', 1, endian='little');
    xPos = readBin(fid, 'double', 1, endian='little');
    yPos = readBin(fid, 'double', 1, endian='little');
    WFrequency = readBin(fid, 'double', 1, endian='little'); # wf sampling fr.
    ADtoMV  = readBin(fid, 'double', 1, endian='little'); # coeff to convert from AD values to Millivolts.
    NPointsWave = readBin(fid, 'int', 1, endian='little'); # number of points in each wave
    NMarkers = readBin(fid, 'int', 1, endian='little'); # how many values are associated with each marker
    MarkerLength = readBin(fid, 'int', 1, endian='little'); # how many characters are in each marker value
    MVOfffset = readBin(fid, 'double', 1, endian='little'); # coeff to shift AD values in Millivolts: mv = raw*ADtoMV+MVOfffset
    filePosition = seek(fid)
    okay = 0;
    if (type == 0) {
      neuronCount = neuronCount+1; okay = 1

      neuron <- list()
      neuron$name = name
      neuron$varvarVersion = varVersion
      if (varVersion > 100) {
        neuron$wireNumber = wireNumber
        neuron$unitNumber = unitNumber
      } else {
        neuron$wireNumber = 0
        neuron$unitNumber = 0
      } 
      neuron$xPos = xPos
      neuron$yPos = yPos
      ## go to the variable data position and read timestamps
      seek(fid, offset, 'start');
      neuron$timestamps = readBin(fid, 'int', n, endian='little') / nexFile$freq

      neurons[[neuronCount]] <- neuron

    }

    if (type == 1) {
      eventCount = eventCount +1; okay = 1
      event = list(name=name)
      event$varVersion = varVersion;
      seek(fid, offset, 'start');
      event$timestamps = readBin(fid, 'int', n, endian='little') / nexFile$freq;
      events[[eventCount]] <- event
    }


    if (type == 2) {
      intervalCount = intervalCount +1; okay = 1
      interval = list(name=name)
      interval$varVersion = varVersion;
      seek(fid, offset, 'start');
      interval$intStarts = readBin(fid, 'int', n, endian='little') / nexFile$freq;
      interval$intEnds   = readBin(fid, 'int', n, endian='little') / nexFile$freq;
      intervals[[intervalCount]] <- interval
    }


    if (type == 3) {
      waveCount = waveCount + 1; okay = 1
      wave = list(name=name)
      wave$varVersion = varVersion;
      wave$NPointsWave = NPointsWave;
      wave$WFrequency = WFrequency;
      if (varVersion > 100) {
        wave$wireNumber = wireNumber;
        wave$unitNumber = unitNumber;
      } else {
        wave$wireNumber = 0;
        wave$unitNumber = 0;
      }
      wave$ADtoMV = ADtoMV
      if (nexFile$version > 104) {
        wave$MVOfffset = MVOfffset
      } else {
        wave$MVOfffset = 0
      }
      seek(fid, offset, 'start')
      wave$timestamps = readBin(fid, 'int', n, endian='little') / nexFile$freq
      wf = readBin(fid, 'int', NPointsWave*n, size=2, endian='little')
      wf = matrix(wf, nrow=NPointsWave, ncol=n)
      wave$waveforms = wf * ADtoMV + wave$MVOfffset;
      waveforms[[waveCount]] <- wave
    }

    if (type == 5) { ## continuous Variable
      contCount = contCount+1; okay = 1;
      cont = list(name = name)
      cont$varVersion = varVersion
      cont$ADtoMV = ADtoMV
      if (nexFile$version > 104) {
        cont$MVOfffset = MVOfffset
      } else {
        cont$MVOfffset = 0
      }
      cont$ADFrequency = WFrequency
      seek(fid, offset, 'start')
      cont$timestamps = readBin(fid, 'int', n, endian='little') / nexFile$freq
      cont$fragmentStarts = readBin(fid, 'int', n, endian='little') + 1
      cont$data = readBin(fid, 'int', size=2, n=NPointsWave, endian='little') * ADtoMV +
        cont$MVOfffset

      conts[[contCount]] <- cont
    }

    ## What follows is the code for the other types of variable.
    ## We can add this if/when we need it.
    
    ##     case 4 # population vector
    ##         popCount = popCount+1;
    ##         nexFile.popvectors{popCount,1}.name = name;
    ##         nexFile.popvectors{popCount,1}.varVersion = varVersion;
    ##         fseek(fid, offset, 'bof');
    ##         nexFile.popvectors{popCount,1}.weights = fread(fid, [n 1], 'double');
    
    
    ##     case 6 # marker
    ##         markerCount = markerCount+1;
    ##         nexFile.markers{markerCount,1}.name = name;
    ##         nexFile.markers{markerCount,1}.varVersion = varVersion;
    ##         fseek(fid, offset, 'bof');
    ##         nexFile.markers{markerCount,1}.timestamps = fread(fid, [n 1], 'int')./nexFile.freq;
    ##         for markerFieldIndex=1:NMarkers
    ##             markerName = fread(fid, 64, '*char')';
    ##             # remove first zero and all characters after the first zero
    ##             markerName(end+1) = 0;
    ##             markerName = markerName(1:min(find(markerName==0))-1);
    ##             nexFile.markers{markerCount,1}.values{markerFieldIndex,1}.name = markerName;
    ##             for markerValueIndex = 1:n
    ##                 markerValue = fread(fid, MarkerLength, '*char')';
    ##                 # remove first zero and all characters after the first zero
    ##                 markerValue(end+1) = 0;
    ##                 markerValue = markerValue(1:min(find(markerValue==0))-1);
    ##                 nexFile.markers{markerCount,1}.values{markerFieldIndex,1}.strings{markerValueIndex, 1} = markerValue;
    ##             end
    ##         end

    if (!okay) {
      stop("Other types not added yet", type)
    }
    

    ## return to file position that was after reading the variable header
    seek(fid, filePosition, 'start');
    dummy = readBin(fid, 'raw', 60, endian='little');
  }

  
  ## Tidy up and return information

  ## check that the number of variables add up.
  stopifnot( nvar == ( neuronCount + eventCount + intervalCount +
                      waveCount + popCount + contCount + markerCount))

  nexFile$neurons <- neurons[1:neuronCount]
  nexFile$events <- events[1:eventCount]
  nexFile$intervals <- intervals[1:intervalCount]
  nexFile$waveforms <- waveforms[1:waveCount]
  nexFile$conts <- conts[1:contCount]
  nexFile

}

