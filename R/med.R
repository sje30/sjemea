## Read in data from the MED64 array.
## This reader was started early 2011, but no data yet available so
## have not pursued.
## ~/proj/carmen/med64/read_med.R

## model on the sql reader.

med64.spike.reader <- function(file, sep=300) {
  ## Read in the MED64 spike trains.
  dat <- read.csv(file,skip=4)

  ## time is in msec, so convert to s
  spikes <- split(dat$within_session_time_ms/1000, dat$channel)

  pos = cbind(x=res[,"x"], y=res[,"y"], electrode.num=res[,"num"])
  layout <- list(xlim=c(0,8)*sep,
                 ylim=c(0,8)*sep,
                 spacing=sep,
                 pos=pos)
  distbreaks <- NULL
  res <- list(layout=layout, spikes=spikes, distbreaks=distbreaks)

}

## channel no and timestamp are really needed. Just to be awkward, the
## MED 64 doesn't have a 200 micron spaced array, they are 150 or 300,
## I'm fairly sure this was a 300. The numbering is a simple
## 1,2,3...64 reading across, and then down (corner electrodes are not
## missing as in MCS, so rows begin 1, 9, 17 etc and 64 is the bottom
## right hand corner.

