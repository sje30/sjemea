## Reader/writer for sql data.



make.sql.file <- function(s, outputdb) {

  require(RSQLite)
  ## Create a database file from an MEA object.

  ## Remove duplicate rows from the data frame; this is where multiple
  ## cells are detected on the same electrode.
  p1 = s$layout$pos
  dups = duplicated(p1)
  p2 = p1[!dups,]

  electrode = data.frame(
    num=as.integer(p2[,"electrode.num"]),
    x=as.integer(p2[,"x"]),
    y=as.integer(p2[,"y"]))
  
  rownames(electrode) <- NULL


  neuron = data.frame(
    id=1:nrow(s$layout$pos),
    name=as.character(rownames(s$layout$pos)),
    electrode=as.integer(s$layout$pos[,"electrode.num"]))
  rownames(neuron) <- NULL


  geometry <- with(s$layout, 
                   data.frame(xlo=xlim[1], xhi=xlim[2],
                              ylo=ylim[1], yhi=ylim[2],
                              spacing=spacing))

  all.spiketimes <- unlist(s$spikes)
  all.spikeids <- rep(1:s$NCells, times=s$nspikes)
  spikes=data.frame(t=all.spiketimes,neuron=all.spikeids)

  distbreaks = data.frame(breaks=s$corr$distance.breaks)

  analysis = data.frame(id=as.integer(0), comments='raw data')


  m = dbDriver("SQLite")

  con = dbConnect(m, dbname = outputdb)

  dbWriteTable(con, "electrode", electrode, overwrite=TRUE, row.names=F)
  dbWriteTable(con, "neuron", neuron, overwrite=TRUE, row.names=F)
  dbWriteTable(con, "spikes", spikes, overwrite=TRUE, row.names=F)
  ## extra tables needed by SJE.
  dbWriteTable(con, "geometry", geometry, overwrite=TRUE, row.names=F)
  dbWriteTable(con, "distbreaks", distbreaks, overwrite=TRUE,row.names=F)

  dbDisconnect(con)

}
  

sql.read.spikes <- function(file,
                            time.interval=1,
                            beg=NULL, end=NULL,
                            min.rate=0,
                            corr.method="ci") {

  ## This could be a general reader... not just for sql.
  s1 = sql.spike.reader(file)

  spikes = s1$spikes
  layout = s1$layout

  
  channels = names(spikes)

  spikes.range <- range(unlist(spikes))
  if (is.null(beg))  beg <-  spikes.range[1]
  if (is.null(end))  end <-  spikes.range[2]
  rec.time <- c(beg, end)
  if (min.rate > 0 ) {
    
    ## Check for inactive channels.
    nspikes <- sapply(spikes,length)
    durn <- diff(rec.time)
    rates <- nspikes/durn
    inactive <- which(rates < min.rate)
    if (any(inactive)) {
      paste("Removing spikes with low firing rates: ",
            paste(inactive, collapse=' '))
      spikes   =   spikes[-inactive]
      channels = channels[-inactive]
    }
  }

  ## Count the number of spikes per channel, and label them.
  nspikes <- sapply(spikes, length)
  names(nspikes) <- channels

  ## meanfiring rate is the number of spikes divided by the (time of
  ## last spike - time of first spike).  

  ## TODO: SUN has this definition
  ##meanfiringrate <- nspikes/ ( sapply(spikes, max) - sapply(spikes, min))
  meanfiringrate <- nspikes/ ( end - beg)
  
  ## TODO; worry about multiple units recorded at the same location?
  
  unit.offsets <- NULL                  #default value.

  ## check that the spikes are monotonic.
  check.spikes.monotonic(spikes)

  rates <- make.spikes.to.frate(spikes, time.interval=time.interval,
                                beg=beg, end=end)

  res <- list( channels=channels,
              spikes=spikes, nspikes=nspikes, NCells=length(spikes),
              meanfiringrate=meanfiringrate,
              file=file,
              layout=layout,
              rates=rates,
              unit.offsets=unit.offsets,
              rec.time=rec.time
              )
  class(res) <- "mm.s"

  if (corr.method == "ci") {
    breaks = s1$distbreaks
    res$corr = corr.index(res, breaks)
  }

  res
  
}

sql.spike.reader <- function(file) {

  require(RSQLite)
  m <- dbDriver("SQLite")
  con <- dbConnect(m, dbname = file)

  rs <- dbSendQuery(con, "select * from electrode")
  electrode <- fetch(rs, n=-1)

  rs <- dbSendQuery(con, "select * from neuron")
  neuron <- fetch(rs, n=-1)
  

  e = neuron[,"electrode"]
  r = match(e, electrode[,1])
  res = electrode[r,]
  ##pos = data.frame(x=res[,"x"], y=res[,"y"], electrode.num=res[,"num"])
  pos = cbind(x=res[,"x"], y=res[,"y"], electrode.num=res[,"num"])
  rownames(pos) = neuron[,"name"]


  rs <- dbSendQuery(con, "select * from geometry")
  geometry <- fetch(rs, n=-1)


  rs <- dbSendQuery(con, "select * from distbreaks")
  distbreaks <- fetch(rs, n=-1)
  
  rs <- dbSendQuery(con, "select * from spikes")
  sp<- fetch(rs, n=-1)

  ## end of querying; close the connection.
  dbClearResult(rs)
  dbDisconnect(con)

  ## keep the order as they appear in the data frame.
  ##spikes = split(sp$t, factor(sp$neuron, levels=unique(sp$neuron)))

  ## extract spikes, and put names back onto the spikes.
  spikes = split(sp$t, sp[,"neuron"])
  names(spikes) <- neuron[,"name"]
  ## TODO, write as
  ## names(spikes) <- neuron[match(names(spikes),...],"name"]
  ## so that it looks up right name, rather than assuming it is in the
  ## right order.
  


  layout <- list(xlim=c(geometry[,"xlo"], geometry[,"xhi"]),
                 ylim=c(geometry[,"ylo"], geometry[,"yhi"]),
                 spacing=geometry[,"spacing"],
                 pos=pos)
  class(layout) <- "mealayout"

  
  res <- list(layout=layout, spikes=spikes, distbreaks=distbreaks[,1])
}


