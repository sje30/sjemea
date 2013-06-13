get.array.info <- function(data) {
  ## Array-specific information; maybe this could go in a file, rather
  ## than be read-in separately.  Useful for the HDF5 functions.
  
  pos <- data$epos;  rownames(pos) <- data$names
  array <- data$array

  if (array == 'stanford_hex_60um') {
    ## e.g. Wong1993
    xlim <- ylim <- c(-300, 300)
    spacing <- 60
    corr.breaks <-  c(0, seq(35, by=70, length=9))
  }

  if (array == 'litke_hex_60um') {
    xlim <- c(-1000, 1000)
    ylim <- c(-500, 500)
    spacing <- 60
    corr.breaks <- c(0, 0.001, seq(from=50, to=2050, by=50))
  }

  if ( any(grep('^Axion', array))) {
    ## e.g. Neurotox ongoing project.
    xlim <- c(0, 8000)
    ylim <- c(0, 6000)
    spacing <- 200
    corr.breaks <-  0                   #TODO; by default, do no breaks!
  }



  ## HACK for HDF5: the arrayname (stored in array) is a string, which
  ## can either be stored as a character vector or as an array.  R was
  ## getting confused about this, and all.equal() was failing...
  ## Need a clearer example...
  
  array <- as.character(array)
  layout <- list(xlim = xlim, ylim = ylim, spacing = spacing, 
                 pos = pos, array=array)
  class(layout) <- "mealayout"
  list(layout=layout, corr.breaks=corr.breaks)
}
