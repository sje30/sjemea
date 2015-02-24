## Code specific to the work with Cambridge PDN 
## Copied from sanger.init()
## 2015-02-23

pdn.init <- function() {
  ## Run for initialisation of analysis of PDN data.

  windows <- .Platform$OS.type !='unix'

  guess.paths <- c("c:/mea2015", "~/Work/PhD/mea2015", "~/proj/mea2015")
  ## The first of the above folders that is readable will be used.

  ## the structure should include
  ## data/timestamps -- raw data
  ## data/hdf5 -- files that can be regenerated from data/timestamps
  valid.paths <- sapply(guess.paths, file.exists)
  first.valid <- which(valid.paths)[1]
  if (is.na(first.valid)) {
    stop("None of the following paths were valid. ", paste(guess.paths, collapse=' '))
  } else {
    path <- guess.paths[first.valid]
    assign("mea.data.dir",  paste0(path,"/data/"),   envir = .GlobalEnv)
    assign("mea.table.dir", paste0(path,"/tables/"), envir = .GlobalEnv)
    assign("mea.op.dir",    paste0(path,"/op/"),     envir = .GlobalEnv)
    assign("mea.R.dir",     paste0(path,"/R/"),      envir = .GlobalEnv)
  }

  if (windows) {
    ## Set up on Windows machine.
    ## Scripts are stored in meadev/scripts.
    setwd(mea.R.dir)
  }

  ## Create the cache of datafiles.
  ## perhaps let the user do this.
  ##assign("mea.data.files",  make.meafile.cache(mea.data.dir),
  ##envir  = .GlobalEnv)
}
