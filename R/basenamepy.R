## basenamepy.R --- Implementation of python's basename function.
## Author: Stephen J Eglen
## Copyright: GPL

basenamepy <- function (f) {
  ## Separate file name F into dir, base, and extn.
  dir <- dirname(f)
  file <- basename(f)

  ## Within file, find the last period (some files, such as a.b.R,
  ## have multiple periods.)
  
  periods <- gregexpr("\\.", file)
  last.periods <- sapply(periods, function(x) { max(x)})

  ## Negative values of last.periods indicates that there is no period
  ## in the filename.
  if ( any (neg <- which(last.periods < 0)) )
    last.periods[neg] <- sapply(file, nchar)[neg] + 1
  
  base <- substring(file, first=1, last=last.periods-1)
  extn <- substring(file, first=last.periods)

  list(dir=dir, base=base, extn=extn)
}

## some test functions.
if (FALSE) {
  basenamepy( c("/some/long/dir/cat.py", "/other/dir/apple.b.py"))
  
  filenames <-  c("/some/long/dir/cat.py", "/some/f/d/noextn", "simple.R",
                  "notmuch", "/home/dir/", "/other/dir/apple.b.py")
  r <- basenamepy(filenames)
  cbind(filenames, r$dir, r$base, r$extn)
}
