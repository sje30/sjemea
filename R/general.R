## File for generally useful functions.


file.or.gz <- function(file) {
  ## Return FILE if it exists, or FILE.gz if that exists.
  ## Otherwise, return NA and generate a warning.
  if (file.exists(file)) {
    file
  } else {
    f2 <- paste(file,".gz", sep="")
    if (file.exists(f2))
      f2
    else {
      warning(paste("File", file,
                    "could not be found, nor its compressed version."))
      NA
    }
  }
}
