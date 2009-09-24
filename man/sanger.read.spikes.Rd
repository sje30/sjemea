\name{sanger.read.spikes}
\alias{sanger.read.spikes}
\title{Read in a nexTimestamps file from Neuroexplorer and create
a "spikes" data structure.}
\description{
The Genes to Cognition Project (G2C) uses Neuroexplorer to output spike
times from an MEA recording.  This function reads those text files and
creates a "spikes" data structure.
}
\usage{
sanger.read.spikes(filename, ids=NULL, time.interval=1,
                   beg=NULL, end=NULL, min.rate=0)
}
\arguments{
  \item{filename}{Name of the text file to be read in.}

  \item{ids}{Which electrodes should be kept.}

  \item{time.interval}{How often to estimate the firing rate.}

  \item{beg}{Start time of the recording.}
  \item{end}{End time of the recording.}
  \item{min.rate}{Minimum firing rate of an electrode.}
}

\value{
  Return the data structure 's'.
}
\details{
  This format is fairly simple; each column represents one channel, with
  the spikes recorded in that channel filling up the column.  Since
  different channels are likely to record different numbers of spikes,
  empty values (just a TAB to next value) are found at the bottom of
  each spike train.

  Text files can be compressed, using gzip, to save space.  (The
  extension .gz is added to compresed files.)  This compression aspect
  should be transparent to the end user -- if a file is compresed, the
  filename can be given with or without the .gz extension.  gzip is
  available for Windows (from \url{http://www.gzip.org}); "gzip
  *.nexTimestamps" will compress all timestamps, and "gzip -d
  *.nexTimestamps" will uncompress them.  }

\seealso{
  \code{\link{jay.read.spikes}}
  \code{\link{compute.ns}}
}
\examples{
data.file <- system.file("examples", "TC89_DIV15_A.nexTimestamps",
                         package = "sjemea")
s <- sanger.read.spikes( data.file, beg=400, end=700)
fourplot(s)
}
\author{Stephen Eglen}
