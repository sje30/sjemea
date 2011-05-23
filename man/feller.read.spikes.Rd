\name{feller.read.spikes}
\alias{feller.read.spikes}
\title{Read in a directory of spike times from Marla Feller and create
a "spikes" data structure.}
\description{
Marla Feller has provided MEA data from an hexagonal array.  The data
from each electrode is stored in a separate file.
This function reads those text files and
creates a "spikes" data structure.
}
\usage{
feller.read.spikes(filename, ids=NULL, time.interval=1,
                   beg=NULL, end=NULL)
}
\arguments{
  \item{filename}{Name of the text file to be read in.}

  \item{ids}{Which electrodes should be kept.}

  \item{time.interval}{How often to estimate the firing rate.}

  \item{beg}{Start time of the recording.}
  \item{end}{End time of the recording.}
}

\value{
  Return the data structure 's'.
}
\details{
  The data from each electrode is not in seconds; divide by 20000 to get
  the time in seconds.

}

\seealso{
  \code{\link{jay.read.spikes}}
  \code{\link{sanger.read.spikes}}
}
\examples{
data.file <- system.file('examples', '2002-11-26-1', package='sjemea')
s <- feller.read.spikes(data.file)
ejcmealayout.check()
fourplot(s)

## This file contains the picture of the array for comparison
system.file("examples/ejc_ElectrodeNumbering.jpg",package="sjemea")

## next two examples show electrodes being rejected/accepted.
s <- feller.read.spikes(data.file, ids=c('-', 'a7a', 'b1a', 'h2a'))
s <- feller.read.spikes(data.file, ids=c(     'a7a', 'b1a', 'h2a'))

}
\author{Stephen Eglen}
