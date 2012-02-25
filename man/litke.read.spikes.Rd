\name{litke.read.spikes}
\alias{litke.read.spikes}
\alias{iit.read.spikes}
\title{
Read in the .mat file containing MEA data from Alan Litke
}
\description{
Read in the matlab file and process as a regular MEA.
}
\usage{
litke.read.spikes(filename, ids, time.interval, beg, end, corr.method)
iit.read.spikes(filename, ids, time.interval, beg, end)
}
\arguments{
  \item{filename}{Name of the text file to be read in.}
  \item{ids}{Optional vector of cell numbers that can be analysed,
    rather than analysing all electrodes in the recording.  Warning: Not
    implemented in all readers.}
  \item{time.interval}{Bin width (in seconds) for estimating firing
    rate.  Defaults to 1 second.}
  \item{beg}{Optional start time.}
  \item{end}{Optional end time.}
  \item{corr.method}{Which correlation method to use?  Only
    'ci' supported for the Likte reader.}
}

\value{
  Return the data structure 's'.
}

\section{METHOD}{
  The data provided by Alan Litke is stored as a matlab file.  There
  should be only one item in the file, which is the matrix containing
  the spikes.

  The matrix contains one row per electrode.  The first item in the row
  is the cell number; the second item is the electrode that recorded the
  cell.  The remaining numbers are the spike times (which must be
  divided by 20000 to get time in seconds).  All spike trains are
  zero-padded to be the same length.

  As the matlab format is used, the R.matlab and Rcompression packages are
  required.  These can be installed from R using:

  install.packages(c('R.utils', 'R.matlab'))
  install.packages('Rcompression',
  repos = "http://www.omegahat.org/R")

  (The IIT data also use a matlab format, and so the above installations
  are also needed for IIT.)

  The location of each electrode is given in the data file
  `litke1layout', see examples below.
  
}

\seealso{
  \code{\link{jay.read.spikes}}
}

\references{
This reader used for the Stafford et al. (2009) data in Neuron.
}

\examples{
## Show the layout of the array.
show.litke.layout()

## The location of the file storing the electrode positions
system.file('data', 'litke1layout.txt', package='sjemea')
}

%% Keywords from /usr/lib/R/doc/KEYWORDS.db
\keyword{math}
