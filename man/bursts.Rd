\name{bursts}
\alias{spikes.to.bursts}
\alias{spikes.to.bursts.surprise}
\title{
  Burst detection of MEA spike trains.
}
\description{For a set of spike trains in an MEA recording, find the
  bursts independently within each spike train.}
\usage{
spikes.to.bursts(s, method="si")
}

\arguments{
  \item{s}{MEA data structure}
  \item{method}{A string, either "si" (surprise method), "mi"
  (maxinterval), "logisi" (Log ISI histogram).}
}

\value{
  Return the "all bursts" data structure.  This is a list of matrices,
  giving the burst information for each electrode.

  Each matrix stores basic information about each burst.  There is one
  row for every burst, with the following columns:

  \tabular{ll}{
    beg \tab index of the first spike in the burst \cr
    len \tab  number of spikes in this burst \cr
    SI \tab  surprise index (calculated only for the surprise method)\cr
    durn \tab  duration (in s) of the burst\cr
    mean.isis \tab  mean of all interspike intervals.\cr
  }

  If no bursts could be found within a spike train, the value NA is used
  rather than an empty matrix.
}

\examples{
data.file <- system.file("examples", "TC89_DIV15_A.nexTimestamps",
package = "sjemea")
s <- sanger.read.spikes(data.file)
s$allb <- spikes.to.bursts.surprise(s)
if (interactive()) spikeview(s)
}
\keyword{Burst analysis, MEA analysis}