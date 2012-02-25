\name{networkspikes}
\alias{compute.ns}
\alias{show.ns}
\title{Compute network spikes}
\description{
  Compute the network spikes in an MEA recording, by averaging over all
  the electrodes in the array.
}
\usage{
compute.ns(s, ns.T, ns.N,sur,plot)
}
\arguments{

  \item{s}{MEA data structure}

  \item{ns.T}{Bin width (in msec) for counting spikes.}

  \item{ns.N}{Threshold number of active electrodes required to make
    network spike}

  \item{sur}{How many bins either side of peak to retain when computing
    the mean network spike}

  \item{plot}{Set to TRUE to plot network spikes.}
  \item{ns}{A network spike data structure, returned by
    \code{\link{compute.ns}}}

  \item{...}{Other plot arguments to pass to \code{\link{show.ns}}}
}
\value{A list with the following elements:
  \item{counts}{vector giving the number of active electrodes in each
    bin; this can be very long!}
  \item{ns.N}{The value of ns.N used.}
  \item{ns.T}{the value of ns.T used.}
  \item{mean}{The profile of the mean network spike (this is a time
    series object)}
  \item{measures}{If N network spikes were found, this is a matrix with N
    rows, one per network spike.}
  \item{brief}{A short vector summarizing the network spikes.}
}
\details{
  To see the mean network spikes after they have computed, just look at
  the mean object.
  
  If you wish to see the individual network spikes, try show.ns(ns, ...)
  where the remaining args are passed to the plot function.
}
\references{
  Eytan and Marom (2006) J Neuroscience.
}
\seealso{\code{\link{sanger.read.spikes}}}
\examples{
data.file <- system.file("examples", "TC89_DIV15_A.nexTimestamps",
                         package = "sjemea")
s <- sanger.read.spikes( data.file, beg=400, end=700)
s$ns <- compute.ns(s, ns.T=0.003, ns.N=10,sur=100)
plot(s$ns)
plot(s$ns, xlim=c(450, 500))

plot(s$ns$mean)                         #show mean NS.
summary(s$ns)
s$ns$brief

}
\author{Stephen Eglen}
\keyword{Network spikes, MEA analysis}
