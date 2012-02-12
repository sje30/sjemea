\name{hist.ab}
\alias{hist.ab}
\alias{histbi.ab}
\alias{count.nab}
\alias{test.count.hist2.nab}
\alias{test.count.hist.nab}
\alias{test.hist.ab}
\alias{test.histograms.versus.r}
\title{
  Histogram routines to help compute the cross-correlation between a
  pair of spike trains.
}
\description{
  For a pair of spike trains, TA (train a) and TB, these related
  routines return the count of the number of spikes in B that occur
  within some time window [-tmax,tmax] of a spike in A.  For histbi.ab,
  we return a histogram of dt values from [-tmax,tmax].  For hist.ab, we
  ignore the sign of each dt and just return a histogram in the range
  [0,tmax].  Finally, for count.nab, we just return the number of dt
  values found in the range [-tmax,tmax], rather than binning dt into a
  histogram.}

\usage{
hist.ab(ta, tb, tmax, nbins)
histbi.ab(ta, tb, tmax, nbins)
count.nab(ta, tb, tmax)
}
\arguments{
  \item{ta}{Vector of spike times, sorted such that lowest is first.}
  \item{tb}{Vector of spike times, sorted such that lowest is first.}
  \item{tmax}{maximum time (in seconds) to bin}
  \item{nbins}{Number of bins in the histogram.  For histbi.ab, each bin
    is of width (2*tmax)/nbins.  For hist.ab, each bin is (tmax)/nbins wide.}
}

\value{
  \code{hist.ab} returns a histogram of counts ignoring sign.
  \code{histbi.ab} returns a histogram of counts respecting sign.
  \code{count.nab} returns the number of dt values.
}

\section{METHOD}{
  For the histogram routines, each bin is of the form [low, high), with
  the exception of the last bin (for +tmax), which is of the form
  [tmax-binwid, tmax].  By assuming the spikes are ordered lowest first,
  the number of spike comparisons is greatly reduced, rather than
  comparing each spike with A with each spike in B.}

%% No examples yet.
\examples{

stopifnot(isTRUE(all.equal.numeric(
  histbi.ab(c(0), c(-2, -2, 0, 0, 1, 1,1, 1.8,2), tmax=2, nbins=4),
  c(2,0,2,5),
  check.attributes=FALSE)))
stopifnot(identical(TRUE, all.equal.numeric(
  hist.ab(c(0), c(-2, -2, 0, 0, 1, 1,1, 1.8,  2), tmax=2, nbins=4),
  c(2,0,3,4),
  check.attributes=FALSE)))
 
test.hist.ab()


%% Following examples are not run since they either take a long time
%% or require an "s" structure.
\dontrun{
test.histograms.versus.r()
test.count.hist.nab()
test.count.hist.nab(s)
test.count.hist2.nab(s)
}
}

\seealso{Nothing else yet...}

\references{
No references here.
}

%% Keywords from /usr/lib/R/doc/KEYWORDS.db
\keyword{math}

