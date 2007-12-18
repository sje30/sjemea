\name{jay.read.spikes}
\alias{jay.read.spikes}
\title{
Read in the .txt file from Neuroexplorer and create a "spikes" data structure.
}
\description{
Read in .txt file and work out array positions...
}
\usage{
jay.read.spikes(filename)
}
\arguments{
  \item{filename}{Name of the text file to be read in.}
}

\value{
  Return the data structure 's'.
}

\section{METHOD}{
  No fancy tricks used here.  If the data file has information about N
  different spike trains, the file has N (tab-separated) columns.  Each
  column then gives the time (in seconds?) of each spike.  Different
  columns are of different lengths since typically each cell will have a
  different number of spikes.

  The txt file of spike times can be compressed (with gzip).
}

%% No examples yet.
\examples{
data.file <- system.file("examples", "P9_CTRL_MY1_1A.txt",
                         package = "sjemea")
s <- jay.read.spikes( data.file)
fourplot(s)
s <- jay.read.spikes( data.file, beg=400, end=700)
fourplot(s)
\dontrun{
s <- jay.read.spikes("/home/stephen/ms/jay/p9data.txt")
fourplot(s)                             #summary plot.
s$mi <- make.mi(s)
show.prob.t.r(s)                        #conditional distributions.
}

\dontrun{crosscorrplots(s, autocorr=T, tmax=3, nbins=100,
               xcorr.nrows=3, xcorr.ncols=3) #plot autocorrs on screen

## Plotting just one cross-correlogram is a slightly different matter:
xcorr.plot( s$spikes[[1]], s$spikes[[2]], "1 v 2")}

}

\seealso{
  \code{\link{sanger.read.spikes}}
}

\references{
No references here.
}

%% Keywords from /usr/lib/R/doc/KEYWORDS.db
\keyword{math}
