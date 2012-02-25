\name{litke1layout}
\alias{litke1layout}
\docType{data}
\title{
Layout of the MEA provided by Alan Litke, as used in the 2009 Neuron paper.
}
\description{
For each electrode name, the (x,y) location of the electrode (in units
of um) is given.
}
\usage{data(litke1layout)}
\format{
  A data frame with 512 observations on the following 3 variables.
  \describe{
    \item{\code{electrode}}{a numeric vector}
    \item{\code{x}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
  }
}
\source{
Stafford et al. (2009) Neuron.
}
\examples{
data(litke1layout)
## maybe str(litke1layout) ; plot(litke1layout) ...
}
\keyword{datasets}

%% THis file created using "promptData(litke1layout)"