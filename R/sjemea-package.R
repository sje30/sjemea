##' Multielectrode array (MEA) analysis in R
#' 
#' Collection of tools for reading in and processing MEA data.
#' ...
#' @name sjemea-package
#' @aliases sjemea
#' @docType package
#' @title Read and process multielectrode array data.
#' @author \email{sje30@@cam.ac.uk}
#' @keywords package
#' @seealso \code{\link{http://damtp.cam.ac.uk/user/eglen/waverepo}}



##' Layout of the MEA provided by Alan Litke, as used in the 2009 Neuron paper.
##' 
##' For each electrode name, the (x,y) location of the electrode (in units of
##' um) is given.
##' 
##' 
##' @name litke1layout
##' @docType data
##' @format A data frame with 512 observations on the following 3 variables.
##' \describe{ \item{list("electrode")}{a numeric vector} \item{list("x")}{a
##' numeric vector} \item{list("y")}{a numeric vector} }
##' @source Stafford et al. (2009) Neuron.
##' @keywords datasets
##' @examples
##' 
##' data(litke1layout)
##' ## maybe str(litke1layout) ; plot(litke1layout) ...
##' 
NULL


