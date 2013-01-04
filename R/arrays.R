get.array.info <- function(data) {

  ## Array-specific information; maybe this could go in a file, rather
  ## than be read-in separately.
  
  pos <- data$epos;  rownames(pos) <- data$names
  array <- data$array
  if (array == 'stanford_hex_60um') {
    xlim <- ylim <- c(-300, 300)
    spacing <- 60
    corr.breaks <-  c(0, seq(35, by=70, length=9))
  }



  layout <- list(xlim = xlim, ylim = ylim, spacing = spacing, 
                 pos = pos)
  class(layout) <- "mealayout"
  list(layout=layout, corr.breaks=corr.breaks)
}
