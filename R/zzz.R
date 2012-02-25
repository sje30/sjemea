## zzz.R --- general functions for loading/unloading package.
## Author: Stephen J Eglen
## Copyright: GPL

.Last.lib <- function (libpath) {
  ## Run when the package is being unloaded.  This allows us to easily test
  ## packages within same session when dynlib is updates.
  library.dynam.unload("sjemea", libpath)
}
  
