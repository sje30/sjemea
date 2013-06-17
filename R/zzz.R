## zzz.R --- general functions for loading/unloading package.
## Author: Stephen J Eglen
## Copyright: GPL


.onUnload <- function (libpath) {
  ## Run when the package is being unloaded.  This allows us to easily test
  ## packages within same session when dynlib is updated.
  ## Detach package using: unloadNamespace("sjemea")
  library.dynam.unload("sjemea", libpath)
}
  
