## Copied from ~/nobackup/stephen/langs/R/tripack/R/zzz.R
.First.lib <- function(lib, pkg) {
  library.dynam("sjemea", pkg, lib)
}

.Last.lib <- function (libpath) {
  ## Run when the package is being unloaded.  This allows us to easily test
  ## packages within same session when dynlib is updates.
  library.dynam.unload("sjemea", libpath)
}
  
