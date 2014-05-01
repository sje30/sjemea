## Suggestion from Laurent for an enviroment
## to store package details.

.myPackageEnv <- new.env(parent = emptyenv(), hash = TRUE)
## default value
assign("myoption", TRUE, envir = .myPackageEnv)
## you might want to lock .myPackageEnv
getMyoption <- function() get("myoption", .myPackageEnv)
setMyoption <- function(x) assign("myoption", x, .myPackageEnv)

