.pps <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname){
     RNGkind("L'Ecuyer-CMRG")
}
