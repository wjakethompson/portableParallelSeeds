.onLoad <- function(libname, pkgname){
     RNGkind("L'Ecuyer-CMRG")
     .pps <- new.env(parent = emptyenv())
}
