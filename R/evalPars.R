evalPars<- function(par.matrix,mod,objective=hydromad.getOption("objective")) {
  apply(par.matrix,1,function(p) {
    thisMod <- update(mod, newpars = p)
    if (!isValidModel(thisMod)) return(NA)
    objFunVal(thisMod,objective)
  })
}


getFreeParsRanges <- function(mod){
  ## identify varying parameters
  par.ranges <- suppressWarnings(coef(mod))
  free <- sapply(par.ranges, function(x) {
    !inherits(x, "AsIs") && length(x) == 2 && (diff(range(x)) > 
                                    0)
  })
  par.ranges.free <- par.ranges[free]
}
