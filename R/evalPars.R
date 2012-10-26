evalPars<- function(par.matrix,object,objective=hydromad.getOption("objective")) {
  stopifnot(inherits(object,"hydromad"))
  apply(par.matrix,1,function(p) {
    thisMod <- update(object, newpars = p)
    if (!isValidModel(thisMod)) return(NA)
    objFunVal(thisMod,objective)
  })
}


getFreeParsRanges <- function(object){
  stopifnot(inherits(object,"hydromad"))
  ## identify varying parameters
  par.ranges <- suppressWarnings(coef(object))
  free <- sapply(par.ranges, function(x) {
    !inherits(x, "AsIs") && length(x) == 2 && (diff(range(x)) > 
                                    0)
  })
  par.ranges[free]
}
