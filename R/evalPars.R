evalPars<- function(par.matrix,object,objective=hydromad.getOption("objective"),
                    parallel=hydromad.getOption("parallel")[["evalPars"]],
                    async=FALSE,export=c()) {
  stopifnot(inherits(object,"hydromad"))
  switch(parallel,
         "foreach"={
           opts=hydromad.options()
           objs <- foreach(p=iter(par.matrix,by="row"),
                           .packages="hydromad",
                           .inorder=TRUE,
                           .export=export,
                           .final=function(x) as.numeric(x),
                           .options.redis=list(async=async)
           ) %dopar% {
             hydromad.options(opts)
             thisMod <- update(object, newpars = p)
             if (!isValidModel(thisMod)) return(NA)
             objFunVal(thisMod,objective=objective)
           }
         },
         "clusterApply"={
           objs <- parApply(cl=cl,par.matrix,1,
                                 function(p,object,objective) {
                                   thisMod <- update(object, newpars = p)
                                   if (!isValidModel(thisMod)) return(NA)
                                   objFunVal(thisMod,objective)
                                 },object=object,objective=objective)
         },
         objs <- apply(par.matrix,1,function(p) {
           thisMod <- update(object, newpars = p)
           if (!isValidModel(thisMod)) return(NA)
           objFunVal(thisMod,objective)
         })
         ) ##switch parallel
  return(objs)
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
