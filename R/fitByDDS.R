fitByDDS<-function (MODEL, objective = hydromad.getOption("objective"), 
    control = hydromad.getOption("dds.control"), save=NULL) 
{
  library(ppso)
  start_time <- proc.time()
  objective <- buildCachedObjectiveFun(objective, MODEL)
  parlist <- as.list(coef(MODEL, warn = FALSE))
  isok <- sapply(parlist, function(x) !any(is.na(x)))
  parlist <- parlist[isok]
  isfixed <- (sapply(parlist, length) == 1)
  if (all(isfixed)) {
    warning("all parameters are fixed, so can not fit")
    return(MODEL)
  }
  parlist <- parlist[!isfixed]
  ##TODO: if (isTRUE(hydromad.getOption("trace"))) 
  lower <- sapply(parlist, min)
  upper <- sapply(parlist, max)
  initpars <- sapply(parlist, mean)
  bestModel <- MODEL
  bestFunVal <- -Inf
  do_dds <- function(pars) {
    thisMod <- update(MODEL, newpars = pars)
    if (!isValidModel(thisMod)) 
      return(NA)
    thisVal <- objFunVal(thisMod, objective = objective)
    if (isTRUE(thisVal > bestFunVal)) {
      bestModel <<- thisMod
      bestFunVal <<- thisVal
    }
    if(!is.null(save)) save(pars,thisVal,thisMod)
    return(-thisVal)
  }
  ans <- do.call("optim_dds",
                 modifyList(control,
                            list(
                                 objective_function = do_dds,
                                 number_of_parameters = length(initpars),
                                 parameter_bounds=cbind(lower,upper),
                                 initial_estimates=as.matrix(initpars)
                                 )
                            ))
  bestModel$msg <- ans$break_flag
  bestModel$funevals <- ans$function_calls
  bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
  bestModel$objective <- objective
  bestModel$fit.call <- match.call()
  bestModel$fit.result <- ans
  return(bestModel)
}
