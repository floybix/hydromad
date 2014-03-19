evalRSM <- function(modx,...,objective = hydromad.getOption("objective")){
  library(rsm)
  ## Sample using design
  bounds <- getFreeParsRanges(modx)
  
  coding=lapply(names(bounds),
    function(var) {
      max.var <- bounds[[var]][1]
      min.var <- bounds[[var]][2]
      as.formula(sprintf("%s.c ~ (2*%s-(%f+%f))/(%f-%f)",var,var,max.var,min.var,max.var,min.var))
    })
  
  evals <- ccd(as.formula(sprintf("~%s",paste(sprintf("%s.c",names(bounds)),collapse="+"))),
               coding=coding,inscribed=TRUE,...)

  cat(sprintf("Running %d model evaluations\n",nrow(evals)))
  evals$Y <- evalPars(decode.data(evals)[,names(bounds)],modx,objective=objective)

################################################################################
  ## Estimate quadratic response surface
  form <- as.formula(sprintf("Y~SO(%s)",paste(sprintf("%s.c",names(bounds)),collapse=",")))
  evals.rsm <- rsm(form, data=evals)

  return(evals.rsm)
}
