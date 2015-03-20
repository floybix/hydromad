## Calculate objective function on a rolling window
## Scalable in length of ts and number of parameters
## Parallelised with results stored on disk
## tempfile is automatically deleted when exiting from R
evalParsRollapply <- function(par.matrix,object,
                              width=30,
                              objective=hydromad.getOption("objective"),
                              parallel=hydromad.getOption("parallel")[["evalParsRollapply"]],
                              filehash.name=tempfile()
){
  
  # Sets default settings for parallelisation if missing
  parallel=hydromad.parallel(parallel)
  
  if(is.null(filehash.name)){
    results=matrix(NA,nrow=nrow(par.matrix),ncol=length(observed(object))-width+1)
    if(parallel$method!="none") 
      warning("setting parallel$method='none', parallelisation requires filehash.name to be non-null")
    parallel$method="none"
  } else {
    library(parallel)
    library(ff)
    results=ff(vmode="double",dim=c(nrow(par.matrix),length(observed(object))-width+1),filename=filehash.name)
  }
  
  ## Do runs, storing all ts
  cat(sprintf("Running %d model evaluations %s parallelisation\n", nrow(par.matrix),
              ifelse(parallel$method=="none","*without*","with")
  ))
  switch(parallel$method,
         "clusterApply"={
           lapply(c("ff",parallel$packages),function(pkg) clusterCall(cl,library,pkg,character.only=TRUE))
           if(length(parallel$export)>0) clusterExport(cl,parallel$export)
           clusterExport(cl,c("par.matrix","object","objective","width","results"),envir=environment())
           parLapply(cl=cl,as.list(1:nrow(par.matrix)),
                     function(ip) {
                       thisMod <- update(object, newpars = par.matrix[ip,])
                       if (!isValidModel(thisMod)) return(NULL)
                       res <- rollapply(cbind(Q=observed(thisMod),X=fitted(thisMod)),
                                        width=width,by.column=FALSE,
                                        FUN=objFunVal,objective=objective
                       )
                       results[ip,] <- res
                       ip
                     })
         },{
           lapply(as.list(1:nrow(par.matrix)),
                  function(ip) {
                    thisMod <- update(object, newpars = par.matrix[ip,])
                    if (!isValidModel(thisMod)) return(NULL)
                    res <- rollapply(cbind(Q=observed(thisMod),X=fitted(thisMod)),
                                     width=width,by.column=FALSE,
                                     FUN=objFunVal,objective=objective
                    )
                    results[ip,] <<- res
                    ip
                  })
         }) ##switch parallel
  ## ff matrix of ts
  results
}