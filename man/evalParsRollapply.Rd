\name{evalParsRollapply}
\alias{evalParsRollapply}
\title{
Calculate an objective function on a rolling time series for a matrix of parameters
}
\description{
For each row of a named matrix of parameters, run a model and calculate an objective function on rolling windows of the resulting time series, i.e. see how the objective function changes over time. Special facilities are provided to evaluate large sample sizes in parallel on a single multicore machine.
}
\usage{
evalParsRollapply(par.matrix, object, width = 30, 
  objective = hydromad.getOption("objective"), 
  parallel = hydromad.getOption("parallel")[["evalParsRollapply"]], 
  filehash.name = tempfile())
}
\arguments{
  \item{par.matrix}{
  Named matrix or data.frame of parameter values, with each row corresponding to a model realisation to evaluate
}
  \item{object}{
    an object of class \code{hydromad}.
}
  \item{width}{
integer specifying window width, aligned center. Passed to \code{\link[zoo]{rollapply}}
}
  \item{objective}{
    the objective function or expression, which can refer to Q and X.  See \code{\link{objFunVal.hydromad}}
}
  \item{parallel}{
If \code{"clusterApply"}, evaluate parameters in parallel using a local cluster. The implementation assumes that the \code{\link{ff}} file \code{filehash.name} can be written to simultaneously all cluster workers.
}
  \item{filehash.name}{
Name of \code{\link{ff}} file in which to store results, allowing large samples that do not fit in memory. Defaults to \code{tempfile()}, which is automatically deleted when exiting from R. To store results in memory, set \code{filehash.name=NULL}.
}
}
\value{
Either a matrix or \code{\link{ff}} file-backed matrix, with each row being a time series of rolling objective functions for the corresponding row of \code{par.matrix}
}
\note{
 When using \code{ff}, performance may be improved by specifying \code{options(ffcaching='mmeachflush')}.
}
\references{
Herman, J. D., P. M. Reed, and T. Wagener. 2013. "Time-Varying Sensitivity Analysis Clarifies the Effects of Watershed Model Formulation on Model Behavior." Water Resources Research 49 (3): 1400-1414. doi:\href{http://dx.doi.org/10.1002/wrcr.20124}{10.1002/wrcr.20124}
}
\author{
Joseph Guillaume
}
\seealso{
\code{\link{evalPars}} to calculate a single objective function for each parameter set
}
\examples{
data(Cotter)
obs <- Cotter[1:1000]

## Define rainfall-runoff model structure
object <- hydromad(obs, sma = "cwi", routing = "expuh",
                   tau_q=c(0,2),tau_s = c(2,100), v_s = c(0,1))

## Set the random seed to obtain replicable results
set.seed(19)

# Draw 10 Latin Hypercube random samples
par.matrix<-parameterSets(getFreeParsRanges(object),samples=10)
# Calculate rolling time series of r.squared for each parameter set,
#   keeping results in memory
runs <- evalParsRollapply(par.matrix,object,
  objective=hmadstat("r.squared"),filehash.name=NULL)


\dontrun{
## Setup parallelisation on three cores
library(parallel)
hydromad.options("parallel"=list("evalParsRollapply"="clusterApply"))
cl <- makeCluster(3)
clusterEvalQ(cl,library(hydromad))

par.matrix<-parameterSets(getFreeParsRanges(object),samples=1000)
# Calculate rolling time series of r.squared for each parameter set,
#  storing result in tempfile()
# Takes about 2 minutes
runs <- evalParsRollapply(par.matrix,object,
  objective=hmadstat("r.squared"))

# Excerpt of results
runs
# Path of backing file - about 7MB
filename(runs)
# ff object can be used like a regular matrix, e.g. plotting the 
#  rolling time series of R.squared for the first parameter set
plot(runs[1,])
}
}
\keyword{ models }
