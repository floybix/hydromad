\name{nseVarTd}
\alias{nseVarTd}
\alias{adjVarTd}
\title{
  Time-delay corrected performance measure
}
\description{
  adjVarTd coalesces modelled flow peaks to
  observed flow peaks for each event separately. nseVarTd calculates
  Nash-Sutcliffe efficiency on the result using \code{\link{nseStat}}.
  Depending on the quality of the coalescing, this better indicate
  performance ignoring timing error.
}
\usage{
nseVarTd(obs, mod, event, ...)
adjVarTd(obs,mod,event,...)
}
\arguments{
  \item{obs}{
    observed data vector
}
  \item{mod}{
    model-predicted data vector corresponding to \code{obs}.
}
  \item{event}{
    zoo object of events, as returned by \code{\link{eventseq}}
}
  \item{\dots}{
    Additional arguments to \code{\link{nseStat}}
}
}
\value{
  For nseVarTd, a single numeric value.
  For adjVarTd, a zoo object with the original modelled and observed
  data, the adjusted model output and the lag estimated for each event.
}
\author{
  Joseph Guillaume
}
\seealso{
  \code{\link{hydromad.stats}},\code{\link{nseStat}},\code{\link{objFunVal}}
}
\examples{

data(Murrindindi)
x <- Murrindindi[1:100]
x <- merge(x,X=lag(x$Q,2))

event <- eventseq(x$P, thresh = 5, inthresh = 3.5, indur = 7, continue = TRUE)

nseStat(x$Q,x$X)
nseVarTd(x$Q,x$X,event)
hmadstat("r.sq.vartd")(x$Q,x$X,event=event)
}
\keyword{ts}
