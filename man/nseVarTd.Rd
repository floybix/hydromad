\name{nseVarTd}
\alias{nseVarTd}
\title{
  Time-delay corrected performance measure
}
\description{
  Performance measure where modelled flow peaks have been coalesced to
  observed flow peaks for each event separately, i.e. gives a better
  indication of performance ignoring timing error
}
\usage{
nseVarTd(obs, mod, event, ...)
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
  a single numeric value
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
