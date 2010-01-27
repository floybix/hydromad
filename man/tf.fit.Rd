\name{tf.fit}
\alias{tf.fit}
\alias{tf.ls.fit}
\alias{tf.sriv.fit}
\title{ Estimate transfer function models from data. }
\description{
  Calibrate transfer function models (see \code{\link{tf}})
  using Least Squares (LS) or Simple Refined Instrumental Variable (SRIV)
  algorithms.
}
\usage{
tf.fit(DATA = list(U=, Q=),
	order = NULL, delay = NULL,
	estimator = c("ls", "sriv"),
	objective,
	warmup,
	na.action = na.pass,
	normalise = TRUE,
	trace, ...)

tf.ls.fit(DATA = list(U=, Q=),
	order = NULL, delay = NULL,
	weights = NULL,
	prefilter = TRUE, a.hat = NULL,
        objective,
	warmup,
	na.action = na.pass,
	trace)

tf.sriv.fit(DATA = list(U=, Q=),
	order = c(n=2, m=1), delay = 0,
	prefilter = TRUE, a.hat = a.init(1, n),
	X,
	warmup = 100, initX = TRUE,
	na.action = na.pass,
	normalise = TRUE,
	epsilon = 0.001, max.iterations = 10,
	trace = FALSE)

}
\arguments{
  \item{DATA}{ ~~Describe \code{U} here~~ }
  \item{order}{ ~~Describe \code{order} here~~ }
  \item{delay}{ ~~Describe \code{delay} here~~ }
  \item{estimator}{ ~~Describe \code{estimator} here~~ }
  \item{objective}{ ~~Describe \code{order} here~~ }
  \item{warmup}{ ~~Describe \code{warmup} here~~ }
  \item{na.action}{ ~~Describe \code{na.action} here~~ }
  \item{normalise}{ ~~Describe \code{normalise} here~~ }
  \item{trace}{ ~~Describe \code{trace} here~~ }
  \item{\dots}{ ~~Describe \code{trace} here~~ }
  \item{weights}{ ~~Describe \code{weights} here~~ }
  \item{prefilter}{ ~~Describe \code{prefilter} here~~ }
  \item{a.hat}{ ~~Describe \code{a.hat} here~~ }
  \item{X}{ This must have the same dimensions and time window as the input \code{U}. }
  \item{epsilon}{ ~~Describe \code{epsilon} here~~ }
  \item{max.iterations}{ ~~Describe \code{max.iterations} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{
Jakeman
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{ \code{\link{tf}}, \code{\link{arima}} }
\examples{

}
\keyword{ ts }

