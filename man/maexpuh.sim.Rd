\name{maexpuh.sim}
\alias{maexpuh.sim}
\title{
  Smoothed exponential stores
}
\description{
  Multiple exponential stores with pre- and post- smoothing on 2nd and 3rd stores
}
\usage{
maexpuh.sim(U, delay = 0, tau_s = 0, tau_q = 0, tau_3 = 0, v_s = 1, v_q = NA, v_3 = 0, series = 0, loss = 0, Xs_0 = 0, Xq_0 = 0, X3_0 = 0, w_s = 7, w_3 = 30, pars = NULL, return_components = FALSE, na.action = na.pass, epsilon = hydromad.getOption("sim.epsilon"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{U}{
%%     ~~Describe \code{U} here~~
}
  \item{delay}{
%%     ~~Describe \code{delay} here~~
}
  \item{tau_s}{
%%     ~~Describe \code{tau_s} here~~
}
  \item{tau_q}{
%%     ~~Describe \code{tau_q} here~~
}
  \item{tau_3}{
%%     ~~Describe \code{tau_3} here~~
}
  \item{v_s}{
%%     ~~Describe \code{v_s} here~~
}
  \item{v_q}{
%%     ~~Describe \code{v_q} here~~
}
  \item{v_3}{
%%     ~~Describe \code{v_3} here~~
}
  \item{series}{
%%     ~~Describe \code{series} here~~
}
  \item{loss}{
%%     ~~Describe \code{loss} here~~
}
  \item{Xs_0}{
%%     ~~Describe \code{Xs_0} here~~
}
  \item{Xq_0}{
%%     ~~Describe \code{Xq_0} here~~
}
  \item{X3_0}{
%%     ~~Describe \code{X3_0} here~~
}
  \item{w_s}{
%%     ~~Describe \code{w_s} here~~
}
  \item{w_3}{
%%     ~~Describe \code{w_3} here~~
}
  \item{pars}{
%%     ~~Describe \code{pars} here~~
}
  \item{return_components}{
%%     ~~Describe \code{return_components} here~~
}
  \item{na.action}{
%%     ~~Describe \code{na.action} here~~
}
  \item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  ## Pre and post filtered Xs
  Us[]=filter(U,rep(1/w_s,w_s),sides=1)
  Us[1:(w_s-1)] <- cumsum(U[1:(w_s-1)])/1:(w_s-1)
  Xstemp[] <- filter_loss(beta_s * U, alpha_s, loss = lossVal,
                      init = Xs_0)
  Xs[]=filter(Xstemp,rep(1/w_s,w_s),sides=1)
  Xs[1:(w_s-1)] <- cumsum(Xstemp[1:(w_s-1)])/1:(w_s-1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
