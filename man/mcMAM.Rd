\name{mcMAM}
\alias{mcMAM}
\alias{mcMAM.default}
\title{
Estimating parameters of minimum age models (3- or 4-parameter) with a Markov chain Monte Carlo method
}
\description{
Constructing Markov chains for the parameters of interest in minimum age models.
}
\usage{
mcMAM(EDdata, ncomp = -1, addsigma = 0, iflog = TRUE, 
      nsim = 50000, inis = list(), control.args = list())
}
\arguments{
  \item{EDdata}{\link{data.frame}(\bold{required}): equivalent doses and associated errors (a total of two columns), minus equivalent dose values can also be analyzed}
  \item{ncomp}{\link{numeric}(with default): the number of components used for optimizing, -1 for the MAM3, -2 for the MAM4}
\item{addsigma}{\link{numeric}(with default): the added spread to the relative (or absolute) standard error of equivalent dose}
  \item{iflog}{\link{logical}(with default): use log-scale or not}
  \item{nsim}{\link{numeric}(with default): the total number of simulations (desired row number of the chains)}
  \item{inis}{\link{list}(with default): initials of the chains, such as inis=list(p=0.2,gamma=20,sigma=0.2) for a 3-parameter minimum age model}
  \item{control.args}{\link{list}(with default): some further arguments (w, m, nstart) used in a Slice Sampling method, see \link{mcFMM} for details}
}
\details{
see \link{mcFMM} for details.
}
\value{
Return an invisible list of class "mcAgeModels", which can be operated with function \link{reportSAM}.
}
\references{
Galbraith, R.F., Roberts, R.G., Laslett, G.M., Yoshida, H. & Olley, J.M., 1999. Optical dating of single grains of quartz from Jinmium rock shelter, northern Australia. Part I: experimental design and statistical models. Archaeometry, 41 (2), pp. 339-364.

Neal, R.M., 2003. "Slice sampling" (with discussion). Annals of Statistics, 31 (3), pp. 705-767.
}
\seealso{
  \link{mcFMM}; \link{reportSAM}
}
\examples{
  data(EDdata)
  # Construct MCMC chains for a MAM3.
  obj<-mcMAM(EDdata$al3,ncomp=-1,addsigma=0.1,nsim=2e4)
  reportSAM(obj,thin=10,burn=5e3)
}
\keyword{Minimum Age Model}
\keyword{equivalent dose}
\keyword{Slice sampling}
