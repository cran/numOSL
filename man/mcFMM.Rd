\name{mcFMM}
\alias{mcFMM}
\alias{mcFMM.default}
\title{
Estimating parameters of finite mixture age models (include the central age model) with a Markov chain Monte Carlo method
}
\description{
Constructing Markov chains for the parameters of interest in finite mixture age models (up to four components).
}
\usage{
mcFMM(EDdata, ncomp = 1, addsigma = 0, iflog = TRUE, 
      nsim = 50000, inis = list(), control.args = list())
}
\arguments{
  \item{EDdata}{\link{data.frame}(\bold{required}): equivalent doses and associated errors (a total of two columns), minus equivalent dose values can also be analyzed}
  \item{ncomp}{\link{numeric}(with default): the number of components used for optimizing, 1 for the central age model, >=2 for a finite mixture age model}
\item{addsigma}{\link{numeric}(with default): the added spread to the relative (or absolute) standard error of equivalent dose}
  \item{iflog}{\link{logical}(with default): use log-scale or not}
  \item{nsim}{\link{numeric}(with default): the total number of simulations (desired row number of the chains)}
  \item{inis}{\link{list}(with default): initials of the chains, such as inis=list(p1=1,p2=2,mu1=10,mu2=20) for a two-parameter Finite Mixture Age Model}
  \item{control.args}{\link{list}(with default): some further arguments (w, m, nstart) used in a Slice Sampling method, see details}
}
\details{
This function is used for constructing Markov chains for the desired parameters of finite mixture age models with the Slice Sampling method (Neal, 2003). The method is based on the fact 
that to sample a random variable one can sample uniformly from the region under the curve of its density function. With given initials and the prior information (a specified model, a number of
observations, feasible intervals of the parameters of interest), each parameter is regarded as an univariate variable and is updated in turns with a accept-reject protocol.
Three arguments (in control.args) are used to control the performance of the sampling:\cr\cr
\bold{w}: the size of the steps for creating an interval from which to sample, default w=1. \cr
\bold{m}: the limit on steps for expanding an interval, the interval is expanded to include more feasible values if possible. m<=1 means no limit on the expandation, m>1 means that the interval
is expanded with a finite number of iterations (depends on the magnitude of m), default m=-100.\cr
\bold{nstart}: the maximum number of trials for updating a variable in each simulation. It can be used to monitor the stability of the chains of a specified model. For example, the chains of a 4-parameter minimum age model is very likely to crash down and sometimes more number of attempts (nstart) are needed for a single simulation, default nstart=1.
}
\value{
Return an invisible list of class "mcAgeModels", which can be operated with function \link{reportSAM}.
}
\references{
Galbraith, R.F., Green, P., 1990. Estimating the component ages in a finite mixture. International Journal of Radiation Applications and Instrumentation. Part D. Nuclear Tracks and Radiation Measurements, 17 (3), pp. 197-206.

Neal, R.M., 2003. "Slice sampling" (with discussion). Annals of Statistics, 31 (3), pp. 705-767.
}

\seealso{
  \link{mcMAM}; \link{reportSAM}
}
\examples{
  data(EDdata)
  # Construct MCMC chains for a 3-component FMM.
  obj<-mcFMM(EDdata$gl11,ncomp=3,nsim=1e4)
  reportSAM(obj,thin=2,burn=1e3)
}
\keyword{Central Age Model}
\keyword{Finite Mixture Age Model}
\keyword{equivalent dose}
\keyword{Slice sampling}
