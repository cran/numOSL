\name{mcFMM}
\alias{mcFMM}
\alias{mcFMM.default}
\title{
Estimating parameters of a finite mixture age model (central age model) with a Markov chain Monte Carlo method
}
\description{
Sampling from the joint-likelihood function of a finite mixture age model (central age model) using a Markov chain Monte Carlo method.
}
\usage{
mcFMM(EDdata, ncomp = 1, addsigma = 0, iflog = TRUE, 
      nsim = 50000, inis = list(), control.args = list())
}
\arguments{
  \item{EDdata}{\link{matrix}(\bold{required}): two columns, e.g., a number of equivalent doses (un-logged) and associated errors}

  \item{ncomp}{\link{numeric}(with default): the number of components (\code{1} means a central age model)}

\item{addsigma}{\link{numeric}(with default): additional uncertainty}

  \item{iflog}{\link{logical}(with default): use log-scale or not}

  \item{nsim}{\link{numeric}(with default): the deseried total number of iterations}

  \item{inis}{\link{list}(with default): initials of the chain. Example, \code{inis=list(p1=1,p2=1,mu1=5,mu2=10)} in FMM2, here the sum of \code{p1} and \code{p2} subject to normalisation}

  \item{control.args}{\link{list}(with default): arguments used in a slice sampling method, see details}
}

\details{
This function is used to sampling from the joint-likelihood function of a finite mixture age model (central age model) using a Markov chain Monte Carlo method called slice sampling (Neal, 2003). Three arguments (\code{control.args}) are used to control the performance of the sampling process:\cr\cr
\emph{w}: the size of the steps for creating an interval from which to sample, default \code{w=1} \cr\cr
\emph{m}: the limit on steps for expanding an interval, the interval is expanded to include more feasible values if possible. \code{m<=1} means no limit on the expandation, \code{m>1} means that the interval is expanded with a finite number of iterations. Default \code{m=-100} \cr\cr
\emph{nstart}: the allowed maximum number of trials for updating a variable in each iteration. It can be used to monitor the stability of the chains of a model. For example, the chain of MAM4 is very likely to crash down and sometimes more number of attempts (\code{nstart>1}) is needed to complete the simulation. Default \code{nstart=1}
}

\value{
Return an invisible list of S3 class \code{"mcAgeModels"} including following elements, which can be operated with function \link{reportSAM}:
  \item{EDdata}{equivalent dose values}
  \item{addsigma}{additional uncertainty}
  \item{model}{the fitting model}
  \item{iflog}{use log-scale or not}
  \item{nsim}{the number of iterations}
  \item{chains}{the simulated quantities of interest}
}

\references{
Galbraith RF, Green P, 1990. Estimating the component ages in a finite mixture. International Journal of Radiation Applications and Instrumentation. Part D. Nuclear Tracks and Radiation Measurements, 17: 197-206.

Neal RM, 2003. "Slice sampling" (with discussion). Annals of Statistics, 31: 705-767. Software is available at \url{http://www.cs.utoronto.ca/~radford/slice.software.html}.
}

\seealso{
  \link{mcMAM}; \link{reportSAM}; \link{RadialPlotter}
}

\examples{
  data(EDdata)
  # Construct a MCMC chain for FMM3.
  obj<-mcFMM(EDdata$gl11,ncomp=3,nsim=5000)
  reportSAM(obj,thin=2,burn=1e3)
}
\keyword{statistical age models}
\keyword{Central Age Model}
\keyword{Finite Mixture Age Model}
\keyword{equivalent dose}
