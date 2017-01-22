\name{mcFMM}
\alias{mcFMM}
\alias{mcFMM.default}
\title{
Finite mixture age model optimization (using a Markov chain Monte Carlo method)
}
\description{
Sampling from the joint-likelihood functions of finite mixture age models (include the central age model) using a Markov chain Monte Carlo (MCMC) method.
}
\usage{
mcFMM(EDdata, ncomp = 1, addsigma = 0, iflog = TRUE, 
      nsim = 50000, inis = list(), control.args = list())
}
\arguments{
  \item{EDdata}{\link{matrix}(\bold{required}):  a two-column matrix (i.e., equivalent dose values and\cr
associated standard errors)}

  \item{ncomp}{\link{integer}(with default): number of components (\code{1} denotes the central age model)}

\item{addsigma}{\link{numeric}(with default): additional uncertainty}

  \item{iflog}{\link{logical}(with default): transform equivalent dose values to log-scale or not}

  \item{nsim}{\link{integer}(with default): deseried number of iterations}

  \item{inis}{\link{list}(with default): initial state of parameters.\cr 
Example: \code{inis=list(p1=1,p2=1,mu1=5,mu2=10)} in FMM2 (the sum of \code{p1} and \code{p2} will be normalized to 1 during the simulation)}

  \item{control.args}{\link{list}(with default): arguments used in the Slice Sampling algorithm, see details}
}

\details{
Function \link{mcFMM} is used for sampling from the joint-likelihood functions of finite mixture age models (include the central age model) using a Markov chain Monte Carlo sampling algorithm called Slice Sampling (Neal, 2003). Three arguments (\code{control.args}) are used for controling the sampling process:\cr
\bold{(1)} \emph{w}: size of the steps for creating an interval from which to sample, default \code{w=1}; \cr
\bold{(2)} \emph{m}: limit on steps for expanding an interval, \code{m<=1} means no limit on the expandation, \code{m>1} means the interval is expanded with a finite number of iterations, default \code{m=-100};\cr
\bold{(3)} \emph{nstart}: maximum number of trials for updating a variable in an iteration. It can be used for monitoring the stability of the simulation. For example, a MAM4 is likely to crash down for data sets with small numbers of data points or less dispersed distributions (see section 8.3 of Galbraith and Roberts, 2012 for a discussion), and sometimes more than one trial (i.e., using \code{nstart>1}) is required to complete the sampling process, default \code{nstart=1}.
}

\value{
Return an invisible \link{list} of S3 class object \code{"mcAgeModels"} including the following elements:
  \item{EDdata}{equivalent dose values}
  \item{addsigma}{additional uncertainty}
  \item{model}{fitting model}
  \item{iflog}{transform equivalent dose values to log-scale or not}
  \item{nsim}{number of iterations}
  \item{chains}{simulated samples}
}

\references{
Galbraith RF, Green P, 1990. Estimating the component ages in a finite mixture. International Journal of Radiation Applications and Instrumentation. Part D. Nuclear Tracks and Radiation Measurements, 17: 197-206.

Neal RM, 2003. "Slice sampling" (with discussion). Annals of Statistics, 31(3): 705-767. Software is freely available at \url{http://www.cs.utoronto.ca/~radford/slice.software.html}.
}

\seealso{
  \link{mcMAM}; \link{reportSAM};\cr
  \link{RadialPlotter}; \link{EDdata}
}

\examples{
  # Not run.
  # data(EDdata)
  # Construct a MCMC chain for FMM3.
  # obj<-mcFMM(EDdata$gl11,ncomp=3,nsim=5000)
  # reportSAM(obj,thin=2,burn=1e3)
}
\keyword{statistical age models}
\keyword{Central Age Model}
\keyword{Finite Mixture Age Model}
\keyword{equivalent dose}
\keyword{Markov chain Monte Carlo}
