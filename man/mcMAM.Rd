\name{mcMAM}
\alias{mcMAM}
\alias{mcMAM.default}
\title{
Minimum age model optimization (with a Markov chain Monte Carlo sampling method)
}
\description{
Sampling from the joint-likelihood function of the minimum age model using a Markov chain Monte Carlo sampling method (MCMC).
}
\usage{
mcMAM(EDdata, ncomp = -1, addsigma = 0, iflog = TRUE, 
      nsim = 50000, inis = list(), control.args = list())
}
\arguments{
  \item{EDdata}{\link{matrix}(\bold{required}): a two-column matrix, equivalent dose values and associated errors are stored in the first and second column, respectively}

  \item{ncomp}{\link{integer}(with default): number of components, \code{-1=MAM3}, \code{-2=MAM4}}

  \item{addsigma}{\link{double}(with default): additional uncertainty}

  \item{iflog}{\link{logical}(with default): transform equivalent dose values to log-scale or not}

  \item{nsim}{\link{integer}(with default): deseried number of iterations}

  \item{inis}{\link{list}(with default): initial state of parameters.\cr
Example: \code{inis=list(p=0.1,gamma=20,sigma=0.5)} in MAM3}

  \item{control.args}{\link{list}(with default): arguments used by the Slice Sampling algorithm, see function \link{mcFMM} for details}
}

\value{
Return an invisible list of S3 class \code{"mcAgeModels"}. See \link{mcFMM} for details.
}

\references{
Galbraith RF, Roberts RG, Laslett GM, Yoshida H, Olley JM, 1999. Optical dating of single grains of quartz from Jinmium rock shelter, northern Australia. Part I: experimental design and statistical models. Archaeometry, 41(2): 339-364.

Neal RM, 2003. "Slice sampling" (with discussion). Annals of Statistics, 31(3): 705-767. Software is freely available at \url{http://www.cs.utoronto.ca/~radford/slice.software.html}.
}

\seealso{
  \link{mcFMM}; \link{reportSAM}; \link{RadialPlotter}
}

\examples{
  data(EDdata)
  # Construct a MCMC chain for MAM3.
  obj<-mcMAM(EDdata$al3,ncomp=-1,addsigma=0.1,nsim=5000)
  reportSAM(obj,burn=1e3,thin=2)
  #
  # The convergence of the simulations may be diagnosed with 
  # the Gelman and Rubin's convergence diagnostic.
  # library(coda)   # Only if package "coda" has been installed.
  # args<-list(nstart=50)
  # inis1<-list(p=0.01,gamma=26,mu=104,sigma=0.01)
  # inis2<-list(p=0.99,gamma=100,mu=104,sigma=4.99)
  # obj1<-mcMAM(EDdata$al3,ncomp=-2,nsim=3000,inis=inis1,control.args=args)
  # obj2<-mcMAM(EDdata$al3,ncomp=-2,nsim=3000,inis=inis2,control.args=args)
  # chain1<-mcmc(obj1$chains)
  # chain2<-mcmc(obj2$chains)
  # chains<-mcmc.list(chain1,chain2)
  # gelman.plot(chains)
}
\keyword{statistical age models}
\keyword{Minimum Age Model}
\keyword{equivalent dose}
