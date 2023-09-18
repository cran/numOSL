\name{mcMAM}
\alias{mcMAM}
\alias{mcMAM.default}
\title{
Optimization of the minimum (maximum) age model  (using a Markov chain Monte Carlo method)
}
\description{
Sampling from the joint-likelihood function of the minimum (maximum) age model using a Markov chain Monte Carlo (MCMC) method.
}
\usage{
mcMAM(EDdata, ncomp = -1, addsigma = 0, iflog = TRUE, 
      nsim = 50000, inis = list(), control.args = list())
}
\arguments{
  \item{EDdata}{\link{matrix}(\bold{required}): a two-column matrix (i.e., equivalent dose values and\cr
associated standard errors)}

  \item{ncomp}{\link{integer}(with default): number of components, \code{ncomp=-1}, \code{ncomp=-2}, \code{ncomp=-3}, or \code{ncomp=-4} indicate fitting the \code{"MAM3"}, \code{"MAM4"}, \code{"MXAM3"}, and \code{"MXAM4"}, respectively}

  \item{addsigma}{\link{numeric}(with default): additional uncertainty, i.e., the sigmab value}

  \item{iflog}{\link{logical}(with default): transform equivalent dose values to log-scale or not}

  \item{nsim}{\link{integer}(with default): deseried number of iterations}

  \item{inis}{\link{list}(with default): initial state of parameters.\cr
Example: \code{inis=list(p=0.1,gamma=20,sigma=0.3)} when \code{ncomp=-1}}

  \item{control.args}{\link{list}(with default): arguments used by the Slice Sampling algorithm, see function \link{mcFMM} for details}
}

\value{
Return an invisible \link{list} of S3 class object \code{"mcAgeModels"}. See \link{mcFMM} for details.
}

\references{
Galbraith RF, Roberts RG, Laslett GM, Yoshida H, Olley JM, 1999. Optical dating of single grains of quartz from Jinmium rock shelter, northern Australia. Part I: experimental design and statistical models. Archaeometry, 41(2): 339-364.

Neal RM, 2003. "Slice sampling" (with discussion). Annals of Statistics, 31(3): 705-767. Software is freely available at \url{https://glizen.com/radfordneal/slice.software.html}.

Peng J, Dong ZB, Han FQ, 2016. Application of slice sampling method for optimizing OSL age models used for equivalent dose determination. Progress in Geography, 35(1): 78-88. (In Chinese).
}

\seealso{
  \link{mcFMM}; \link{reportMC}; \link{RadialPlotter}; \link{EDdata}; \link{optimSAM}; \link{sensSAM}
}

\examples{
  # Not run.
  # data(EDdata)
  # Construct a MCMC chain for MAM3.
  # obj<-mcMAM(EDdata$al3,ncomp=-1,addsigma=0.1,nsim=5000)
  # reportMC(obj,burn=1e3,thin=2)
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
\keyword{Markov chain Monte Carlo}
