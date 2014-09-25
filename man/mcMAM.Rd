\name{mcMAM}
\alias{mcMAM}
\alias{mcMAM.default}
\title{
Estimating parameters of a minimum age model with a Markov chain Monte Carlo method
}
\description{
Sampling from the joint-likelihood function of a minimum age model using a Markov chain Monte Carlo method.
}
\usage{
mcMAM(EDdata, ncomp = -1, addsigma = 0, iflog = TRUE, 
      nsim = 50000, inis = list(), control.args = list())
}
\arguments{
  \item{EDdata}{\link{matrix}(\bold{required}): two columns, i.e., equivalent doses and associated errors}

  \item{ncomp}{\link{numeric}(with default): the number of components, \code{-1=MAM3}, \code{-2=MAM4}}

  \item{addsigma}{\link{numeric}(with default): an additional uncertainty}

  \item{iflog}{\link{logical}(with default): use log-scale or not}

  \item{nsim}{\link{numeric}(with default): the total number of iteration}

  \item{inis}{\link{list}(with default): initials of the chains. 
Example: \code{inis=list(p=0.1,gamma=20,sigma=0.5)} in MAM3}

  \item{control.args}{\link{list}(with default): arguments used by slice sampling method, see function \link{mcFMM} for details}
}

\value{
Return an invisible list of S3 class \code{"mcAgeModels"}, which can be operated with function \link{reportSAM}. See \code{value} of \link{mcFMM} for details.
}

\references{
Galbraith RF, Roberts RG, Laslett GM, Yoshida H, Olley JM, 1999. Optical dating of single grains of quartz from Jinmium rock shelter, northern Australia. Part I: experimental design and statistical models. Archaeometry, 41: 339-364.

Neal RM, 2003. "Slice sampling" (with discussion). Annals of Statistics, 31: 705-767. Software is available at \url{http://www.cs.utoronto.ca/~radford/slice.software.html}.
}

\seealso{
  \link{mcFMM}; \link{reportSAM}; \link{RadialPlotter}
}

\examples{
  data(EDdata)
  # Construct MCMC chains for a MAM3.
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
