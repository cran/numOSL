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
Return an invisible list of class "mcAgeModels" (has the same members as that returned by function \link{mcFMM}), which can be operated with function \link{reportSAM}.
}
\references{
Galbraith, R.F., Roberts, R.G., Laslett, G.M., Yoshida, H. & Olley, J.M., 1999. Optical dating of single grains of quartz from Jinmium rock shelter, northern Australia. Part I: experimental design and statistical models. Archaeometry, 41 (2), pp. 339-364.

Neal, R.M., 2003. "Slice sampling" (with discussion). Annals of Statistics, 31 (3), pp. 705-767.
}
\seealso{
  \link{mcFMM}; \link{reportSAM}; \link{RadialPlotter}
}
\examples{
  data(EDdata)
  # Construct MCMC chains for a MAM3.
  obj<-mcMAM(EDdata$al3,ncomp=-1,addsigma=0.1,nsim=1e4)
  reportSAM(obj,burn=1e3)
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
\keyword{Minimum Age Model}
\keyword{equivalent dose}
