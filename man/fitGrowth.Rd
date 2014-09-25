\name{fitGrowth}
\alias{fitGrowth}
\alias{fitGrowth.default}
\title{
Fitting a growth curve
}
\description{
Fitting a growth curve using the Levenberg-Marquardt algorithm.
}
\usage{
fitGrowth(Curvedata, model = c("line","exp","lexp","dexp"), 
          origin = FALSE, nstart = 100, upb = 0.5, plot = TRUE)
}

\arguments{
  \item{Curvedata}{\link{matrix}(\bold{required}): three columns, i.e., regenerative doses, standardised OSL signals, and standard errors of 
standardised OSL signals
}
  \item{model}{\link{character}(with default): a fitting model, default \code{model="line"}, see details
}
  \item{origin}{\link{logical}(optional): force the fitting to pass the origin or not
}
  \item{nstart}{\link{numeric}(with default): the maximum number of attempts of initializing the non-linear growth curve fitting, see details
}
  \item{upb}{\link{numeric}(with default): the upper limit on b value (the reciprocal of the saturating dose), initial b value is generated uniformly from the space \code{(0,upb)}, see details
}
  \item{plot}{\link{logical}(with default): draw a plot or not
}
}

\details{
For fitting a growth curve, four models are available: \cr
\bold{(1)} \code{"line"}: a linear model, \code{y=a*x+b};\cr
\bold{(2)} \code{"exp"}: an exponential model, \code{y=a*[1-exp(-b*x)]+c};\cr
\bold{(3)} \code{"lexp"}: an exponential plus linear model, \code{y=a*[1-exp(-b*x)]+c*x+d};\cr
\bold{(4)} \code{"dexp"}: an double exponential model, \code{y=a*[1-exp(-b*x)]+c*[1-exp(-d*x)]+e}.\cr\cr
The weighted curve fitting is performed using the Levenberg-Marquardt algorithm by J.J. More (available at \url{http://ftp.netlib.org/minpack}). 
When fitting a non-linear model, a series of \code{b(d)} values (i.e., the reciprocal of the saturating dose) are initialized randomly with a uniform distribution in the space \code{(0,upb)}, then other parameters are estimated with a Linear Algebra method. For a growth curve that shows high saturating dose, it is more reasonable to use a smaller \code{upb} value. For example, a double saturating exponential model with a \code{upb} value less than 0.05 or 0.01 usually works better. The fitting process will be performed repeatedly with \code{nstart} times to find an optimal estimate. User is advised to set argument \code{plot} to be TRUE if possible in order to see if the fitted growth curve is correct.
}

\value{
Return an list that contains following elements:
  \item{LMpars}{parameters and standard errors of the growth curve}
  \item{value}{the minimized chi-square value of the fitted growth curve}
  \item{fit.value}{observations and fitted values}
}

\references{
More JJ, 1978. "The Levenberg-Marquardt algorithm: implementation and theory," in Lecture Notes in Mathematics: Numerical Analysis, 
Springer-Verlag: Berlin. 105-116.
}

\note{
For a given model, this routine will return an error if any standard errors of parameters of the growth curve cannot be estimated by numerical approximation.\cr\cr
}

\seealso{
  \link{analyst}; \link{calED}; \link{fastED}
}

\examples{
 data(Signaldata) 
 dat<-analyst(Signaldata[[1]], Redose=c(80,160,240,320,0,80)*0.13, 
              model="line", origin=TRUE, plot=FALSE, nsim=100)
 fitGrowth(dat$Curvedata, model="exp", origin=FALSE, plot=TRUE)
}

\keyword{growth curve}