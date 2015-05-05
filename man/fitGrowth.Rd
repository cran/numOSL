\name{fitGrowth}
\alias{fitGrowth}
\alias{fitGrowth.default}
\title{
Growth curve fitting
}
\description{
Fitting a growth curve using the Levenberg-Marquardt algorithm.
}
\usage{
fitGrowth(Curvedata, model = c("line","exp","lexp","dexp"), 
          origin = FALSE, nstart = 100, upb = 0.5, 
          weight = TRUE, plot = TRUE)
}

\arguments{
  \item{Curvedata}{\link{matrix}(\bold{required}): three columns (regenerative doses, sensitivity-corrected OSLs, and standard errors of sensitivity-corrected OSLs)}

  \item{model}{\link{character}(with default): a fitting model, default \code{model="line"}, see details}

  \item{origin}{\link{logical}(optional): force the growth curve to pass the origin or not}

  \item{nstart}{\link{integer}(with default): maximum number of trials used for initializing the non-linear growth curve, see details}

  \item{upb}{\link{double}(with default): upper limit on b value (i.e., the reciprocal of the saturating dose), initial b value is generated uniformly from the space \code{(0,upb)}, see details}

  \item{weight}{\link{logical}(with default): fit the growth curve using a weighted or an un-weighted procedure, see details}

  \item{plot}{\link{logical}(with default): draw a plot or not}
}

\details{
In growth curve fitting, four models are available: \cr
\bold{(1)} \emph{"line"}: a linear model, \code{y=a*x+b};\cr
\bold{(2)} \emph{"exp"}: an exponential model, \code{y=a*[1-exp(-b*x)]+c};\cr
\bold{(3)} \emph{"lexp"}: an exponential plus linear model, \code{y=a*[1-exp(-b*x)]+c*x+d};\cr
\bold{(4)} \emph{"dexp"}: a double exponential model, \code{y=a*[1-exp(-b*x)]+c*[1-exp(-d*x)]+e}.\cr\cr
The fitting process is performed using the Levenberg-Marquardt algorithm (minpack: Fortran 90 version by John Burkardt, freely available at \url{http://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.html}). If \code{weight=TRUE}, a weighted procedure will be performed through weighting each data point by its inverse variance. In non-linear growth curve fitting, a series of random numbers (i.e., \code{b} in \emph{"exp"} and \emph{"lexp"}, or \code{c(b, d)} in \emph{"dexp"}) are initialized with a uniform generator from the space \code{(0,upb)}, then other parameters (i.e., \code{c(a, c)} in \emph{"exp"}, \code{c(a, c, d)} in \emph{"lexp"}, or \code{c(a, c, e)} in \emph{"dexp"}) are estimated with a Linear Algebra method. For a growth curve that shows high saturating dose, it is more reasonable to use a small \code{upb} value. For example, a double saturating exponential model with a \code{upb} value less than 0.05 or 0.01 usually works better. The fitting process will be performed repeatedly with \code{nstart} times to find an optimal estimate. User is advised to set argument \code{plot=TRUE} if possible to visualize the fitting quality.
}

\value{
Return a list that contains the following elements:
  \item{LMpars}{parameters of the growth curve}
  \item{value}{minimized objective of the growth curve}
  \item{fit.value}{observations and fitted values}
}

\references{
More JJ, 1978. "The Levenberg-Marquardt algorithm: implementation and theory" in Lecture Notes in Mathematics: Numerical Analysis, 
Springer-Verlag: Berlin. 105-116.
}

\note{
The model to be optimized should not be underdetermined, which means that the number of data points should exceed (or at least be equal to) the number of parameters. For a given model, the procedure will return an error if any standard error of a parameter cannot be estimated by numerical difference-approximation.\cr\cr
}

\seealso{
  \link{analyst}; \link{calED}; \link{fastED}
}

\examples{
 Curvedata<-cbind(c(0, 18, 36, 54, 72, 0, 18),               
                  c(0.026, 1.55, 2.39, 3.46, 4.13, 0.023, 1.61),  
                  c(0.005, 0.11, 0.27, 0.22, 0.20, 0.008, 0.24)) 
 fitGrowth(Curvedata, model="exp", origin=FALSE, plot=TRUE)
}

\keyword{growth curve}
