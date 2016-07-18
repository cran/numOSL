\name{fitGrowth}
\alias{fitGrowth}
\alias{fitGrowth.default}
\title{
Growth curve fitting
}
\description{
Fitting growth curves using the Levenberg-Marquardt algorithm.
}
\usage{
fitGrowth(Curvedata, model = "exp", origin = FALSE, 
          weight = TRUE, plot = TRUE)
}

\arguments{
  \item{Curvedata}{\link{matrix}(\bold{required}): a three-column matrix (i.e., regenerative doses, standardised regenerative dose signals, and measurement errors of standardised regenerative dose signals)}

  \item{model}{\link{character}(with default): model used for growth curve fitting, see details}

  \item{origin}{\link{logical}(optional): force the fitted growth curve to pass the origin or not}

  \item{weight}{\link{logical}(with default): perform growth curve fitting using a weighted or an un-weighted procedure, see details}

  \item{plot}{\link{logical}(with default): plot output or not}
}

\details{
In growth curve fitting using function \link{fitGrowth}, five models are available: \cr
\bold{(1)} \emph{"line"}: a linear model, \code{y=a*x+b};\cr
\bold{(2)} \emph{"exp"}: a single saturation exponential model, \code{y=a*[1-exp(-b*x)]+c};\cr
\bold{(3)} \emph{"lexp"}: a single saturation exponential plus linear model, \code{y=a*[1-exp(-b*x)]+c*x+d};\cr
\bold{(4)} \emph{"dexp"}: a double saturation exponential model, \code{y=a*[1-exp(-b*x)]+c*[1-exp(-d*x)]+e};\cr
\bold{(5)} \emph{"gok"}: a general order kinetic model (Guralnik et al., 2015), \code{y=a*[1-(1+b*c*x)^(-1/c)]+d}.\cr\cr
The fitting process is performed using the Levenberg-Marquardt algorithm (minpack: Fortran 90 source code by John Burkardt, freely available at \url{http://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.html}). If \code{weight=TRUE}, a weighted procedure will be performed through weighting each data point by its inverse variance. User is advised to set argument \code{plot=TRUE} if possible to visualize the quality of fit.
}

\value{
Return a list that contains the following elements:
  \item{LMpars}{optimized parameters for the fitted growth curve}
  \item{value}{minimized objective for the fitted growth curve}
  \item{fit.value}{observed standardised regenerative dose signals and fitted counterparts}
}

\references{
Guralnik B, Li B, Jain M, Chen R, Paris RB, Murray AS, Li SH, Pagonis V, Valla PG, Herman F, 2015. Radiation-induced growth and isothermal decay of infrared-stimulated luminescence from feldspar. Radiation Measurements, 81: 224-231. 

More JJ, 1978. "The Levenberg-Marquardt algorithm: implementation and theory" in Lecture Notes in Mathematics: Numerical Analysis, 
Springer-Verlag: Berlin. 105-116.
}

\note{
The model to be optimized should not be underdetermined. This means that the number of data points should exceed (or at least be equal to) the number of parameters. For a given model, the procedure will return an error if any standard errors of parameters cannot be estimated by numerical difference-approximation.\cr\cr
}

\seealso{
  \link{analyst}; \link{calED}; \link{fastED}; \link{calSARED}; \link{calSGCED}
}

\examples{
 Curvedata<-cbind(c(0, 18, 36, 54, 72, 0, 18),               
                  c(0.026, 1.55, 2.39, 3.46, 4.13, 0.023, 1.61),  
                  c(0.005, 0.11, 0.27, 0.22, 0.20, 0.008, 0.24)) 
 fitGrowth(Curvedata, model="gok", origin=FALSE, plot=TRUE)
}

\keyword{growth curve}
