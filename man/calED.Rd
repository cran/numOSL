\name{calED}
\alias{calED}
\alias{calED.default}

\title{
Equivalent dose calculation and error assessment
}

\description{
Calculating an equivalent dose and assessing its standard error.
}

\usage{
calED(Curvedata, Ltx, model = "exp", origin = FALSE, 
      ErrorMethod = "mc", nsim = 1000, weight = TRUE, 
      plot = TRUE)
}

\arguments{
  \item{Curvedata}{\link{matrix}(\bold{required}): a three-column matrix (i.e., regenerative doses, standardised regenerative dose signals, and measurement errors of standardised regenerative dose signals)}

  \item{Ltx}{\link{vector}(\bold{required}): a two-element vector consisting of standardised natural dose signal and its measurement error}

  \item{model}{\link{character}(with default): model used for growth curve fitting, see \link{fitGrowth} for available models}

  \item{origin}{\link{logical}(with default): force the fitted growth curve to pass the origin or not}

  \item{ErrorMethod}{\link{character}(with default): method used for assessing the standard error of an equivalent dose, see details}

  \item{nsim}{\link{integer}(with default): number of Monte Carlo simulations used for equivalent dose error estimation if \code{ErrorMethod="mc"}}

  \item{weight}{\link{logical}(with default): perform growth curve fitting using a weighted or an un-weighted procedure, see function \link{fitGrowth} for details}

  \item{plot}{\link{logical}(with default): plot output or not}
}

\details{
Function \link{calED} is used for calculating an equivalent dose and assessing its standard error.\cr\cr
The standard errors of an equivalent dose can be assessd using two methods outlined by Duller (2007): \cr
\bold{(1)} \code{"sp"}: a simple transformation of the standard error of the standardised natural dose signal to the standard error of an equivalent dose;\cr
\bold{(2)} \code{"mc"}: a Monte Carlo method (also called a "parametric bootstrap" protocol).\cr\cr
Interpolation is performed using a combination of golden section search and successive parabolic interpolation (\bold{\emph{R}} function \link{optimize} in package \bold{\emph{stats}}) (freely available Fortran 77 source code at \url{http://www.netlib.org/fmm/fmin.f}). See function \link{fitGrowth} for more details on growth curve fitting.
}

\value{
Return a list that contains the following elements:
  \item{mcED}{simulated random equivalent doses if \code{ErrorMethod="mc"}}
  \item{LMpars}{optimized parameters for the fitted growth curve}
  \item{value}{minimized objective for the fitted growth curve}
  \item{fit.value}{observed standardised regenerative dose signals and fitted counterparts}
  \item{ED}{calculated equivalent dose and standard error}
}

\references{
Duller GAT, 2007. Assessing the error on equivalent dose estimates derived from single aliquot regenerative dose measurements. Ancient TL, 25(1): 15-24.

\bold{Further reading}

Galbraith RF, Roberts RG, 2012. Statistical aspects of equivalent dose and error calculation and display in OSL dating: an overview and some recommendations. Quaternary Geochronology, 11: 1-27.
}

\seealso{
  \link{analyst}; \link{calSARED}; \link{calSGCED}; \link{fastED}; \link{fitGrowth}
}


\note{
For a given model, the curve fitting procedure will return an error if any standard errors of parameters cannot be estimated by numerical difference-approximation. Function \code{sgcED} in previous versions was bundled to function \link{calSGCED}.
}

\examples{
  Curvedata<-cbind(c(0, 18, 36, 54, 72, 0, 18),               
                   c(0.026, 1.55, 2.39, 3.46, 4.13, 0.023, 1.61),  
                   c(0.005, 0.11, 0.27, 0.22, 0.20, 0.008, 0.24))                         
  Ltx<-c(3.1,0.31)
  calED(Curvedata, Ltx, model="exp", origin=FALSE)
}
\keyword{growth curve}
\keyword{equivalent dose}
