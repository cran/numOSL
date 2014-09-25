\name{calED}
\alias{calED}
\alias{calED.default}

\title{
Calculating an equivalent dose value and estimating its standard error
}

\description{
Calculating an equivalent dose value and assessing its standard error using methods described in Duller (2007).
}

\usage{
calED(Curvedata, Ltx, model = c("line","exp","lexp","dexp"), 
      origin = FALSE, nstart = 100, upb = 0.5, ErrorMethod = c("mc","sp"),
      nsim = 1000, plot = TRUE, outfile = NULL)
}

\arguments{
  \item{Curvedata}{\link{matrix}(\bold{required}): three columns, i.e., regenerative doses, standardised OSL signals, and standard errors of 
standardised OSL signals}

  \item{Ltx}{\link{vector} or \link{matrix}(\bold{required}): natural standardised OSL signal(s) and standard error(s), see details}

  \item{model}{\link{character}(with default): a fitting model, default \code{model="line"}, see \link{fitGrowth} for available models}

  \item{origin}{\link{logical}(optional): force the fitting to pass the origin or not}

  \item{nstart}{\link{numeric}(with default): the maximum number of attempts of initializing the non-linear growth curve fitting, see function \link{fitGrowth} for details}

  \item{upb}{\link{numeric}(with default): the upper limit on \code{b} value (i.e., the reciprocal of the saturating dose), initial \code{b} value is generated uniformly from the space \code{(0,upb)}, see function \link{fitGrowth} for details}

  \item{ErrorMethod}{\link{character}(with default): method used for assessing the standard error of an equivalent dose, default ErrorMethod="mc", see details}

  \item{nsim}{\link{numeric}(with default): the allowed maximum number of simulations when applying a Monte Carlo method to assess the standard error of an equivalent dose}

  \item{plot}{\link{logical}(with default): draw a plot or not}

  \item{outfile}{\link{character}(with default): if specified, simulated ED values will be written to a file named \code{"outfile"} in \code{".csv"} format and saved to the current work directory}
}

\details{
This function can be used for calculating a single ED value (if \code{Ltx} is a two-element \link{vector}) or a number of ED values using the "Standardised Growth Curve" method (if \code{Ltx} is a two-column multi-row \link{matrix}).\cr\cr
Standard error of an equivalent dose can be assessd using two method outlined by Duller (2007): \cr
\bold{(a)} \code{"sp"}: a simple transformation of \code{s(Lx/Tx)} to \code{s(ED)}; \cr
\bold{(b)} \code{"mc"}: a Monte Carlo method \cr\cr
Interpolation is performed using a combination of golden section search and successive parabolic interpolation (available at \url{http://www.netlib.org/fmm/fmin.f}).\cr\cr
See function \link{fitGrowth} for more details.
}

\value{
Return an list that contains following elements:
  \item{LMpars}{parameters and standard errors of the growth curve}
  \item{value}{the minimized chi-square value of the fitted growth curve}
  \item{fit.value}{observations and fitted values}
  \item{ED}{equivalent dose(s) and standard error(s)}
}

\references{
Duller GAT, 2007. Assessing the error on equivalent dose estimates derived from single aliquot regenerative dose measurements. 
Ancient TL, 25: 15-24.

Roberts HM, Duller GAT, 2004. Standardised growth curves for optical dating of sediment using multiple-grain aliquots. Radiation Measurements, 38: 241-252.
}

\seealso{
  \link{analyst}; \link{fastED}; \link{fitGrowth}
}


\note{
Interpolation is performed within the maximum regenerative dose, so the value of \code{Ltx} that used to calculate an equivalent dose cannot be larger than the maximum standardised OSL that is used to construct the growth curve when applying a non-linear model.\cr\cr
For a given model, this routine will return an error if any standard errors of parameters of the growth curve cannot be estimated by numerical approximation.\cr\cr
Function \emph{sgcED} in previous version was bundled to function \emph{calED}.
}

\examples{
  # Calculating a number of ED values (a SGC method).
  ltx1<-c(0.032,1.61,2.55,3.21,3.87,0.031,1.55) # Lx/Tx for the first aliquot
  ltx2<-c(0.025,1.44,2.47,3.35,4.17,0.033,1.47) # Lx/Tx for the second aliquot
  ltx3<-c(0.027,1.51,2.68,3.52,4.41,0.021,1.39) # Lx/Tx for the third aliquot
  ltx4<-c(0.018,1.71,2.28,3.81,4.03,0.017,1.62) # Lx/Tx for the four aliquot
  ltx5<-c(0.026,1.49,1.99,3.43,4.17,0.015,2.01) # Lx/Tx for the five aliquot
  ltx<-cbind(ltx1,ltx2,ltx3,ltx4,ltx5)
  ltx<-cbind(apply(ltx,MARGIN=1,mean),
             apply(ltx,MARGIN=1,sd))  # means and standard deviations
  redose<-c(0,18,36,54,72,0,18)       # the same ReDose for the five aliquots
  Curvedata<-data.frame(redose,ltx)
  Ltx<-cbind(c(0.5,1.0,1.8,2.3,2.8,3.1,3.6,4.0),rep(0.1,8))         
  calED(Curvedata,Ltx,model="exp",origin=FALSE) # fitting y=a*(1-exp(-b*x))+c
}
\keyword{growth curve}
\keyword{equivalent dose}
