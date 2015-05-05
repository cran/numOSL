\name{calED}
\alias{calED}
\alias{calED.default}

\title{
Equivalent dose calculation and error assessment
}

\description{
Calculating an equivalent dose value and assessing its standard error using methods described in Duller (2007).
}

\usage{
calED(Curvedata, Ltx, model = c("line","exp","lexp","dexp"), 
      origin = FALSE, nstart = 100, upb = 0.5, ErrorMethod = c("mc","sp"),
      nsim = 1000, weight = TRUE, plot = TRUE, outfile = NULL)
}

\arguments{
  \item{Curvedata}{\link{matrix}(\bold{required}): three columns, i.e., regenerative doses, sensitivity-corrected OSLs, and standard errors of sensitivity-corrected OSLs}

  \item{Ltx}{\link{vector} or \link{matrix}(\bold{required}): sensitivity-corrected natural OSL(s) and associated standard error(s)}

  \item{model}{\link{character}(with default): fitting model for a growth curve, default \code{model="line"}, see \link{fitGrowth} for available models}

  \item{origin}{\link{logical}(optional): force the growth curve to pass the origin or not}

  \item{nstart}{\link{integer}(with default): maximum number of attempts for initializing the non-linear growth curve, see function \link{fitGrowth} for details}

  \item{upb}{\link{double}(with default): upper limit on \code{b} value (i.e., the reciprocal of the saturating dose), initial \code{b} value is generated uniformly from the space \code{(0,upb)}, see function \link{fitGrowth} for details}

  \item{ErrorMethod}{\link{character}(with default): method used for assessing the standard error of an equivalent dose, default \code{ErrorMethod="mc"}, see details}

  \item{nsim}{\link{integer}(with default): maximum number of simulations if a Monte Carlo method is applied to assessing the standard error of an equivalent dose}

  \item{weight}{\link{logical}(with default): fit the growth curve using a weighted or an un-weighted procedure, see function \link{fitGrowth} for details}

  \item{plot}{\link{logical}(with default): draw a plot or not}

  \item{outfile}{\link{character}(with default): if specified, simulated equivalent dose values will be written to a file named \code{"outfile"} in \code{".csv"} format and saved to the current work directory}
}

\details{
Function \emph{calED} can be used for calculating a single equivalent dose value (Duller, 2007) (if \code{Ltx} is a two-element \link{vector}) or a series of equivalent dose values using the "Standardised Growth Curve" method (i.e. "SGC" used by Roberts and Duller, 2004; or "gSGC" proposed by Li et al., 2015) (if \code{Ltx} is a row \link{matrix}).\cr\cr
Standard error of an equivalent dose can be assessd using two methods outlined by Duller (2007): \cr
\bold{(1)} \code{"sp"}: a simple transformation from error of a sensitivity-corrected OSL (i.e., \code{s(Lx/Tx)}) to error of an equivalent dose (i.e., \code{s(ED)}); \cr
\bold{(2)} \code{"mc"}: a Monte Carlo method (a "parametric bootstrap" protocol).\cr\cr
Interpolation is performed using a combination of golden section search and successive parabolic interpolation (using freely available subroutine at \url{http://www.netlib.org/fmm/fmin.f}). See function \link{fitGrowth} for more details on growth curve fitting.
}

\value{
Return a list that contains the following elements:
  \item{LMpars}{parameters of the growth curve}
  \item{value}{minimized objective of the growth curve}
  \item{fit.value}{observations and fitted values}
  \item{ED}{equivalent dose(s) and standard error(s)}
}

\references{
Duller GAT, 2007. Assessing the error on equivalent dose estimates derived from single aliquot regenerative dose measurements. Ancient TL, 25(1): 15-24.

Li B, Roberts RG, Jacobs Z, Li SH, 2015. Potential of establishing a "global standardised growth curve" (gSGC) for optical dating of quartz from sediments. Quaternary Geochronology, 27: 94-104.

Roberts HM, Duller GAT, 2004. Standardised growth curves for optical dating of sediment using multiple-grain aliquots. Radiation Measurements, 38(2): 241-252.
}

\seealso{
  \link{analyst}; \link{fastED}; \link{fitGrowth}
}


\note{
Interpolation is performed within the range of regenerative dose values, which means that the sensitivity-corrected natural OSL (i.e., \code{Ltx}) used for calculating an equivalent dose should not excced the maximum sensitivity-corrected regenerative OSL used for constructing the growth curve if a non-linear model is applied.\cr\cr
For a given model, the curve fitting procedure will return an error if any standard error of a parameter cannot be estimated by numerical difference-approximation. Function \emph{sgcED} in previous versions was bundled to function \link{calED}.
}

\examples{
  # In SGC method application, if aliquots have different test 
  # dose values, their sensitivity-corrected OSL (Lx/Tx) need be 
  # normalized using their test dose Dt so as to reduce the between-aliquot
  # scatter. But in the following example, normalization is unnecessary, 
  # as the two aliquots share a common test dose. 
  data1<-cbind(c(0, 18, 36, 54, 72, 0, 18),               
               c(0.026, 1.55, 2.39, 3.46, 4.13, 0.023, 1.61),  
               c(0.005, 0.11, 0.27, 0.22, 0.20, 0.008, 0.24)) 
  data2<-cbind(c(0, 18, 36, 54, 72, 0, 18),                
               c(0.021, 1.47, 2.51, 3.37, 4.31, 0.028, 1.52),  
               c(0.003, 0.09, 0.35, 0.21, 0.29, 0.01, 0.17))  
  data <- rbind(data1, data2)  
  # data3, data4, ... can be added in the same way.                         
  Ltx<-cbind(c(0.5, 1.0, 1.8, 2.3, 2.8, 3.1, 3.6, 4.0),         
             c(0.02, 0.08, 0.11, 0.12, 0.23, 0.31, 0.33, 0.51)) 
  calED(Curvedata=data, Ltx=Ltx, model= "exp", origin=FALSE)
}
\keyword{growth curve}
\keyword{equivalent dose}
