\name{fastED}
\alias{fastED}
\alias{fastED.default}
\title{
Estimating a fast-component equivalent dose
}
\description{
A number of decay curves obtained through a single aliquot regenerative (SAR) protocol (Murray and Wintle, 2000) are decomposed and fast-component decaying signals are extracted (if possible) to build a dose response curve to estimate a fast-component equivalent dose.
}
\usage{
fastED(Sigdata, Redose, ncomp = 2, constant=TRUE,  
       control.args = list(), typ = "cw", nstart = 100, upb = 1, 
       ErrorMethod = c("mc", "sp"), model = NULL, origin = NULL)
}
\arguments{
  \item{Sigdata}{\link{data.frame}(\bold{required}): a number of decay curves obtained with the SAR method, the first column should be the time values}
  \item{Redose}{\link{vector}(\bold{required}): the regenerative dose values (redose1, redose2, etc)}
  \item{ncomp}{\link{numeric}(with default): the number of components used to decompose a decay curve, which should be an integer in the space [2,4]}
  \item{constant}{\link{logical}(with default): subtract a constant component when fitting the curve or not, see \link{decompc}}
  \item{control.args}{\link{list}(with default): some further arguments used in the differential evolution algorithm, see \link{decomp}}
  \item{typ}{\link{character}(with default): the type of an OSL decay curve, only type of "CW" is available currently}
  \item{nstart}{\link{numeric}(with default): the allowed maximum number of attempts that used to initialize parameters during the fitting of a growth curve, see \link{calED}}
  \item{upb}{\link{numeric}(with default): the upper limit on b value, see \link{calED}}
  \item{ErrorMethod}{\link{character}(with default): a method ("sp" or "mc") for assessing the standard error of an equivalent dose. Default ErrorMethod="mc". See \link{calED}}
  \item{model}{\link{character}(optional): the fitting model ("line", "exp", or "line+exp"). See \link{calED}}
  \item{origin}{\link{logical}(optional): force the fitting to pass the origin (x=0,y=0) or not, see \link{calED}}
}
\details{
Function fastED() is used to estimated a fast-component equivalent dose with data obtained by the SAR protocol (Murray and Wintle, 2000), only CW-OSL data can be analyzed currently. The routine trys to decompose a series of decay curves to a specified number of components, then the fast-component signal from each decay curve is extracted to build a dose response curve, from which a fast-component equivalent dose can be obtained.
}
\value{
Return a list containing following elements:
  \item{pars}{the estimated parameters and standard errors for each decay curve}
  \item{decayRate}{a summarized list for decay rates of decay curves (except the decay curve that corresponds to the 0-regenerative dose and its test dose)}
  \item{Curvedata}{data used to build the fast-component dose response curve, see \link{calED}}
  \item{Ltx}{the estimated standardised natural OSLs and their standard error, see \link{calED}}
  \item{model}{the model used for fitting the dose-response curve, which is chosen automatically between "line", "exp" and "line+exp" with the principle of the minimum of the sum of squared residuals if model=NULL}
  \item{LMpars}{the estimated characteristic parameters and standard errors of the growth curve}
  \item{value}{the sum of squared residuals for the fitted growth curve}
  \item{ED}{the estimated fast-component equivalent dose (and its standard error)}
  \item{RecyclingRatio}{the estimated fast-component recycling ratio}
  \item{Recuperation}{the estimated fast-component recuperation}
}
\references{
Murray, A.S., Wintle, A.G., 2000. Luminescence dating of quartz using improved single-aliquot regenerative-dose protocol. Radiation Measurements, 32, pp.57-73.

Li, S.H., Li, B., 2006. Dose measurement using the fast component of LM-OSL signals from quartz. Radiation Measurements, 41, pp.534-541.
}
\seealso{
  \link{calED}; \link{decomp}; \link{decompc}
}
\note{
The number of trapped electrons that corresponds to the largest decay rate will be regarded as the fast-component signal, which cannot always ensure that a pure fast-component signal can be extracted if an ultra-fast component appears. Hopefully the arguments (factor, f, cr, maxiter, tol) in "control.args" may be employed to achieve better performance for decay curves decomposition. For example, larger number of parents (NP) (controlled by "factor") may perform better at the cost of slowdowning the running speed. The author thank Sheng-hua Li and Geoff Duller for their helpful discussions concerning the fast-component equivalent dose calculation.
}
\examples{
 data(Signaldata)
 ### Use only the first two regenerative dose data. 
 fastED(Signaldata$cw[,-(8:15)],Redose=c(80,160)*0.13,
        ncomp=2,constant=TRUE,model="line",origin=TRUE)
 ### Use all regenerative dose data (except the zero-dose).
 # fastED(Signaldata$cw[,c(-12,-13)],Redose=c(80,160,240,320,80)*0.13,
 #        ncomp=3,constant=FALSE)
}
\keyword{fitting OSL decay curve}
\keyword{fast-component equivalent dose}
