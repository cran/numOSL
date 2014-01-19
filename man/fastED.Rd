\name{fastED}
\alias{fastED}
\alias{fastED.default}
\title{
Estimating a fast-component equivalent dose in OSL dating
}
\description{
With the data obtained by the single aliquot regenerative (SAR) protocol (Murray and Wintle, 2000), a series decay curves are decomposed and signals of fast-component are extracted (if possible) to build a dose response curve to approximate a fast-component equivalent dose.
}
\usage{
fastED(Sigdata, Redose, ncomp = 2, constant=TRUE,  
       control.args = list(), typ = "cw", nstart = 100,  
       upb = 1, ErrorMethod = c("mc", "sp"), origin = NULL)
}
\arguments{
  \item{Sigdata}{\link{data.frame}(\bold{required}): a series of decay curves obtained with SAR method, the first column of which must be time values and at least a total of five columns must be provided
}
  \item{Redose}{\link{vector}(\bold{required}): regenerative dose values, such as redose1, redose2, etc
}
  \item{ncomp}{\link{numeric}(with default): number of components used to decompose each decay curve, which must be an integer in the space [2,4]
}
  \item{constant}{\link{logical}(with default): subtract a constant in curve decomposition or not
}
  \item{control.args}{\link{list}(with default): further arguments used in differential evolution algorithm, see \link{decomp}
}
  \item{typ}{\link{character}(with default): type of the OSL decay curve, only type of "CW" is available currently
}
  \item{nstart}{\link{numeric}(with default): maximum number of attempts that used to initialize parameters in dose response curve fitting, see \link{calED}
}
  \item{upb}{\link{numeric}(with default): upper boundary of b value, initial b values will be generated uniformly from (0, upb), see \link{calED}
}
  \item{ErrorMethod}{\link{character}(with default): method ("sp" or "mc") for estimating the standard error of equivalent dose. Default ErrorMethod="mc". See \link{calED}
}
  \item{origin}{\link{logical}(optional): whether force the fitting to pass the origin (x=0,y=0) or not, see \link{calED}
}
}
\details{
Function fastED() is used to estimated a fast-component equivalent dose with data obtained by the SAR protocol (Murray and Wintle, 2000), only CW-OSL data can be analyzed currently. The routine trys to decompose a series of decay curves to a specified number of components, then the fast-component signal from each decay curve will be extracted to build the dose response curve, from which a fast-component equivalent dose could be obtained.
}
\value{
Return a list containing following elements:\cr
  \item{pars}{estimated parameters and associated standard errors for each decay curve}
  \item{decayRate}{summarized list about all decay rates (if zero regenerative dose appears, the corresponded decay rates will be precluded before calculating), such as the uniformity of decay rates between each component and each decay curve; the average and associated standard error of the fast-component decay rates}
  \item{Curvedata}{data used to build the fast-component dose response curve, see \link{calED}}
  \item{Lxt}{standardized natural signal and its standard error, see \link{calED}}
  \item{model}{model used for fitting the dose-response curve, which is chosen automatically from "line", "exp" or "line+exp" during the calculation with the principle of the minimum of the sum of squared residuals}
  \item{LMpars}{estimated characteristic parameters and standard errors of the dose-response curve}
  \item{residual}{sum of the squared residuals of the fitted dose response curve}
  \item{ED}{estimated fast-component equivalent dose (and corresponded standard error)}
  \item{RecyclingRatio}{estimated recycling ratio}
  \item{Recuperation}{estimated recuperation}
}
\references{
Murray, A.S., Wintle, A.G., 2000. Luminescence dating of quartz using improved single-aliquot regenerative-dose protocol. Radiation Measurements, 32, pp.57-73.

Li, S.H., Li, B., 2006. Dose measurement using the fast component of LM-OSL signals from quartz. Radiation Measurements, 41, pp.534-541.
}
\seealso{
  \link{calED}; \link{decomp}; \link{decompc}
}
\note{
The number of trapped electrons that corresponds to the largest decay rate will be regarded as the fast-component signal, which can not always ensure a pure fast-component signal to be extracted if ultra-fast component appears. To extract purer fast-component signal, the background from experimental data should be subtracted before the analysis if possible, or it might be more reasonable to keep the argument "constant" being active. Hopefully the arguments (factor, f, cr, maxiter, tol) in  "control.args" could be employed to achieve better performance in decay curves decomposition. For example, larger number of parents (NP) (controlled by "factor") might performs better at the cost of slowdown of the running speed. The author thank Sheng-hua Li and Geoff Duller for their helpful discussions concerning fast-component equivalent dose calculation.
}
\examples{
# Loading signal data
 data(Signaldata)
# Calculating a fast-component dose with three components after removing the 
# decay curves that correspond to zero-redose (12th) and associated test-dose (13th)
 fastED(Signaldata$cw[,c(-12,-13)],Redose=c(80,160,240,320,80)*0.13,
        ncomp=3,constant=FALSE)
# Trying with two component plus a constant
 fastED(Signaldata$cw,Redose=c(80,160,240,320,0,80)*0.13,ncomp=2)
}
\keyword{CW-OSL decay curve}
\keyword{fast-component De}
