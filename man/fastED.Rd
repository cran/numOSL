\name{fastED}
\alias{fastED}
\alias{fastED.default}
\title{
Estimating a fast-component equivalent dose
}
\description{
Estimating a fast-component equivalent dose using decay curves from a SAR sequence.
}
\usage{
fastED(Sigdata, Redose, ncomp = 2, constant = TRUE,  
       control.args = list(), typ = "cw", nstart = 100, upb = 0.5, 
       ErrorMethod = c("mc","sp"), model = NULL, origin = NULL)
}
\arguments{
  \item{Sigdata}{\link{data.frame}(\bold{required}): a number of decay curves, the first column should be the time values}

  \item{Redose}{\link{vector}(\bold{required}): the regenerative dose values. Example, \code{Redose=c(1,2,3,4,0,1)}}

  \item{ncomp}{\link{numeric}(with default): the number of components to be decomposed}

  \item{constant}{\link{logical}(with default): subtract a constant component when fitting a decay curve or not}

  \item{control.args}{\link{list}(with default): arguments used by the differential evolution algorithm, see function \link{decomp} for details}

  \item{typ}{\link{character}(with default): the type of an OSL decay curve, only type of \code{"CW"} is available currently}

  \item{nstart}{\link{numeric}(with default): the maximum number of attempts of initializing the non-linear growth curve fitting, see function \link{fitGrowth} for details}

  \item{upb}{\link{numeric}(with default): the upper limit on \code{b} value, see function \link{fitGrowth} for details}

  \item{ErrorMethod}{\link{character}(with default): method used for assessing the standard error of an equivalent dose, default \code{ErrorMethod="mc"}, see function \link{calED} for details}

  \item{model}{\link{character}(optional): a model used for fitting the growth curve, see function \link{fitGrowth} for available models. It will find a model automatically if \code{model=NULL}}

  \item{origin}{\link{logical}(optional): force the fitting to pass the origin or not. It will choose whether pass the origin or not automatically if \code{origin=NULL}}
}

\details{
Function \emph{fastED} is used to estimated a fast-component equivalent dose using a dataset from the SAR protocol (Murray and Wintle, 2000). The routine trys to decompose a series of decay curves to a specified number of components, then the numbers of trapped electrons from the fast-component are used to construct a growth curve to estimate a fast-component equivalent dose.
}

\value{
Return a list containing following elements:
  \item{pars}{the decay rates and trapped electrons for a number of decay curves}
  \item{decayRateMat}{the decay rates for a number of decay curves (a matrix)}
  \item{Curvedata}{data used to construct the fast-component growth curve}
  \item{Ltx}{the natural fast-component standardised OSL and its standard error}
  \item{model}{the fitting  model}
  \item{LMpars}{the parameters of the growth curve}
  \item{value}{the minimized chi-square value of the fitted growth curve}
  \item{ED}{the fast-component equivalent dose and its standard error}
  \item{RecyclingRatio}{the fast-component recycling ratio}
  \item{Recuperation}{the fast-component recuperation}
}
\references{
Li SH, Li B, 2006. Dose measurement using the fast component of LM-OSL signals from quartz. Radiation Measurements, 41: 534-541.

Murray AS, Wintle AG, 2000. Luminescence dating of quartz using improved single-aliquot regenerative-dose protocol. Radiation Measurements, 32: 57-73.
}

\seealso{
  \link{analyst}; \link{calED}; \link{decomp}; \link{fitGrowth}
}

\note{
The number of trapped electrons that corresponds to the largest decay rate will be regarded as the fast-component signal, which cannot always ensure that a pure fast-component signal can be extracted if an ultra-fast component appears. If \code{ErrorMethod="mc"}, then the standard error of an ED value will be assessed using \code{1000} Monte Carlo iterations (i.e., \code{nsim=1000}).\cr\cr
The author thank Sheng-Hua Li and Geoff Duller for their helpful discussions concerning the fast-component equivalent dose calculation.
}
\examples{
 data(Signaldata)
 ### Use only the first two regenerative dose data. 
 fastED(Signaldata$cw[,-(8:15)],Redose=c(80,160)*0.13,
        control.args=list(factor=5),model="line",origin=TRUE)
 ### Use all regenerative dose data (except the zero-dose).
 # fastED(Signaldata$cw[,c(-12,-13)],Redose=c(80,160,240,320,80)*0.13,
 #        ncomp=3,constant=FALSE)
}
\keyword{OSL decay curve}
\keyword{fast-component equivalent dose}
\keyword{growth curve}
