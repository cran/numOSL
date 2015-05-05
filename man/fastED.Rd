\name{fastED}
\alias{fastED}
\alias{fastED.default}
\title{
Fast-component equivalent dose calculation
}
\description{
Estimating a fast-component equivalent dose using decay curves obtained from a SAR sequence.
}
\usage{
fastED(Sigdata, Redose, ncomp = 2, constant = TRUE,  
       control.args = list(), typ = "cw", nstart = 100, 
       upb = 0.5, ErrorMethod = c("mc","sp"), nsim = 1000,
       model = NULL, origin = NULL, weight = c("g","d","b","n"))
}
\arguments{
  \item{Sigdata}{\link{matrix}(\bold{required}): a series of decay curves stored column by column, the structure is different from that in function \link{analyst}, here the first column represents stimulation time values}

  \item{Redose}{\link{vector}(\bold{required}): regenerative dose values. 
Example, \code{Redose=c(1,2,3,4,0,1)}}

  \item{ncomp}{\link{integer}(with default): number of components to be decomposed}

  \item{constant}{\link{logical}(with default): subtract a constant component or not}

  \item{control.args}{\link{list}(with default): arguments used by the differential evolution algorithm, see function \link{decomp} for details}

  \item{typ}{\link{character}(with default): type of an OSL decay curve, only CW-OSL decay curve can be analyzed currently}

  \item{nstart}{\link{integer}(with default): maximum number of trials used for initializing the non-linear growth curve, see function \link{fitGrowth} for details}

  \item{upb}{\link{double}(with default): upper limit on \code{b} value (i.e., the reciprocal of the saturating dose), see function \link{fitGrowth} for details}

  \item{ErrorMethod}{\link{character}(with default): a method used for assessing the standard error of an equivalent dose, default \code{ErrorMethod="mc"}, see function \link{calED} for details}

  \item{nsim}{\link{integer}(with default): maximum number of simulations if a Monte Carlo method is applied to assessing the standard error of an equivalent dose}

  \item{model}{\link{character}(optional): model used for fitting the growth curve, see function \link{fitGrowth} for available models. The procedure will find an optimal model between \code{"line"} and \code{"exp"} automatically if \code{model=NULL}}

  \item{origin}{\link{logical}(optional): force the fast-component growth curve to pass the origin or not. The procedure will choose whether pass the origin or not automatically if \code{origin=NULL}}

  \item{weight}{\link{character}(with default): type of curves to be fitted by a weighted procedure, \code{"g"} means applying a weighted procedure to fitting the growth curve, \code{"d"} means applying a weigthed procedure to fitting the decay curves, \code{"b"} means using a weighted procedure to fit both types of curves, and \code{"n"} means applying an un-weighted procedure to fitting both types of curves. Default \code{weight="g"}}

}

\details{
Function \emph{fastED} is used to estimate a fast-component equivalent dose using a data set obtained from a SAR protocol (Murray and Wintle, 2000). The routine trys to decompose a series of decay curves to a specified number of components, then the numbers of trapped electrons from the fast-component are used to construct a growth curve to estimate a fast-component equivalent dose. See function \link{decomp}, \link{fitGrowth}, and \link{calED} for more details concerning decay curve decomposition, growth curve fitting, and equivalent dose calculation.
}

\value{
Return a list containing the following elements:
  \item{pars}{estimated decay rates and trapped electrons for each decay curve}
  \item{decayRateMat}{decay rates for each decay curve stored in a matrix}
  \item{Curvedata}{data used for constructing the fast-component growth curve}
  \item{Ltx}{sensitivity-corrected natural OSL for the fast-component and associated standard error}
  \item{model}{fitting model of the growth curve}
  \item{LMpars}{optimizaed parameters of the growth curve}
  \item{value}{minimized objective of the growth curve}
  \item{ED}{fast-component equivalent dose and its standard error}
  \item{RecyclingRatio}{fast-component recycling ratio}
  \item{Recuperation}{fast-component recuperation}
}
\references{
Li SH, Li B, 2006. Dose measurement using the fast component of LM-OSL signals from quartz. Radiation Measurements, 41(5): 534-541.

Murray AS, Wintle AG, 2000. Luminescence dating of quartz using improved single-aliquot regenerative-dose protocol. Radiation Measurements, 32(1): 57-73.
}

\seealso{
  \link{analyst}; \link{calED}; \link{decomp}; \link{fitGrowth}
}

\note{
Argument \code{Sigdata} is a column matrix filled with stimulation time values in the first column and a \cr
number of decay curves (i.e., \emph{natural_signal}, \emph{test_natural_signal}, \emph{Redose[1]_signal}, \cr
\emph{test_Redose[1]_signal},....etc) stored in the rest columns.\cr\cr 
The number of trapped electrons that corresponds to the largest decay rate will be regarded as the fast-component signal, which cannot always ensure that a pure fast-component signal be extracted if an ultra-fast decaying component appears.\cr\cr
The author thank Professor Sheng-Hua Li and Professor Geoff Duller for their helpful discussions concerning fast-component equivalent dose calculation.
}
\examples{
 data(Signaldata)
 # Use only the first two regenerative dose values. 
 fastED(Signaldata$cw[,-(8:15)],Redose=c(80,160)*0.13,model="line",
        control.args=list(factor=5,maxiter=50),origin=TRUE)
 # Use all regenerative dose values (except the zero-dose).
 # fastED(Signaldata$cw[,c(-12,-13)],Redose=c(80,160,240,320,80)*0.13,
 #        ncomp=3,constant=FALSE)
}
\keyword{OSL decay curve}
\keyword{fast-component equivalent dose}
\keyword{growth curve}
