\name{fastED}
\alias{fastED}
\alias{fastED.default}
\title{
Fast-component equivalent dose calculation
}
\description{
Estimating a fast-component equivalent dose using decay curves obtained from the single aliquot regenerative-dose (SAR) sequence.
}
\usage{
fastED(Sigdata, Redose, ncomp = 2, constant = TRUE,  
       control.args = list(), typ = "cw", model = "exp", 
       origin = FALSE, ErrorMethod = "mc", 
       nsim = 1000, weight = TRUE)
}
\arguments{
  \item{Sigdata}{\link{matrix}(\bold{required}): a series of decay curves stored in a matrix column by column, the first column denotes stimulation time values}

  \item{Redose}{\link{vector}(\bold{required}): regenerative dose values. Example, \code{Redose=c(1,2,3,4,0,1)}}

  \item{ncomp}{\link{integer}(with default): number of decomposed components}

  \item{constant}{\link{logical}(with default): subtract a constant background or not}

  \item{control.args}{\link{list}(with default): arguments used in the differential evolution algorithm, see function \link{decomp} for details}

  \item{typ}{\link{character}(with default): type of an OSL decay curve, only CW-OSL decay curve can be analyzed currently}

  \item{model}{\link{character}(with default): model used for growth curve fitting, see function \link{fitGrowth} for available models}

  \item{origin}{\link{logical}(with default): force the fitted growth curve to pass the origin or not}

  \item{ErrorMethod}{\link{character}(with default): method used for assessing the standard error of an equivalent dose, see function \link{calED} for details}

  \item{nsim}{\link{integer}(with default): number of Monte Carlo simulations used for equivalent dose error estimation if \code{ErrorMethod="mc"}}

  \item{weight}{\link{character}(with default): perform growth curve fitting using a weighted or an un-weighted procedure, see function \link{fitGrowth} for details}

}

\details{
Function \link{fastED} is used to estimate a fast-component equivalent dose using data sets obtained from the SAR protocol (Murray and Wintle, 2000). The routine trys to decompose a series of decay curves to a specified number of components, then the numbers of trapped electrons from the fast-component will be used to construct the growth curve to estimate a fast-component equivalent dose. See function \link{decomp}, \link{fitGrowth}, and \link{calED} for more details concerning decay curve decomposition, growth curve fitting, and equivalent dose calculation.
}

\value{
Return a list containing the following elements:
  \item{decayPars}{estimated parameters of fitted decay curves}
  \item{Curvedata}{data used for constructing the fast-component growth curve}
  \item{Ltx}{standardised natural dose signal of the fast-component and its standard error}
  \item{LMpars}{optimizaed parameters for the fitted growth curve}
  \item{value}{minimized objective for the fitted growth curve}
  \item{fastED}{fast-component equivalent dose and its standard error}
  \item{RecyclingRatio}{fast-component recycling ratio and its standard error}
  \item{Recuperation}{fast-component recuperation and its standard error (in percent)}
}
\references{
Li SH, Li B, 2006. Dose measurement using the fast component of LM-OSL signals from quartz. Radiation Measurements, 41(5): 534-541.

Murray AS, Wintle AG, 2000. Luminescence dating of quartz using improved single-aliquot regenerative-dose protocol. Radiation Measurements, 32(1): 57-73.
}

\seealso{
  \link{analyst}; \link{calED}; \link{decomp}; \link{fitGrowth}
}

\note{
Argument \code{Sigdata} is a matrix made up with stimulation time values in the first column and a \cr
number of decay curves (i.e., \emph{natural_signal}, \emph{test_natural_signal}, \emph{Redose[1]_signal}, \cr
\emph{test_Redose[1]_signal},....etc) in the rest columns.\cr\cr 
The number of trapped electrons that corresponds to the largest decay rate will be regarded as the fast-component signal, which cannot always ensure that a pure fast-component signal be extracted if an ultra-fast decaying component appears.\cr\cr
The authors thank Professor Sheng-Hua Li and Professor Geoff Duller for their helpful discussions concerning fast-component equivalent dose calculation.
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
\keyword{decay curve}
\keyword{fast-component equivalent dose}
\keyword{growth curve}
