\name{analyst}
\alias{analyst}
\alias{analyst.default}

\title{
Equivalent dose analysis
}

\description{
Determing an equivalent dose using data sets obtained from the single aliquot regenerative-dose (SAR) sequence. 
}

\usage{
analyst(Sigdata, Redose, sig.channel = NULL, back.channel = NULL, 
        mr = 0.01, typ = "cw", model = "exp", origin = FALSE, 
        ErrorMethod = "mc", nsim = 1000, weight = TRUE, plot = TRUE)
}

\arguments{
  \item{Sigdata}{\link{matrix}(\bold{required}): a series of decay curves stored in a column-matrix}

  \item{Redose}{\link{vector}(\bold{required}): regenerative dose values. Example, \code{Redose=c(1,2,3,4,0,1)}}

  \item{sig.channel}{\link{vector}(optional): channels used for signal integration, default \code{sig.channel=1:4}}

  \item{back.channel}{\link{vector}(optional): channels used for background subtraction, default \code{back.channel} is set equal to the last 20 channels}

  \item{mr}{\link{double}(with default): instrumental measurement error}

  \item{typ}{\link{character}(with default): type of an OSL decay curve, only CW-OSL decay curves can be analyzed currently}

  \item{model}{\link{character} (with default): model used for growth curve fitting, see function \link{fitGrowth} for available models}

  \item{origin}{\link{logical}(with default): force the fitted growth curve to pass the origin or not}

  \item{ErrorMethod}{\link{character}(with default): method used for assessing the standard error of an equivalent dose, see function \link{calED} for details}

  \item{nsim}{\link{integer}(with default): number of Monte Carlo simulations used for equivalent dose error estimation if \code{ErrorMethod="mc"}}

 \item{weight}{\link{logical}(with default): perform growth curve fitting using a weighted or an un-weighted procedure, see function \link{fitGrowth} for details}

  \item{plot}{\link{logical}(with default): plot output or not}
}

\details{
Function \link{analyst} is used for determing an equivalent dose using data sets obtained from the SAR sequence. It works in a manner similar to software \bold{\emph{Analyst}} written by Duller (2007). See function \link{fitGrowth} and \link{calED} for more details concerning growth curve fitting and equivalent dose calculation.
}

\value{
Return a list containing the following elements:
  \item{Curvedata}{data used for constructing the growth curve}
  \item{Ltx}{standardised natural OSL signal and its standard error}
  \item{LMpars}{optimized parameters for the fitted growth curve}
  \item{value}{minimized objective for the fitted growth curve}
  \item{ED}{calculated equivalent dose and its standard error}
  \item{RecyclingRatio}{recycling ratio and its standard error}
  \item{Recuperation}{recuperation and its standard error (in percent)}
}

\references{
Galbraith R, 2002. A note on the variance of a backround-corrected OSL count. Ancient TL, 20(2): 49-51.

Duller, G., 2007. Software Analyst, user mannual. Freely available at \url{http://www.nutech.dtu.dk/english/Products-and-Services/Dosimetry/Radiation-Measurement-Instruments/TL_OSL_reader/Software}.
}

\note{
Argument \code{Sigdata} is a matrix made up with a number of decay curves column by column\cr
(i.e., \emph{natural_signal},\emph{test_natural_signal}, \emph{Redose[1]_signal},\cr
\emph{test_Redose[1]_signal},....etc).\cr\cr
Variation in background counts is assumed to follow a Poisson distribution, \bold{Eqn.(3)} of Galbraith (2002) is used to estimate the standard error of OSL signals.\cr\cr
Note that the curve fitting procedure will return an error if any standard errors of parameters cannot be estimated by numerical difference-approximation.\cr\cr
}

\seealso{
 \link{calED}; \link{fastED}; \link{fitGrowth}
}

\examples{
 # load package "numOSL".
 library(numOSL)
 # load data.
 data(Signaldata)
 # Early background subtraction.
 analyst(Signaldata$cw[,-1L], Redose=c(80,160,240,320,0,80)*0.13,
         sig.channel=1:2, back.channel=3:8, model="exp")
}
\keyword{growth curve}
\keyword{equivalent dose}
