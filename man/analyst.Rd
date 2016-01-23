\name{analyst}
\alias{analyst}
\alias{analyst.default}

\title{
Equivalent dose determination
}

\description{
Calculating an equivalent dose using dataset obtained from a SAR sequence. 
}

\usage{
analyst(Sigdata, Redose, sig.channel = NULL, back.channel = NULL, 
        mr = 0.01, typ = "cw", nstart = 100, upb = 0.5, ErrorMethod = c("mc","sp"), 
        nsim = 1000, weight = TRUE, plot = TRUE, model = NULL, origin = NULL)
}

\arguments{
  \item{Sigdata}{\link{matrix}(\bold{required}): a series of decay curves stored column by column}

  \item{Redose}{\link{vector}(\bold{required}): regenerative dose values. Example, \code{Redose=c(1,2,3,4,0,1)}}

  \item{sig.channel}{\link{vector}(with default): channels used for signal integration, default \code{sig.channel=1:4}}

  \item{back.channel}{\link{vector}(with default): channels used for background counts subtraction, default \code{back.channel} is set equal to the last 20 channels}

  \item{mr}{\link{double}(with default): measurement errror}

  \item{typ}{\link{character}(with default): type of an OSL decay curve, only CW-OSL decay curve can be analyzed currently}

  \item{nstart}{\link{integer}(with default): maximum number of trials used for initializing parameters of a non-linear growth curve, see function \link{fitGrowth} for details}

  \item{upb}{\link{double}(with default): upper limit on \code{b} value (i.e., the reciprocal of the saturating dose), see function \link{fitGrowth} for details}

  \item{ErrorMethod}{\link{character}(with default): method for assessing the standard error of an equivalent dose, default \code{ErrorMethod="mc"}, see function \link{calED} for details}

  \item{nsim}{\link{integer}(with default): maximum number of simulations if a Monte Carlo method is applied to assessing the standard error of an equivalent dose}

 \item{weight}{\link{logical}(with default): fit the growth curve using a weighted or an un-weighted procedure, see function \link{fitGrowth} for details}

  \item{plot}{\link{logical}(with default): draw a plot or not}

  \item{model}{\link{character} (optional): model used for fitting the growth curve, see function \link{fitGrowth} for available models. The procedure will try to find an optimal model between \code{"line"} and \code{"exp"} automatically if \code{model=NULL}}

  \item{origin}{\link{logical}(optional): force the growth curve to pass the origin or not. The procedure will choose whether to pass the origin or not automatically if \code{origin=NULL}}
}

\details{
Function \emph{analyst} is used for estimating an equivalent dose using dataset obtained from the SAR sequence. It works in a manner similar to software \bold{\emph{Analyst}} written by Duller (2007). See function \link{fitGrowth} and \link{calED} for more details concerning growth curve fitting and equivalent dose calculation.
}

\value{
Return a list containing the following elements:
  \item{Curvedata}{data used for constructing the growth curve}
  \item{Ltx}{sensitivity-corrected natural OSL and its standard error}
  \item{model}{fitting model for the growth curve}
  \item{LMpars}{optimized parameters for the growth curve}
  \item{value}{minimized objective for the growth curve}
  \item{ED}{equivalent dose and its standard error}
  \item{RecyclingRatio}{recycling ratio}
  \item{Recuperation}{recuperation}
}

\references{
Galbraith R, 2002. A note on the variance of a backround-corrected OSL count. Ancient TL, 20(2): 49-51.

Duller, G., 2007. Software Analyst, user mannual. Freely available at \url{http://www.nutech.dtu.dk/english/Products-and-Services/Dosimetry/Radiation-Measurement-Instruments/TL_OSL_reader/Software}.

\bold{Further reading}

Li B, 2007. A note on estimating the error when subtracting background counts from weak OSL signals. Ancient TL, 25(1): 9-14.

Adamiec G, Heer AJ,  Bluszcz A, 2012. Statistics of count numbers from a photomultiplier tube and its implications for error estimation. Radiation Measurements, 47(9): 746-751.
}

\note{
Argument \code{Sigdata} is a column matrix filled with a number of decay curves (i.e., \emph{natural_signal},\cr 
\emph{test_natural_signal}, \emph{Redose[1]_signal}, \emph{test_Redose[1]_signal},....etc) column by column.\cr\cr
Variation in background counts is assumed to follow a Poisson distribution, \bold{Eqn.(3)} of Galbraith (2002) is used to estimate the standard error of a standardised OSL.\cr\cr
Note that the curve fitting procedure will return an error if any standard errors of parameters for the growth curve cannot be estimated by numerical difference-approximation.\cr\cr
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
         sig.channel=1:2, back.channel=3:8)
}
\keyword{growth curve}
\keyword{equivalent dose}
