\name{analyst}
\alias{analyst}
\alias{analyst.default}

\title{
Analysing equivalent dose values.  
}

\description{
Calculating an equivalent dose value and estimating its standard error using dataset from a SAR sequence. 
}

\usage{
analyst(Sigdata, Redose, sig.channel = NULL, back.channel = NULL, 
        mr = 0.01, disb = c("p","op"), typ = "cw", nstart = 100,
        upb = 0.5, ErrorMethod = c("mc","sp"), nsim = 1000,
        weight = TRUE, plot = TRUE, model = NULL, origin = NULL)
}

\arguments{
  \item{Sigdata}{\link{data.frame}(\bold{required}): a number of decay curves, the first column should be the stimulation time values}

  \item{Redose}{\link{vector}(\bold{required}): the regenerative dose values. 
Example, \code{Redose=c(1,2,3,4,0,1)}}

  \item{sig.channel}{\link{vector}(with default): channels used for signal integration, 
default \code{sig.channel=1:4}}

  \item{back.channel}{\link{vector}(with default): channels used for background counts subtraction, 
default \code{back.channel} is set equal to the last 20 channels}

  \item{mr}{\link{numeric}(with default): the measurement errror}

  \item{disb}{\link{character}(with default): the assumed distribution of background counts, \code{"p"} means Possion, and \code{"op"} means over-Possion}

  \item{typ}{\link{character}(with default): the type of an OSL decay curve, only type of \code{"CW"} is available currently}

  \item{nstart}{\link{numeric}(with default): the maximum number of attempts for initializing the non-linear growth curve, 
see function \link{fitGrowth} for details
}

  \item{upb}{\link{numeric}(with default): the upper limit on \code{b} value (e.g., the reciprocal of the saturating dose), see function \link{fitGrowth} for details
}

  \item{ErrorMethod}{\link{character}(with default): a method for assessing the standard error of an equivalent dose, default \code{ErrorMethod="mc"}, see function \link{calED} for details}

  \item{nsim}{\link{numeric}(with default): the allowed maximum number of simulations when applying a Monte Carlo method to assess the standard error of an equivalent dose}

 \item{weight}{\link{logical}(with default): fit the growth curve using a weighted or an un-weighted procedure, see function \link{fitGrowth} for details}

  \item{plot}{\link{logical}(with default): draw a plot or not}

  \item{model}{\link{character}(optional): a model used for fitting the growth curve, see function \link{fitGrowth} for available models. It will find a model automatically if \code{model=NULL}}

  \item{origin}{\link{logical}(optional): force the fitting to pass the origin or not. It will choose whether pass the origin or not automatically if \code{origin=NULL}}
}

\details{
Function \emph{analyst} is used for estimating an equivalent dose and assessing its standard error using a SAR dataset. It analyses an ED value in a manner similar to software \bold{\emph{Analyst}} by Duller (2007). See function \link{fitGrowth} for more details concerning growth curve fitting.
}

\value{
Return a list containing following elements:
  \item{Curvedata}{data used to construct the growth curve}
  \item{Ltx}{the natural standardised OSL and its standard error}
  \item{model}{the fitting  model}
  \item{LMpars}{the parameters of the growth curve}
  \item{value}{the minimized objective of the growth curve}
  \item{ED}{the equivalent dose and its standard error}
  \item{RecyclingRatio}{the recycling ratio}
  \item{Recuperation}{the recuperation}
}

\references{
Galbraith R, 2002. A note on the variance of a backround-corrected OSL count. Ancient TL, 20: 49-51.

Duller, G., 2007. Analyst. Use mannual. Available at \url{http://www.nutech.dtu.dk/english/Products-and-Services/Dosimetry/Radiation-Measurement-Instruments/TL_OSL_reader/Software}.

\bold{Further reading}

Li B, 2007. A note on estimating the error when subtracting background counts from weak OSL signals. Ancient TL, 25: 9-14.

Adamiec G, Heer AJ,  Bluszcz A, 2012. Statistics of count numbers from a photomultiplier tube and its implications for error estimation. Radiation Measurements, 47: 746-751.
}

\note{
Argument \code{Sigdata} should be a data structure that seems like: \emph{stimulation_time}, \emph{natural_signal},
\emph{test_natural_signal}, \emph{Redose[1]_signal}, \emph{test_Redose[1]_signal}, \emph{...}, ect.\cr\cr
If the variation in background counts is assumed to follow a Poisson distribution, then \bold{Eqn.(3)} of Galbraith (2002) will be used to estimate the standard error of a standardised OSL. If the variation in background counts is assumed to be over-dispersed, then \bold{Eqn.(6)}  of Galbraith (2002) will be used.\cr\cr
Note that the number of background channels \code{(m)} need be a multiple of the number of integration channels \code{(n)}, and that it will return an error if any standard errors of parameters of the growth curve cannot be estimated by numerical approximation.\cr\cr
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
 analyst(Signaldata[[1]], Redose=c(80,160,240,320,0,80)*0.13,
         sig.channel=1:2, back.channel=3:8)
}
\keyword{growth curve}
\keyword{equivalent dose}
