\name{calSARED}
\alias{calSARED}
\alias{calSARED.default}
\alias{summary.calSARED}
\alias{summary.calSARED.default}

\title{
SAR equivalent doses calculation
}

\description{
Calculating a series of equivalent doses using the single aliquot regenerative-dose (SAR) method (Murray and Wintle, 2000) using a wrapper of function \link{calED}.
}

\usage{
calSARED(Data, model = "gok", origin = FALSE, 
         ErrorMethod = "mc", nsim = 1000, weight = TRUE, 
         trace = TRUE, outpdf = NULL, outfile = NULL)
}

\arguments{
  \item{Data}{\link{data.frame}(\bold{required}): data used for SAR equivalent dose calculation, it should contain five columns (i.e., Grain.NO, SAR.Cycle, Dose, Signal, and Signal.Err), see \link{SARdata} for details}

  \item{model}{\link{character}(with default): model used for growth curve fitting, see \link{fitGrowth} for available models}

  \item{origin}{\link{logical}(with default): force the fitted growth curve to pass the origin or not}

  \item{ErrorMethod}{\link{character}(with default): method used for assessing the standard error of an equivalent dose, see function \link{calED} for details}

  \item{nsim}{\link{integer}(with default): number of Monte Carlo simulations used for equivalent dose error estimation if \code{ErrorMethod="mc"}, see function \link{calED} for details}

  \item{weight}{\link{logical}(with default): perform growth curve fitting using a weighted or an un-weighted procedure, see function \link{fitGrowth} for details}

  \item{trace}{\link{logical}(with default): logical value indicating if a trace of the calculating progress should be printed}

  \item{outpdf}{\link{character}(optional): if specified, results of SAR equivalent dose calculation will be drawn to a PDF file named \code{"outpdf"} and saved to the current work directory}

  \item{outfile}{\link{character}(optional): if specified, calculated equivalent doses and relevant quantities will be written to a CSV file named \code{"outfile"} and saved to the current work directory}
}

\value{
Return an invisible list of S3 class "calSARED" that contains the following elements:
  \item{LMpars}{parameters of fitted growth curves for each Grain.NO}
  \item{N}{number of analyzed Grain.NO}
  \item{failFit.NO}{Grain.NO failed in growth curve fitting}
  \item{saturate.NO}{Grain.NO with saturated standardised natural dose signal}
  \item{failED.NO}{Grain.NO failed in equivalent dose (error) calculation}
  \item{extrapolate.NO}{Grain.NO from which equivalent dose is calculated by extrapolation}
  \item{largeRcy.NO}{Grain.NO with recycling ratio falls outside the range 0.9-1.1}
  \item{largeRcp5.NO}{Grain.NO with recuperation (i.e., ratio of the standardised zero regenerative dose signal to natural dose signal) larger than 5 percent}
  \item{largeRcp10.NO}{Grain.NO with recuperation larger than 10 percent}
  \item{tab}{calculated SAR equivalent doses and relevant quantities, such as recycling ratio, recuperation, etc}
}

\references{
Murray AS, Wintle AG, 2000. Luminescence dating of quartz using improved single-aliquot regenerative-dose protocol. Radiation Measurements, 32(1): 57-73.
}

\seealso{
  \link{fitGrowth}; \link{calED}; \link{calSGCED}
}

\examples{
  data(SARdata)
  # Use only the first two aliquots of SARdata.
  Data <- SARdata[1:16,]
  res <- calSARED(Data, model="gok", trace=FALSE)
  summary(res)
  print(res$LMpars)
}

\keyword{growth curve}
\keyword{equivalent dose}
