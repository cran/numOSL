\name{calSARED}
\alias{calSARED}
\alias{calSARED.default}

\title{
SAR equivalent doses calculation and selection
}

\description{
Calculating and selecting a series of equivalent doses in a batch mode according to the single aliquot regenerative-dose (SAR) method (Murray and Wintle, 2000).
}

\usage{
calSARED(obj_analyseBIN, model = "gok", origin = FALSE, 
         errMethod = "sp", nsim = 500, weight = TRUE,  
         trial = TRUE, Tn.above.3BG = TRUE, TnBG.ratio.low = NULL, 
         rseTn.up = NULL, FR.low = NULL, rcy1.range = NULL, 
         rcy2.range = NULL, rcy3.range = NULL, rcp1.up = NULL, 
         rcp2.up = NULL, fom.up = NULL, rcs.up = NULL,
         calED.method = NULL, rseED.up = NULL, use.se = TRUE, 
         outpdf = NULL, outfile = NULL) 
}

\arguments{
  \item{obj_analyseBIN}{\link{list}(\bold{required}): an object of S3 class "analyseBIN" produced by\cr 
function \link{analyseBINdata}}

  \item{model}{\link{character}(with default): model used for growth curve fitting, see \link{fitGrowth} for available models}

  \item{origin}{\link{logical}(with default): logical value indicating if the growth curve should be forced to pass the origin}

  \item{errMethod}{\link{character}(with default): method used for equivalent dose error assessment. See function \link{calED} for details}

  \item{nsim}{\link{integer}(with default): desired number of randomly simulated equivalent dose obtained by Monte Carlo simulation}

  \item{weight}{\link{logical}(with default): logical value indicating if the growth curve should be fitted using a weighted procedure, see function \link{fitGrowth} for details}

   \item{trial}{\link{logical}(with default): logical value indicating if the growth curve should be fitted using other models if the given model fails, see function \link{fitGrowth} for details}

   \item{Tn.above.3BG}{\link{logical}(with default): logical value indicating if aliquot (grain) with Tn below 3 sigma BG should be rejected}

  \item{TnBG.ratio.low}{\link{numeric}(optional): lower limit on ratio of initial Tn signal to BG}

  \item{rseTn.up}{\link{numeric}(optional): upper limit on relative standard error of Tn in percent}

  \item{FR.low}{\link{numeric}(optional): lower limit on fast ratio of Tn}

  \item{rcy1.range}{\link{vector}(optional): a two-element vector indicating the lower and upper limits on recycling ratio 1, 
Example: \code{rcy1.range=c(0.9,1.1)}}

  \item{rcy2.range}{\link{vector}(optional): a two-element vector indicating the lower and upper limits on recycling ratio 2}

  \item{rcy3.range}{\link{vector}(optional): a two-element vector indicating the lower and upper limits on recycling ratio 3}

  \item{rcp1.up}{\link{numeric}(optional): upper limit on recuperation 1 (i.e., ratio of the\cr
sensitivity-corrected zero-dose signal to natural-dose signal) in percent}

  \item{rcp2.up}{\link{numeric}(optional): upper limit on recuperation 2 (i.e., ratio of the\cr
sensitivity-corrected zero-dose signal to maximum regenerative-dose signal)\cr
in percent}

  \item{fom.up}{\link{numeric}(optional): upper limit on figure of merit (FOM) values of growth curves in percent}

  \item{rcs.up}{\link{numeric}(optional): upper limit on reduced chi-square (RCS) values of growth curves}

  \item{calED.method}{\link{character}(optional): method used for equivalent dose calculation, i.e.,\cr
\code{"Interpolation"} or \code{"Extrapolation"}}

  \item{rseED.up}{\link{numeric}(optional): upper limit on the relative standard error of equivalent dose in percent}

  \item{use.se}{\link{logical}(with default): logical value indicating if standard errors of values should be used during application of rejection criteria}

  \item{outpdf}{\link{character}(optional): if specified, results of SAR equivalent dose calculation will be written to a PDF file named \code{"outpdf"} and saved to the current work directory}

  \item{outfile}{\link{character}(optional): if specified, SAR equivalent doses related quantities will be written
to a CSV file named \code{"outfile"} and saved to the current work directory}
}

\value{
Return an invisible \link{list} that contains the following elements:
  \item{LMpars}{a \link{list} containing optimized parameters for successfully fitted growth curves}
  \item{sarED}{calculated (selected) SAR equivalent doses and associated standard errors}
  \item{ConfInt}{68 percent (one sigma) and 95 percent (two sigma) confidence intervals of SAR equivalent doses}
  \item{Tn}{values and standard errors of Tn of calculated (selected) SAR equivalent doses}
  \item{agID}{aliquot (grain) ID of calculated (selected) SAR equivalent doses}
}

\references{
Duller GAT, 2016. Analyst (v4.31.9), User Mannual (\url{http://www.nutech.dtu.dk/english/Products-and-Services/Dosimetry/Radiation-Measurement-Instruments/TL_OSL_reader/Manuals}).

Murray AS, Wintle AG, 2000. Luminescence dating of quartz using improved single-aliquot regenerative-dose protocol. Radiation Measurements, 32(1): 57-73.

Wintle AG, Murray AS, 2006. A review of quartz optically stimulated luminescence characteristics and their relevance in single-aliquot regeneration dating protocols. Radiation Measurements, 41(4): 369-391.
}

\seealso{
  \link{analyseBINdata}; \link{fitGrowth};\cr
  \link{calED}; \link{calSGCED}; \link{pickSARdata}
}

\examples{
  data(BIN)
  obj_pickBIN <- pickBINdata(BIN, Position=c(2,4,6,8,10), Grain=0, 
                             LType="OSL", view=FALSE)
  obj_analyseBIN <- analyseBINdata(obj_pickBIN, nfchn=3, nlchn=20) 
  res_SARED <- calSARED(obj_analyseBIN, model="exp", origin=FALSE)
  # plot(res_SARED$sarED[,1], res_SARED$Tn[,1], xlab="ED (Gy|s)", ylab="Tn")
}

\keyword{growth curve}
\keyword{equivalent dose}
\keyword{data extraction}
