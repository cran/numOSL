\name{pickSARdata}
\alias{pickSARdata}
\alias{pickSARdata.default}

\title{
SAR data set selection
}
\description{
Selecting SAR data sets (growth curves) in a batch model according to specified rejection criteria.
}
\usage{
pickSARdata(obj_analyseBIN, model = "gok", origin = FALSE, 
            weight = TRUE, trial = TRUE, nofit.rgd = NULL, 
            Tn.above.3BG = TRUE, TnBG.ratio.low = NULL, 
            rseTn.up = NULL, FR.low = NULL, rcy1.range = NULL, 
            rcy2.range = NULL, rcy3.range = NULL, 
            rcp1.up = NULL, rcp2.up = NULL, fom.up = NULL, 
            rcs.up = NULL, use.se = TRUE, norm.dose = NULL, 
            outpdf = NULL, outfile = NULL)
}

\arguments{
  \item{obj_analyseBIN}{\link{list}(\bold{required}): an object of S3 class "analyseBIN" produced by\cr 
function \link{analyseBINdata} or \link{as_analyseBIN}}

  \item{model}{\link{character}(with default): model used for growth curve fitting, see \link{fitGrowth} for available models}

  \item{origin}{\link{logical}(with default): logical value indicating if the growth curve should be forced to pass the origin}

  \item{weight}{\link{logical}(with default): logical value indicating if the growth curve should be fitted using a weighted procedure, 
see function \link{fitGrowth} for details}

  \item{trial}{\link{logical}(with default): logical value indicating if the growth curve should be fitted using other models if the given model fails, see function \link{fitGrowth} for details}

   \item{nofit.rgd}{\link{integer}(optional): regenerative doses that will not be used during the fitting. 
For example, if \code{nofit.rgd=2} then the second regenerative dose will not be used during growth curve fitting}

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

  \item{use.se}{\link{logical}(with default): logical value indicating if standard errors of values should be used during application of rejection criteria}

  \item{norm.dose}{\link{numeric}(optional): dose value used for SAR data set re-normalisation, 
for example, if \code{norm.dose=100}, then sensitivity-corrected signal for \code{Redose=100} obtained 
through growth curve fitting will be used to re-normalise a SAR data set}

  \item{outpdf}{\link{character}(optional): if specified, results of growth curve fitting will be written to 
a PDF file named \code{"outpdf"} and saved to the current work directory}

  \item{outfile}{\link{character}(optional): if specified, SAR data related quantities will be written
to a CSV file named \code{"outfile"} and saved to the current work directory}
}

\value{
  Return an invisible \link{list} that contains the following elements:
  \item{LMpars}{a \link{list} containing optimized parameters of growth curves of selected SAR data sets}
  \item{SARdata}{a \link{data.frame} containing selected SAR data sets}
  \item{norm.SARdata}{a \link{data.frame} containing re-normalised SAR data sets,\cr
it returns \code{NULL} if \code{norm.dose=NULL}}
  \item{agID}{aliquot or grain ID (i.e., \code{c("NO","Position","Grain")}) of selected SAR data}
  \item{summary.info}{a summary of the SAR data selection}
}

\seealso{
   \link{analyseBINdata}; \link{fitGrowth}; \link{lsNORM}; \link{calSARED}
}

\examples{
 # Not run.
 # data(BIN)
 # obj_pickBIN <- pickBINdata(BIN, Position=1:48, Grain=0, 
 #                            LType="OSL", view=FALSE)
 # obj_analyseBIN <- analyseBINdata(obj_pickBIN, nfchn=3, nlchn=20) 
 # pickSARdata(obj_analyseBIN, model="gok", fom.up=3, outpdf="SARdata")
}

\keyword{growth curve}
\keyword{data extraction}
\keyword{batch analysis}
