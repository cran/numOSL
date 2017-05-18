\name{calSGCED}
\alias{calSGCED}
\alias{calSGCED.default}

\title{
SGC Equivalent dose calculation and selection
}
\description{
Calculating and selecting equivalent doses according to the "standardised growth curve" (SGC) method suggested 
by Roberts and Duller (2004) or the "global standardised growth curve" (gSGC) method suggested by Li et al. (2015).
}

\usage{
calSGCED(obj_analyseBIN, SGCpars, model, origin, avgDev, 
         method = "SGC", SAR.Cycle = "N", errMethod = "sp", 
         Tn.above.3BG = TRUE, TnBG.ratio.low = NULL, rseTn.up = NULL,
         FR.low = NULL, rseED.up = NULL, use.se = TRUE, outpdf = NULL,
         outfile = NULL)
}

\arguments{
   \item{obj_analyseBIN}{\link{list}(\bold{required}): an object of S3 class "analyseBIN" produced by\cr 
function \link{analyseBINdata} or \link{as_analyseBIN}}

  \item{SGCpars}{\link{vector}(\bold{required}): optimized parameters of the SGC obtained using function \link{fitGrowth} or \link{lsNORM}}

  \item{model}{\link{character}(\bold{required}): fitting model used for obtaining \code{SGCpars}}

  \item{origin}{\link{logical}(\bold{required}): logical value indicating if established SGC passes the origin}

  \item{avgDev}{\link{numeric}(\bold{required}): average deviation (i.e., average fit error) of the SGC obtained using function \link{fitGrowth} or \link{lsNORM}. This quantity stands for the uncertainty of established SGC when assessing the equivalent dose error using the simple transformaion method}

  \item{method}{\link{character}(with default): method used for equivalent dose calculation, i.e.,\cr 
\code{method="SGC"} (for the original SGC method) or \code{method="gSGC"} (for the improved SGC method)}

  \item{SAR.Cycle}{\link{character}(with default): SAR cycles used for SGC equivalent dose calculation.\cr 
Example: \code{SAR.Cycle=c("N","R3")}}

  \item{errMethod}{\link{character}(with default): method used for equivalent dose error assessment}

   \item{Tn.above.3BG}{\link{logical}(with default): logical value indicating if aliquot (grain) with Tn below 3 sigma BG should be rejected}

  \item{TnBG.ratio.low}{\link{numeric}(optional): lower limit on ratio of initial Tn signal to BG}

  \item{rseTn.up}{\link{numeric}(optional): upper limit on relative standard error of Tn in percent}

  \item{FR.low}{\link{numeric}(optional): lower limit on fast ratio of Tn}

  \item{rseED.up}{\link{numeric}(optional): upper limit on the relative standard error of equivalent dose in percent}

  \item{use.se}{\link{logical}(with default): logical value indicating if standard errors of values should be used during application of rejection criteria}

  \item{outpdf}{\link{character}(optional): if specified, results of SGC equivalent dose calculation will be written to a PDF file named \code{"outpdf"} and saved to the current work directory}

  \item{outfile}{\link{character}(optional): if specified, SGC equivalent doses related quantities will be written to a CSV file named \code{"outfile"} and saved to the current work directory}
}

\value{
  Return an invisible \link{list} that contains the following elements:
  \item{scale.Ltx}{scaled standardised natural-dose signals and associated standard errors used for SGC equivalent dose calculation. Note that standardised natural-dose signals will remain un-scaled if \code{method="SGC"}} 
  \item{sgcED}{calculated SGC equivalent doses}
  \item{ConfInt}{68 percent (one sigma) and 95 percent (two sigma) confidence intervals of SGC equivalent doses}
  \item{agID}{aliquot (grain) ID of calculated (selected) SGC equivalent doses}
}

\references{
Li B, Roberts RG, Jacobs Z, Li SH, 2015. Potential of establishing a "global standardised growth curve" (gSGC) for optical dating of quartz from sediments. Quaternary Geochronology, 27: 94-104.

Li B, Jacobs Z, Roberts RG, 2016. Investigation of the applicability of standardised growth curves for OSL dating of quartz from Haua Fteah cave, Libya. Quaternary Geochronology, 35: 1-15.

Roberts HM, Duller GAT, 2004. Standardised growth curves for optical dating of sediment using multiple-grain aliquots. Radiation Measurements, 38(2): 241-252.
}

\seealso{
  \link{fitGrowth}; \link{lsNORM}; \link{SARdata};\cr 
  \link{scaleSGCN}; \link{calED}; \link{calSARED}
}

\examples{
 data(SARdata)
 ### (1) gSGC ED calculation:
 ### gSGCpars were obtained using function "lsNORM".
 gSGCpars <- c(137.440874251, 0.007997863, 2.462035263, -0.321536177)
 avg.error2 <- 0.1111623
 res <- calSGCED(as_analyseBIN(SARdata), gSGCpars, method="gSGC", 
                 model="gok", origin=FALSE, avgDev=avg.error2,
                 SAR.Cycle=c("N","R3"))
 print(res$sgcED)

 ### (2) SGC ED calculation (not run): 
 ### SGCpars were obtained using function "fitGrowth".
 # SGCpars <- c(183.474322547,  0.007038048,  4.254287622, -0.337734151)
 # avg.error <- 0.3156259
 # calSGCED(as_analyseBIN(SARdata), SGCpars, method="SGC", model="gok", 
 #          origin=FALSE, avgDev=avg.error, SAR.Cycle="N", outpdf="SGCED")

 ### (3) gSGC ED calculation and signal-related 
 ###     rejection criteria application (not run):
 # data(BIN)
 # res_pickBIN <-pickBINdata(BIN, LType="OSL")
 # res_analyseBIN <- analyseBINdata(res_pickBIN, nfchn=4, nlchn=30)
 # res_lsNORM <- lsNORM(res_analyseBIN$SARdata, model="gok", origin=FALSE)

 # calSGCED(res_analyseBIN, SGCpars=res_lsNORM$LMpars2[,1], 
 #         model="gok", origin=FALSE, avgDev=res_lsNORM$avg.error2,
 #         method="gSGC", SAR.Cycle=c("N","R3"), Tn.above.3BG=TRUE, 
 #         TnBG.ratio.low=4, rseTn.up=10, outpdf="foo", outfile="foo")


}

\keyword{growth curve}
\keyword{LS-normalisation}
\keyword{equivalent dose}
