\name{calSGCED}
\alias{calSGCED}
\alias{calSGCED.default}

\title{
Equivalent dose calculation using SGC or gSGC method
}
\description{
Calculating equivalent doses according to the "standardised growth curve" (SGC) method suggested by Roberts and Duller (2004) or 
the "global standardised growth curve" (gSGC) method suggested by Li et al. (2015).
}

\usage{
calSGCED(SGCdata, pars, model, origin, avgDev, method = "gSGC", 
         errMethod = "sp", nsim = 500, outpdf = NULL)
}

\arguments{
  \item{SGCdata}{\link{data.frame}(\bold{required}): standardised natural-dose signal and regenerative-dose signal data used for SGC equivalent dose calculation, it should contain five columns (i.e., \code{NO}, \code{SAR.Cycle}, \code{Dose}, \code{Signal}, and \code{Signal.Err}), see \link{SARdata} for details. If \code{method="SGC"}, only SAR cycle of "N" need be provided, if \code{method="gSGC"}, both
SAR cycle of "N" and "R" need be supplied}

  \item{pars}{\link{vector}(\bold{required}): optimized parameters of the SGC obtained using function \link{fitGrowth} or \link{lsNORM}}

  \item{model}{\link{character}(\bold{required}): fitting model used for obtaining \code{pars}}

  \item{origin}{\link{logical}(\bold{required}): logical value indicating if established SGC passes the origin}

  \item{avgDev}{\link{numeric}(\bold{required}): average deviation (i.e., average fit error) of the SGC obtained using function \link{fitGrowth} or \link{lsNORM}}

  \item{method}{\link{character}(with default): method used for equivalent dose calculation, i.e.,\cr 
\code{method="SGC"} or \code{method="gSGC"}}

  \item{errMethod}{\link{character}(with default): method used for equivalent dose error assessment. See function \link{calED} for details}

  \item{nsim}{\link{numeric}(with default): desired number of randomly simulated equivalent dose obtained by Monte Carlo simulation}

  \item{outpdf}{\link{character}(optional): if specified, results of SGC equivalent dose calculation will be written to a PDF file named \code{"outpdf"} and saved to the current work directory}
}

\value{
  Return an invisible \link{list} that contains the following elements:
  \item{saturate.NO}{aliquot (grain) number saturated in standardised natural-dose signal}
  \item{failED.NO}{aliquot (grain) number failed in equivalent dose calculation}
  \item{failEDError.NO}{aliquot (grain) number failed in equivalent dose error assessment}
  \item{scale.Ltx}{scaled standardised natural-dose signals and associated standard errors. Note that standardised natural-dose signals will remain un-scaled if \code{method="SGC"}} 
  \item{sgcED}{calculated SGC equivalent doses}
  \item{ConfInt}{68 percent (one sigma) and 95 percent (two sigma) confidence intervals of SGC equivalent doses}
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

\note{
Argument \code{avgDev} will affect estimated standard errors of SGC equivalent doses if the Simple Transformation method is applied. 
However, modification of its values will not affect the standard errors of SGC equivalent doses if the Monte Carlo method is applied. 
}

\examples{
 data(SARdata)
 ### (1) gSGC ED calculation:
 gSGCdata <- SARdata[is.element(SARdata[,2], c("N","R3")),]
 ### gSGCpars were obtained using function "lsNORM".
 gSGCpars <- c(137.440874251, 0.007997863, 2.462035263, -0.321536177)
 avg.error2 <- 0.1111623
 res <- calSGCED(gSGCdata, gSGCpars, method="gSGC", 
                 model="gok", origin=FALSE, avgDev=avg.error2)
 print(res$sgcED)

 ### (2) SGC ED calculation (not run): 
 # SGCdata <- SARdata[is.element(SARdata[,2], "N"),]
 ### SGCpars were obtained using function "fitGrowth".
 # SGCpars <- c(183.474322547,  0.007038048,  4.254287622, -0.337734151)
 # avg.error <- 0.3156259
 # calSGCED(SGCdata, SGCpars, method="SGC", model="gok", 
 #          origin=FALSE, avgDev=avg.error, outpdf="SGCED")

}

\keyword{growth curve}
\keyword{LS-normalisation}
\keyword{equivalent dose}
