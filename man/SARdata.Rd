\name{SARdata}
\alias{SARdata}
\docType{data}
\title{
Data sets used for SAR equivalent dose calculation
}
\description{
SAR data sets for individual aliquots from the "later" group of sample HF11 from Haua Fteah cave, Libya (Li et al., 2016).
}
\usage{data(SARdata)}
\format{
  A \link{data.frame} with 840 observations containing the following 5 variables:
  \describe{
    \item{NO}{aliquot (grain) number}
    \item{SAR.Cycle}{SAR cycles}
    \item{Dose}{regenerative doses}
    \item{Signal}{OSL signals}
    \item{Signal.Err}{standard error of OSL signals}
  }
}

\references{
Li B, Jacobs Z, Roberts RG, 2016. Investigation of the applicability of standardised growth curves for OSL dating of quartz from Haua Fteah cave, Libya. Quaternary Geochronology, 35: 1-15.
}

\seealso{
  \link{fitGrowth}; \link{lsNORM}; \link{calSGCED};\cr
  \link{as_analyseBIN}
}

\examples{
 # Not run.
 # data(SARdata)
 # head(SARdata)
}
\keyword{OSL dating}
\keyword{growth curve}
\keyword{equivalent dose}
