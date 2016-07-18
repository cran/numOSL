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
  A data frame with 840 observations on the following 5 variables.
  \describe{
    \item{Grain.NO}{Grain number for each aliquot (grain)}
    \item{SAR.Cycle}{natural (N) or regenerative (R) measurement cycles for each Grain.NO}
    \item{Dose}{regenerative doses for each Grain.NO}
    \item{Signal}{standardised OSL signals for each Grain.NO}
    \item{Signal.Err}{measurement error of standardised OSL signals for each Grain.NO}
  }
}

\references{
Li B, Jacobs Z, Roberts RG, 2016. Investigation of the applicability of standardised growth curves for OSL dating of quartz from Haua Fteah cave, Libya. Quaternary Geochronology, 35: 1-15.
}
\examples{
data(SARdata)
head(SARdata)
}
\keyword{OSL dating}
\keyword{growth curve}
\keyword{equivalent dose}
