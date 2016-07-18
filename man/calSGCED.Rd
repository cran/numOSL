\name{calSGCED}
\alias{calSGCED}
\alias{calSGCED.default}

\title{
Equivalent dose calculation using the standardised growth curve (SGC) method
}
\description{
Calculating equivalent doses according to the standardised growth curve (SGC) method using known parameters.
}

\usage{
calSGCED(Data1, pars, model, origin, 
         ErrorMethod = "sp", outpdf = NULL)
}

\arguments{
  \item{Data1}{\link{data.frame}(\bold{required}): standardised natural dose signal and additional regenerative dose signal data used for SGC equivalent dose calculation, it should contain five columns (i.e., Grain.NO, SAR.Cycle, Dose, Signal, and Signal.Err), see \link{SARdata} for details}

  \item{pars}{\link{vector}(\bold{required}): optimized parameters of the SGC obtained using function \link{lsNORM}}

  \item{model}{\link{character}(\bold{required}): fitting model used for obtaining \code{pars}}

  \item{origin}{\link{logical}(\bold{required}): if established SGC passes the origin or not}

  \item{ErrorMethod}{\link{character}(with default): method used for assessing the standard error of an equivalent dose, only \code{ErrorMethod="sp"} is available}

  \item{outpdf}{\link{character}(optional): if specified, results of SGC equivalent dose calculation will be written to a PDF file named \code{"outpdf"} and saved to the current work directory}
}

\value{
  Return a list that contains the following elements:
  \item{scaleLtx}{scaled standardised natural dose signals and associated standard errors} 
  \item{sgcED}{calculated SGC equivalent doses}
  \item{saturate.NO}{Saturated Grain.NO from which SGC equivalent doses cannot be calculated}
}

\references{
Li B, Roberts RG, Jacobs Z, Li SH, 2015. Potential of establishing a "global standardised growth curve" (gSGC) for optical dating of quartz from sediments. Quaternary Geochronology, 27: 94-104.

Li B, Jacobs Z, Roberts RG, 2016. Investigation of the applicability of standardised growth curves for OSL dating of quartz from Haua Fteah cave, Libya. Quaternary Geochronology, 35: 1-15.

Roberts HM, Duller GAT, 2004. Standardised growth curves for optical dating of sediment using multiple-grain aliquots. Radiation Measurements, 38(2): 241-252.
}

\seealso{
  \link{lsNORM}; \link{fitGrowth}; \link{calED}; \link{calSARED}
}

\examples{
 data(SARdata)
 Data1 <- SARdata[is.element(SARdata[,2], c("N","R3")),]
 pars <- c(137.440874251, 0.007997863, 2.462035263, -0.321536177)
 calSGCED(Data1, pars, model="gok", origin=FALSE)
}

\keyword{growth curve}
\keyword{LS-normalisation procedure}
\keyword{equivalent dose}
