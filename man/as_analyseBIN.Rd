\name{as_analyseBIN}
\alias{as_analyseBIN}
\alias{as_analyseBIN.default}

\title{
Transfom SAR data sets into S3 class object "analyseBIN"
}

\description{
Transfom SAR data sets into S3 class object "analyseBIN".
}

\usage{
as_analyseBIN(SARdata)
}

\arguments{
  \item{SARdata}{\link{matrix}(\bold{required}): SAR data set, it should contain five columns\cr
(i.e., \code{NO}, \code{SAR.Cycle}, \code{Dose}, \code{Signal}, and \code{Signal.Err}), see \link{SARdata} for details}
}


\value{
Return an invisible \link{list} of S3 class object "analyseBIN" containing the following elements:
  \item{SARdata}{a \link{data.frame} containing SAR data sets}

  \item{criteria}{values used as rejection criteria, here it is set equal to \code{NULL}}

  \item{Tn}{values of Tn and associated standard errors, here it is set equal to \code{NULL}}

  \item{LnTn.curve}{decay curves of Ln and Tn for different aliquots (grains), here it is set equal to \code{NULL}}

  \item{TxTn}{ratios of Tx to Tn for various SAR cycles, here it is set equal to \code{NULL}}

  \item{agID}{aliquot or grain ID (i.e., \code{NO}, \code{Position}, and \code{Grain}), here both \code{Position} and \code{Grain} are set equal to 0}

\code{SARdata} is a \link{data.frame} containing the following elements:
\tabular{lll}{
\bold{Element}      \tab  \bold{Description}                             \cr
\code{NO}           \tab  aliquot (grain) number                         \cr
\code{SAR.Cycle}    \tab  SAR cycle (N, R1, R2, R3, ...)                 \cr 
\code{Dose}         \tab  regenerative dose                              \cr
\code{Signal}       \tab  OSL signal                                     \cr
\code{Signal.Err}   \tab  standard error of OSL signal                                  
}
}

\seealso{
  \link{analyseBINdata}; \link{SARdata};\cr
  \link{calSARED}; \link{pickSARdata}
}

\note{
Function \link{as_analyseBIN} transforms SAR data sets (see \link{SARdata}) into S3 class object "analyseBIN". 
The returned elements such as \code{criteria}, \code{Tn}, \code{LnTn.curve}, and \code{TxTn} are set equal to \code{NULL}.
}

\examples{
  ### Example 1:
  data(SARdata)
  obj_analyseBIN <- as_analyseBIN(SARdata[1:8,,drop=FALSE])
  res_calSARED <- calSARED(obj_analyseBIN)
  res_calSARED$sarED

  ### Example 2 (not run):
  # obj_analyseBIN <- as_analyseBIN(SARdata)
  # res_calSARED <- calSARED(obj_analyseBIN, rcy1.range=c(1,1), outpdf="SARED")

  ### Example 3 (not run):
  # obj_analyseBIN <- as_analyseBIN(SARdata)
  # res_pickSARdata <- pickSARdata(obj_analyseBIN, fom.up=6, outpdf="SARdata")
  # res_pickSARdata$SARdata
}

\keyword{growth curve}
