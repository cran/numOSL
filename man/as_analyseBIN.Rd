\name{as_analyseBIN}
\alias{as_analyseBIN}
\alias{as_analyseBIN.default}

\title{
Transfom data sets into S3 class object "analyseBIN"
}

\description{
Transfom data sets into S3 class object "analyseBIN"
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

  \item{criteria}{values used as rejection criteria, values are set as: \code{Tn3BG=1}, \code{TnBG.ratio=99999},\cr
\code{seTnBG.ratio=0}, \code{rseTn=0}, \code{FR=99999}, \code{seFR=0}}

  \item{Tn}{values of Tn and associated standard errors, values are set as: \code{Tn=99999}, \code{seTn=0}}

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
S3 class object "analyseBIN" produced using function \link{analyseBINdata} contains five elements (i.e., \code{SARdata}, \code{criteria}, \code{Tn}, \code{TxTn}, \code{agID}).
However, function \link{as_analyseBIN} uses only \code{SARdata} as input and other elements are lost. Assume that user-supplied \code{SARdata} are selected beforehand using rejection criteria such as 
\code{Tn3BG}, \code{TnBG.ratio}, \code{rseTn}, \code{FR}, values in \code{criteria} are reset so as to ensure that result produced by function \link{as_analyseBIN} can pass those rejection criteria 
if it is used as input of function \link{calSARED} or \link{pickSARdata}.
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
