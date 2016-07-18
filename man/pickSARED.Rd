\name{pickSARED}
\alias{pickSARED}
\alias{pickSARED.default}

\title{
Single aliquot regenerative-dose (SAR) equivalent dose record selection
}

\description{
Extracting SAR equivalent dose records according to certain criteria (recycling ratio, recuperation, figure-of-merit and reduced-chi-square values of fitted growth curves).
}

\usage{
pickSARED(obj, rsdED.limit = NULL, rcy.interval = NULL, 
          rcp1.limit = NULL, rcp2.limit = NULL, method = NULL, 
          fom.limit = NULL, rcs.limit = NULL, outfile = NULL)
}

\arguments{
  \item{obj}{\link{list}(\bold{required}): an object of S3 class \code{"calSARED"}, which is produced by function \link{calSARED}}

  \item{rsdED.limit}{\link{double}(optional): upper limit on the relative standard error of equivalent dose in percent}

  \item{rcy.interval}{\link{vector}(optional): a two-element vector indicating the lower and upper limits on recycling ratio, for example, \code{rcy.interval=c(0.9,1.1)}}

  \item{rcp1.limit}{\link{double}(optional): upper limit on recuperation1 (i.e., the ratio of standardised zero regenerative dose signal to standardised natural dose signal) in percent}

  \item{rcp2.limit}{\link{double}(optional): upper limit on recuperation2 (i.e., the ratio of standardised zero regenerative dose signal to maximum standardised regenerative dose signal) in percent}

  \item{method}{\link{character}(optional): method used for SAR equivalent dose calculation (interpolation or extrapolation)}

  \item{fom.limit}{\link{double}(optional): upper limit on figure-of-merit (FOM) values of fitted growth curves in percent}

  \item{rcs.limit}{\link{double}(optional): upper limit on reduced-chi-square (RCS) values of fitted growth curves}

  \item{outfile}{\link{character}(optional): if specified, selected SAR equivalent doses and relevant quantities will be written to a CSV file named \code{"outfile"} and saved to the current work directory}
}

\value{
 Return a list that contains the following elements:
 \item{sarED}{selected SAR equivalent doses}
 \item{reject.NO}{rejected Grain.NO}
}

\seealso{
  \link{calSARED}; \link{pickSARdata}; \link{select}
}

\examples{
  # Do not run.
  # data(SARdata)
  # res <- calSARED(SARdata, model="gok", trace=TRUE)
  # pickSARED(res, rcy.interval=c(0.9, 1.1))
}

\keyword{equivalent dose}
\keyword{data extraction}
