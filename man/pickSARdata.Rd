\name{pickSARdata}
\alias{pickSARdata}
\alias{pickSARdata.default}

\title{
SAR data set selection
}
\description{
Extracting SAR data sets according to certain criteria (recycling ratio, recuperation, figure-of-merit and reduced-chi-square values of fitted growth curves).
}
\usage{
pickSARdata(Data, model = "gok", origin = FALSE, weight = TRUE, 
            rcy.interval = NULL, rcp1.limit = NULL, rcp2.limit = NULL, 
            fom.limit = NULL, rcs.limit = NULL, outfile = NULL)
}

\arguments{
  \item{Data}{\link{data.frame}(\bold{required}): data used for constructing growth curves, it should contain five columns (i.e., Grain.NO, SAR.Cycle, Dose, Signal, and Signal.Err), see \link{SARdata} for details}

  \item{model}{\link{character}(with default): model used for growth curve fitting, see \link{fitGrowth} for available models}

  \item{origin}{\link{logical}(with default): force the fitted growth curve to pass the origin or not}

  \item{weight}{\link{logical}(with default): perform growth curve fitting using a weighted or an un-weighted procedure, see function \link{fitGrowth} for details}

  \item{rcy.interval}{\link{vector}(optional): a two-element vector indicating the lower and upper limits on recycling ratio, for example, \code{rcy.interval=c(0.9,1.1)}}

  \item{rcp1.limit}{\link{double}(optional): upper limit on recuperation1 (i.e., the ratio of standardised zero regenerative dose signal to standardised natural dose signal) in percent}

  \item{rcp2.limit}{\link{double}(optional): upper limit on recuperation2 (i.e., the ratio of standardised zero regenerative dose signal to maximum standardised regenerative dose signal) in percent}

  \item{fom.limit}{\link{double}(optional): upper limit on figure-of-merit (FOM) values of fitted growth curves in percent}

  \item{rcs.limit}{\link{double}(optional): upper limit on reduced-chi-square (RCS) values of fitted growth curves}

  \item{outfile}{\link{character}(optional): if specified, quantities related to selected Grain.NO will be written to a CSV file named \code{"outfile"} and saved to the current work directory}
}

\value{
  Return a list that contains the following elements:
  \item{Data}{Selected SAR data sets}
  \item{tab}{quantities related to seleted Grain.NO, such as recycling ratio, recuperation, FOM, RCS, etc}
  \item{failFit.NO}{Grain.NO failed in growth curve fitting}
  \item{reject.NO}{rejected Grain.NO}
}

\seealso{
   \link{fitGrowth}; \link{lsNORM}; \link{pickSARED}; \link{select}
}

\examples{
 # Do not run.
 # data(SARdata)
 # pickSARdata(SARdata, rcy.interval=c(0.9,1.1), fom.limit=10)
}

\keyword{growth curve}
\keyword{data extraction}
