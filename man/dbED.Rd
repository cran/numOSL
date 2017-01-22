\name{dbED}
\alias{dbED}
\alias{dbED.default}

\title{
De distribution summarization 
}

\description{
Calculating statistical parameters (skewness, kurtosis, quantiles) for a number of equivalent dose values.
}

\usage{
dbED(EDdata, plot = TRUE, typ = "pdf", from = NULL, 
     to = NULL, step = NULL, nbin = 15, pcolor = "grey", 
     psize = 1.5, outfile = NULL)
}

\arguments{
  \item{EDdata}{\link{matrix}(\bold{required}): a two-column matrix (i.e., equivalent dose values and\cr
associated standard errors)}
 
  \item{plot}{\link{logical}(with default): draw a plot or not}
  
  \item{typ}{\link{character}(with default): type of plot, \code{typ="pdf"} means a probability density plot and \code{typ="hist"} means a histogram plot}

  \item{from}{\link{numeric}(optional): desired lower limit on x-axis}

  \item{to}{\link{numeric}(optional): desired upper limit on x-axis}

  \item{step}{\link{numeric}(optional): a step-size used for constructing the probability density plot (if \code{typ="pdf"}). Smaller step value gives smoother density curve}

  \item{nbin}{\link{integer}(with default): desired number of intervals for the histogram (if \code{typ="hist"})}

  \item{pcolor}{\link{character}(with default): color of data points, input \link{colors}() to see available colors}

  \item{psize}{\link{numeric}(with default): size of data points}
  
  \item{outfile}{\link{character}(optional): if specified, calculated probability densities (if \code{typ="pdf"}) will be written to a CSV file named \code{"outfile"} and saved to the current work directory}

}

\value{
Return a \link{list} that contains the following elements:
  \item{weight.ED}{weigthed mean of equivalent dose values and associated standard error}
  \item{skewness}{weighted skewness of equivalent dose values and associated standard error}
  \item{kurtosis}{kurtosis of equivalent dose values and associated standard error}
  \item{quantile.ED}{quantiles of equivalent dose values}
}

\references{
Galbraith RF, 2010. On plotting OSL equivalent doses. Ancient TL, 28(1): 1-10.

Galbraith RF, Roberts RG, 2012. Statistical aspects of equivalent dose and error calculation and display in OSL dating: an overview and some recommendations. Quaternary Geochronology, 11: 1-27.
}

\seealso{
  \link{psRadialPlot}; \link{RadialPlotter}; \link{EDdata}
}

\examples{
 data(EDdata)
 dbED(EDdata$gl11,typ="pdf")
}

\keyword{equivalent dose}
