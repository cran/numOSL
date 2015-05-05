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
dbED(EDdata, plot = TRUE, typ = c("pdf","hist","scatter"), 
     from = NULL, to = NULL, step = 0.1, nbin = 15, 
     pcolor = "grey", psize = 1.5)
}

\arguments{
  \item{EDdata}{\link{matrix}(\bold{required}): a two-column matrix, equivalent dose values and associated errors are stored in the first and second column, respectively}
 
  \item{plot}{\link{logical}(with default): draw a plot or not}
  
  \item{typ}{\link{character}(with default): type of plot, \code{"pdf"} means a probability density plot, \code{"hist"} means a histogram, and \code{"scatter"} means a scatterplot, default \code{typ="pdf"}}

  \item{from}{\link{double}(optional): desired lower limit on x-axis}

  \item{to}{\link{double}(optional): desired upper limit on x-axis}

  \item{step}{\link{double}(with default): a step-size used for constructing the probability density plot. Smaller step value gives smoother density curve}

  \item{nbin}{\link{integer}(with default): desired number of intervals of a histogram}

  \item{pcolor}{\link{character}(with default): color of a data point, input \link{colors}() to see available colors}

  \item{psize}{\link{double}(with default): size of a data point}
}

\value{
Return a list that contains the following elements:
  \item{weight.ED}{weigthed mean of equivalent dose values and associated standard error}
  \item{skewness}{weighted skewness of equivalent dose values and associated standard error}
  \item{kurtosis}{kurtosis of equivalent dose values and associated standard error}
  \item{quantile.ED}{quantiles of equivalent dose values}
}

\references{
Bailey RM, Arnold LJ, 2006. Statistical modelling of single grain quartz De distributions and an assessment of procedures for estimating burial dose. Quaternary Science Reviews, 25(19): 2475-2502.

Galbraith RF, 2010. On plotting OSL equivalent doses. Ancient TL, 28(1): 1-10.

Galbraith RF, Roberts RG, 2012. Statistical aspects of equivalent dose and error calculation and display in OSL dating: an overview and some recommendations. Quaternary Geochronology, 11: 1-27.
}

\seealso{
  \link{psRadialPlot}; \link{RadialPlotter}
}

\examples{
 data(EDdata)
 dbED(EDdata$gl11,typ="pdf")
}

\keyword{equivalent dose}
