\name{dbED}
\alias{dbED}
\alias{dbED.default}

\title{
Summarizing the distribution of equivalent dose values
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
  \item{EDdata}{\link{matrix}(\bold{required}): two columns, i.e., equivalent dose values and their standard errors}
 
  \item{plot}{\link{logical}(with default): draw a plot or not}
  
  \item{typ}{\link{character}(with default): type of plot, \code{"pdf"} means a probability density plot, \code{"hist"} means a histogram, and \code{"scatter"} means a scatterplot, default \code{typ="pdf"}}

  \item{from}{\link{numeric}(optional): the desired lower limit on x-axis}

  \item{to}{\link{numeric}(optional): the desired upper limit on x-axis}

  \item{step}{\link{numeric}(with default): a step-size used for constructing the probability density plot. Smaller step value gives smoother density curve}

  \item{nbin}{\link{numeric}(with default): the desired number of intervals of a histogram}

  \item{pcolor}{\link{character}(with default): the color of a data point, input \link{colors}() to see more available colors}

  \item{psize}{\link{numeric}(with default): the size of a data point}
}

\value{
Return a list that contains following elements:
  \item{weight.ED}{the weigthed mean of equivalent dose values and its standard error}
  \item{skewness}{the weighted skewness of equivalent dose values and its standard error}
  \item{kurtosis}{the kurtosis of equivalent dose values and its standard error}
  \item{quantile.ED}{the quantiles of ED values}
}

\references{
Bailey RM, Arnold LJ, 2006. Statistical modelling of single grain quartz De distributions and an assessment of procedures for estimating burial dose. Quaternary Science Reviews, 25: 2475-2502.

Galbraith RF, 2010. On plotting OSL equivalent doses. Ancient TL, 28: 1-10.

Galbraith RF, Roberts RG, 2012. Statistical aspects of equivalent dose and error calculation and display in OSL dating: An overview and some recommendations. Quaternary Geochronology, 11: 1-27.
}

\seealso{
\link{RadialPlotter}
}

\examples{
 data(EDdata)
 dbED(EDdata$gl11,typ="pdf")
}

\keyword{equivalent dose}
