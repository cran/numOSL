\name{sensSAM}
\alias{sensSAM}
\alias{sensSAM.default}
\title{
Investigate of the sensitivity of a statistical age model to the additional uncertainty (sigmab)
}
\description{
Estimate of the parameters of a statistical age model using a number of sigmab values.
}
\usage{
sensSAM(EDdata, model, sigmaVEC = NULL, iflog = TRUE, 
        maxcomp = 8, plot = TRUE, outfile = NULL)
}
\arguments{
  \item{EDdata}{\link{matrix}(\bold{required}): a two-column matrix (i.e., equivalent dose values and\cr
associated standard errors)}

  \item{model}{\link{character}(with default): the fitting model, one of \code{"com"}, \code{"cam"}, \code{"mam3"}, \code{"mam4"},
\code{"mxam3"}, \code{"mxam4"}, \code{"fmm0"}, \code{"fmm1"}, \code{"fmm2"}, ..., \code{"fmm9"}}

  \item{sigmaVEC}{\link{vector}(with default): a series of sigmab values that will be used as inputs for the model. For example, \code{sigmaVEC=seq(from=0,to=0.3,by=0.01)}}

  \item{iflog}{\link{logical}(with default): transform equivalent dose values to log-scale or not}

  \item{maxcomp}{\link{integer}(with default): the maximum number of components in the FMM}

  \item{plot}{\link{logical}(with default): logical value indicating if the results should be plotted}

  \item{outfile}{\link{character}(optional): if specified, the results will be written to a CSV file named \code{"outfile"} and saved to the current work directory}
}

\value{
Return an invisible \link{list} that contains the following elements:
  \item{pars}{a list that contains the optimized parameters for each sigmab value}
  \item{mat}{a matrix that contains the optimized parameters, the maximum logged likelihood value, and the corresponding Bayesian Information Criterion (BIC) value}
}

\references{
Peng J, Li B, Jacobs Z, 2020. Modelling heterogeneously bleached single-grain equivalent dose distributions: Implications for the reliability of burial dose determination. Quaternary Geochronology, 60: 101108.

Peng J, Li B, Jacobs Z, Gliganic LA, 2023. Optical dating of sediments affected by post-depositional mixing: Modelling, synthesizing and implications. Catena, 232: 107383.
}

\seealso{
 \link{RadialPlotter}; \link{EDdata}; \link{optimSAM}
}

\examples{

  # Not run.
  # data(EDdata)
  # sensSAM(EDdata$al3, model="mam4", iflog=TRUE)

}
\keyword{statistical age models}
\keyword{equivalent dose}
\keyword{Central Age Model}
\keyword{Finite Mixture Age Model}
\keyword{Minimum Age Model}
\keyword{Maximum Age Model}
