\name{reportSAM}
\alias{reportSAM}
\alias{reportSAM.default}
\title{
Report the results of statistical age models
}
\description{
Return a list that summaries the density distributions of the quantities of interest. 
}
\usage{
reportSAM(obj, burn = 10000, thin = 5, 
          plot = TRUE, outfile = NULL, ...)
}
\arguments{
  \item{obj}{\bold{required}: an object of S3 class \code{"mcAgeModels"}, which is produced by function \link{mcFMM} or \link{mcMAM} }

  \item{burn}{\link{numeric}(with default): the number of iterations to be discarded, i.e., the initial, non-stationary portion of the chain}

  \item{thin}{\link{numeric}(with default): only one sample from each portion that constains every thin values is retained}

  \item{plot}{\link{logical}(with default): plot the results or not}

  \item{outfile}{\link{character}(optional): if specified, simulated chain will be written to a file named \code{"outfile"} in \code{".csv"} format and saved to the current work directory}

  \item{\dots}{do not use}
}

\details{
This function summarizes the output of a Markov chain (the mean values, the standard deviations, and the 95 percent confidence intervals). The numbers of "burn-in" and "thining" are used to refine the simulated chain. The initial portion of the chains may not sufficiently close to stationarity and need to be discarded. The "thining" protocol is used to reduce the posterior correlation of parameters.
}

\value{
Return a list that contains following elements:
\item{pars}{parameters and their standard errors}
\item{quantile}{quantiles of simulated quantities}
\item{maxlik}{the logged maximum likelihood value}
}

\references{
Lunn D, Jackson C, Best N, Thomas A, Spiegelhalter D, 2013. The BUGS Book: A Practical Introduction to Bayesian Analysis. Boca Raton, Florida: CRC Press.

Robert CP, Casella G, 2010. Introducing Monte Carlo Methods with R. Springer.
}

\seealso{
  \link{mcFMM}; \link{mcMAM}
}

\keyword{statistical age models}
\keyword{equivalent dose}
