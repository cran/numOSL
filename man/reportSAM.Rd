\name{reportSAM}
\alias{reportSAM}
\alias{reportSAM.default}
\title{
Reporting MCMC outputs for statistical age models
}
\description{
Summarizing distributions of parameters simulated from statistical age models. 
}
\usage{
reportSAM(obj, burn = 10000, thin = 5, 
          plot = TRUE, outfile = NULL, ...)
}
\arguments{
  \item{obj}{\bold{required}: an object of S3 class \code{"mcAgeModels"}, which is produced by function \link{mcFMM} or \link{mcMAM} }

  \item{burn}{\link{integer}(with default): number of iterations (i.e., the initial, non-stationary portion of the chain) to be discarded}

  \item{thin}{\link{integer}(with default): take every \code{thin-th} iteration}

  \item{plot}{\link{logical}(with default): plot the MCMC output or not}

  \item{outfile}{\link{character}(optional): if specified, simulated parameters will be written to a file named \code{"outfile"} in \code{CSV} format and saved to the current work directory}

  \item{\dots}{do not use}
}

\details{
Function \emph{reportSAM} summarizes the output of a Markov chain (the mean values, the standard deviations, and the 95 percent confidence intervals). The initial \emph{i} (\code{burn=i}) samples may have been affected by the inital state and has to be discarded ("burn-in"). Autocorrelation of simulated samples can be reduced by taking every \emph{j-th} (\code{thin=j}) iteration ("thining").
}

\value{
Return a list that contains the following elements:
\item{pars}{parameters and relevant statistics}
\item{quantile}{quantiles of simulated parameters}
\item{maxlik}{calculated maximum logged likelihood value}
}

\references{
Lunn D, Jackson C, Best N, Thomas A, Spiegelhalter D, 2013. The BUGS book: a practical introduction to bayesian analysis. Chapman & Hall/CRC Press.

Gelman A, Carlin JB, Stern HS, Dunson DB, Vehtari A, Rubin DB, 2013. Bayesian data analysis. Chapman & Hall/CRC Press.
}

\seealso{
  \link{mcFMM}; \link{mcMAM}
}

\keyword{statistical age models}
\keyword{equivalent dose}
