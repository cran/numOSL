\name{reportMC}
\alias{reportMC}
\alias{reportMC.default}
\title{
Reporting MCMC outputs for statistical age models
}
\description{
Summarizing distributions of parameters simulated from statistical age models using a Markov Chain Monte Carlo method. 
}
\usage{
reportMC(obj, burn = 10000, thin = 5, 
         plot = TRUE, outfile = NULL, ...)
}
\arguments{
  \item{obj}{\link{list}(\bold{required}): an object of S3 class \code{"mcAgeModels"}, which is produced by function \link{mcFMM} or \link{mcMAM} }

  \item{burn}{\link{integer}(with default): number of iterations (i.e., the initial, non-stationary\cr
portion of the chain) to be discarded}

  \item{thin}{\link{integer}(with default): take every \code{thin-th} iteration}

  \item{plot}{\link{logical}(with default): plot the MCMC output or not}

  \item{outfile}{\link{character}(optional): if specified, simulated parameters will be written to a CSV file named \code{"outfile"} and saved to the current work directory}

  \item{\dots}{do not use}
}

\details{
Function \link{reportMC} summarizes the output of a Markov Chain (the mean values, the standard deviations, the mode values, and the 2.5, 25, 50, 75, 97.5 quantiles of the simulated parameters). The initial \emph{i} (\code{burn=i}) samples may have been affected by the inital state and has to be discarded ("burn-in"). Autocorrelation of simulated samples can be reduced by taking every \emph{j-th} (\code{thin=j}) iteration ("thining").
}

\value{
Return a \link{list} that contains the following elements:
\item{pars}{means, standard deviations, and modes of simulated parameters}
\item{quantile}{quantiles of simulated parameters}
\item{maxlik}{maximum logged likelihood values calculated using the means and modes of simulated parameters}
\item{bic}{Bayesian Information Criterion (BIC) values calculated using the means and modes of simulated parameters}
}

\references{
Lunn D, Jackson C, Best N, Thomas A, Spiegelhalter D, 2013. The BUGS book: a practical introduction to bayesian analysis. Chapman & Hall/CRC Press.

Gelman A, Carlin JB, Stern HS, Dunson DB, Vehtari A, Rubin DB, 2013. Bayesian data analysis. Chapman & Hall/CRC Press.

Peng J, Dong ZB, Han FQ, 2016. Application of slice sampling method for optimizing OSL age models used for equivalent dose determination. Progress in Geography, 35(1): 78-88. (In Chinese).
}

\seealso{
  \link{mcFMM}; \link{mcMAM}
}

\keyword{statistical age models}
\keyword{Central Age Model}
\keyword{Finite Mixture Age Model}
\keyword{Minimum Age Model}
\keyword{equivalent dose}
\keyword{Markov chain Monte Carlo}
