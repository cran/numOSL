\name{reportSAM}
\alias{reportSAM}
\alias{reportSAM.default}
\title{
Report the results of statistical age models
}
\description{
Return a plot and a list that summaries the posterior density distributions of the parameters, the history of simulations and the auto-correlations of those parameters. Also return the inferences.
}
\usage{
reportSAM(x, burn = NULL, thin = NULL, plot = TRUE, outfile = NULL, ...)
}
\arguments{
  \item{x}{\bold{required}: an object of class "mcAgeModels", which is produced by function \link{mcFMM} or \link{mcMAM} }
  \item{burn}{\link{numeric}(with default): the number of simulations to be discarded, i.e. the initial, non-stationary portion of the chains, default burn=floor(nsim/5L)}
  \item{thin}{\link{numeric}(with default): only one value from each portion that constains every thin values is retained, default thin=5}
  \item{plot}{\link{logical}(with default): plot the results or not}
  \item{outfile}{\link{character}(optional): if specified, simulated chains will be written to a file named "outfile" in ".csv" format and will be saved to the current work directory}
  \item{\dots}{do not use}
}
\details{
The numbers of "burn-in" and "thining" are used to refine the simulated chains in order to obtain more reliable inference. The initial portion of the chains may not sufficiently close to stationarity and
need to be discarded. Where the "thining" protocol is used to reduce the posterior correlation between parameters.
}
\value{
Return list that contains following elements:
\item{pars}{estimates and their standard errors}
\item{maxlik}{the corresponding logged maximum likelihood value}
}
\references{
Lunn, D., Jackson, C., Best, N., Thomas, A., Spiegelhalter, D., 2013. The BUGS Book: A Practical Introduction to Bayesian Analysis. Journal of Statistical Software,
52, pp. 1-4.
}
\seealso{
  \link{mcFMM}; \link{mcMAM}
}
\keyword{statistical age models}
\keyword{Slice sampling}
