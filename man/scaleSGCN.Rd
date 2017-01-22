\name{scaleSGCN}
\alias{scaleSGCN}
\alias{scaleSGCN.default}

\title{
Natural-dose signal re-scaling
}

\description{
Re-scaling standardised natural-dose signals according to the "global standardised growth curve" (gSGC) method suggested by Li et al. (2015).
}

\usage{
scaleSGCN(SGCdata, pars, model, origin)
}

\arguments{
    \item{SGCdata}{\link{data.frame}(\bold{required}): standardised natural-dose signal and regenerative-dose signal data used for 
natural-dose signal re-scaling, it should contain five columns (i.e., \code{NO}, \code{SAR.Cycle}, \code{Dose}, \code{Signal}, and \code{Signal.Err}), see \link{SARdata} for details}

  \item{pars}{\link{vector}(\bold{required}): optimized parameters of the SGC obtained using\cr 
function \link{lsNORM}}

  \item{model}{\link{character}(\bold{required}): fitting model used for obtaining \code{pars}}

  \item{origin}{\link{logical}(\bold{required}): logical value indicating if established SGC passes the origin}
}

\details{
Standardised natural-dose signals are re-scaled according to \bold{Eqn.(10)} of Li et al. (2015).
}

\value{
Return a \link{matrix} containing re-scaled standardised natural-dose signals.
}

\references{
Li B, Roberts RG, Jacobs Z, Li SH, 2015. Potential of establishing a "global standardised growth curve" (gSGC) for optical dating of quartz from sediments. Quaternary Geochronology, 27: 94-104.

Li B, Jacobs Z, Roberts RG, 2016. Investigation of the applicability of standardised growth curves for OSL dating of quartz from Haua Fteah cave, Libya. Quaternary Geochronology, 35: 1-15.
}

\seealso{
  \link{lsNORM}; \link{calSGCED}
}

\examples{
 # Not run.
 # data(SARdata)
 # gSGCdata <- SARdata[is.element(SARdata[,2], c("N","R3")),]
 # gSGCpars <- c(137.440874251, 0.007997863, 2.462035263, -0.321536177)
 # scaleSGCN(gSGCdata, gSGCpars, model="gok", origin=FALSE)
}

\keyword{LS-normalisation}
\keyword{equivalent dose}
