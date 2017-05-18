\name{scaleSGCN}
\alias{scaleSGCN}
\alias{scaleSGCN.default}

\title{
Natural-dose signal re-scaling
}

\description{
Re-scaling Sensitivity-corrected natural-dose signals according to the "global standardised growth curve" (gSGC) method 
suggested by Li et al. (2015).
}

\usage{
scaleSGCN(obj_analyseBIN, SGCpars, model, origin, SAR.Cycle, 
          Tn.above.3BG = TRUE, TnBG.ratio.low = NULL, rseTn.up = NULL, 
          FR.low = NULL, use.se = TRUE, outfile = NULL)
}

\arguments{
  \item{obj_analyseBIN}{\link{list}(\bold{required}): an object of S3 class "analyseBIN" produced by\cr 
function \link{analyseBINdata} or \link{as_analyseBIN}}

  \item{SGCpars}{\link{vector}(\bold{required}): optimized parameters of the SGC obtained using function \link{fitGrowth} or \link{lsNORM}}

  \item{model}{\link{character}(\bold{required}): fitting model used for obtaining \code{SGCpars}}

  \item{origin}{\link{logical}(\bold{required}): logical value indicating if established SGC passes the origin}

  \item{SAR.Cycle}{\link{character}(with default): a two-element character vector containing SAR cycles used for 
natural-dose signal re-scaling. Example: \code{SAR.Cycle=c("N","R3")}}

   \item{Tn.above.3BG}{\link{logical}(with default): logical value indicating if aliquot (grain) with Tn below 3 sigma BG should be rejected}

  \item{TnBG.ratio.low}{\link{numeric}(optional): lower limit on ratio of initial Tn signal to BG}

  \item{rseTn.up}{\link{numeric}(optional): upper limit on relative standard error of Tn in percent}

  \item{FR.low}{\link{numeric}(optional): lower limit on fast ratio of Tn}

  \item{use.se}{\link{logical}(with default): logical value indicating if standard errors of values should be used during application of rejection criteria}

  \item{outfile}{\link{character}(optional): if specified, scaled SGC data related quantities will be written to a CSV file named \code{"outfile"} and saved to the current work directory}
}

\details{
Sensitivity-corrected natural-dose signals are re-scaled according to \bold{Eqn.(10)} of Li et al. (2015).
}

\value{
Return an invisible \link{list} that contains the following elements:
  \item{scale.Ltx}{scaled natural-dose signals and associated standard errors}
  \item{agID}{aliquot (grain) ID of scaled natural-dose signals}
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
 data(SARdata)
 gSGCpars <- c(137.440874251, 0.007997863, 2.462035263, -0.321536177)
 scaleSGCN(as_analyseBIN(SARdata), SGCpars=gSGCpars, model="gok", 
           origin=FALSE, SAR.Cycle=c("N","R3"))
}

\keyword{LS-normalisation}
\keyword{equivalent dose}
