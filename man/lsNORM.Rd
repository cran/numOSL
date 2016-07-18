\name{lsNORM}
\alias{lsNORM}
\alias{lsNORM.default}

\title{
Regenerative dose signal data optimization using the least-squares normalisation (LS-normalisation) procedure
}
\description{
Optimizing regenerative dose signal data according to the least-squares normalisation (LS-normalisation) procedure using an iterative scaling and fitting procedure proposed by Li et al. (2016).
}
\usage{
lsNORM(Data, model = "gok", origin = FALSE,
       maxiter = 10, weight = TRUE, plot = TRUE)
}

\arguments{
  \item{Data}{\link{matrix}(\bold{required}): data used for performing the LS-normalisation procedure, it should contain five columns (i.e., Grain.NO, SAR.Cycle, Dose, Signal, and Signal.Err), see \link{SARdata} for details}

  \item{model}{\link{character}(with default): model used for growth curve fitting, see \link{fitGrowth} for available models}

  \item{origin}{\link{logical}(optional): force the growth curve to pass the origin or not}

  \item{maxiter}{\link{integer}(with default): allowed maximum number of iterations during the LS-normalisation optimization process}

  \item{weight}{\link{logical}(with default): fit the growth curve using a weighted or an un-weighted procedure, see function \link{fitGrowth} for details}

  \item{plot}{\link{logical}(with default): plot the output or not}
}

\details{
Function \link{lsNORM} is used for optimizing regenerative dose signal data from a number of grains (aliquots) according to the least-squares normalisation (LS-normalisation) procedure outlined by Li et al. (2016) using an iterative scaling and fitting procedure.\cr\cr
The LS-normalisation procedure for growth curve optimization involves the following steps:\cr
(1) Fit standardised regenerative dose signals from all of the aliquots;\cr
(2) Re-scale the individual growth curve from each aliquot using a scaling factor. The scaling factor for each aliquot is determined in a way such that the sum of squared residuals is minimized. Each aliquots is treated individually, and different scaling factors are calculated for different aliquots.\cr
(3) Iterate the fitting (step 1) and re-scaling (step 2). The iterative procedure is performed repeatedly until there is negligible change in the relative standard deviation of the normalised growth curve data.
}
\value{
  Return a list that contains the following elements:
  \item{optData}{regenerative dose signals optimized using the LS-normalisation procedure}
  \item{sf}{total scaling factor for regenerative dose signals for each aliquot (grain)}
  \item{iter}{number of iterations required}
  \item{LMpars1}{optimized parameters for the fitted growth curve using original regenerative dose signals}
  \item{value1}{minimized objective for the fitted growth curve using original regenerative dose signals}
  \item{LMpars2}{optimized parameters for the fitted growth curve using LS-normalised regenerative dose signals}
  \item{value2}{minimized objective for the fitted growth curve using LS-normalised regenerative dose signals}
}

\references{
Li B, Jacobs Z, Roberts RG, 2016. Investigation of the applicability of standardised growth curves for OSL dating of quartz from Haua Fteah cave, Libya. Quaternary Geochronology, 35: 1-15.
}

\seealso{
  \link{fitGrowth}; \link{calSGCED}
}

\examples{
  data(SARdata)
  # Use only the first five aliquots of SARdata.
  Data <- SARdata[1:40,]
  lsNORM(Data, model="gok")
}
\keyword{growth curve}
\keyword{LS-normalisation procedure}
