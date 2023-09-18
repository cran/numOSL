\name{lsNORM}
\alias{lsNORM}
\alias{lsNORM.default}

\title{
Regenerative-dose signal optimization using the LS-normalisation procedure
}
\description{
Optimizing standardised regenerative-dose signals according to the least-squares normalisation (LS-normalisation) procedure using an iterative scaling and fitting procedure proposed by Li et al. (2016).
}
\usage{
lsNORM(SARdata, model = "gok", origin = FALSE, 
       weight = FALSE, natural.rm = TRUE, 
       norm.dose = NULL, maxiter = 10, 
       plot = TRUE)
}

\arguments{
  \item{SARdata}{\link{matrix}(\bold{required}): SAR data used for performing the LS-normalisation\cr
procedure, it should contain five columns (i.e., \code{NO}, \code{SAR.Cycle}, \code{Dose},\cr
\code{Signal}, and \code{Signal.Err}), see \link{SARdata} for details}

  \item{model}{\link{character}(with default): model used for growth curve fitting, see \link{fitGrowth} for available models}

  \item{origin}{\link{logical}(with default): logical value indicating if the growth curve should be forced to pass the origin}

  \item{weight}{\link{logical}(with default): logical value indicating if the growth curve should be fitted using a weighted procedure, see function \link{fitGrowth} for details}

  \item{natural.rm}{\link{logical}(with default): logical value indicating if the standardised natural-dose signal 
should be removed from \code{SARdata}}

  \item{norm.dose}{\link{numeric}(optional): regenerative-dose used for re-scaling established gSGC. 
For example, if \code{norm.dose=100}, then the signal value for a dose value of 100 (Gy|s) will be re-scaled to unity}

  \item{maxiter}{\link{integer}(with default): allowed maximum number of iterations during the\cr
LS-normalisation optimization process}

  \item{plot}{\link{logical}(with default): logical value indicating if the results should be plotted}
}

\details{
Function \link{lsNORM} is used for optimizing regenerative-dose signal data from a number of grains (aliquots) according to the least-squares normalisation (LS-normalisation) procedure outlined by Li et al. (2016) using an iterative scaling and fitting procedure.\cr\cr
The LS-normalisation procedure for growth curve optimization involves the following steps:\cr
(1) Fit standardised regenerative-dose signals from all of the aliquots;\cr
(2) Re-scale the individual growth curve from each aliquot using a scaling factor. The scaling factor for each aliquot is determined in a way such that the sum of squared residuals is minimized. Each aliquots is treated individually, and different scaling factors are calculated for different aliquots.\cr
(3) Iterate the fitting (step 1) and re-scaling (step 2). The iterative procedure is performed repeatedly until there is negligible change in the relative standard deviation of the normalised growth curve data.
}
\value{
  Return an invisible \link{list} that contains the following elements:
  \item{norm.SARdata}{SAR data sets optimized using the LS-normalisation procedure}
  \item{sf}{scaling factor of standardised regenerative-dose signals}
  \item{iter}{number of iterations required}
  \item{LMpars1}{optimized parameters for the growth curve before LS-normalisation}
  \item{value1}{minimized objective for the growth curve before LS-normalisation}
  \item{avg.error1}{average fit error for the growth curve before LS-normalisation}
  \item{RCS1}{reduced chi-square value for the growth curve before LS-normalisation}
  \item{FOM1}{figure of merit value for the growth curve before LS-normalisation in percent}
  \item{LMpars2}{optimized parameters for the growth curve after LS-normalisation}
  \item{value2}{minimized objective for the growth curve after LS-normalisation}
  \item{avg.error2}{average fit error for the growth curve after LS-normalisation}
  \item{RCS2}{reduced chi-square value for the growth curve after LS-normalisation}
  \item{FOM2}{figure of merit value for the growth curve after LS-normalisation in percent} 
}

\references{
Li B, Jacobs Z, Roberts RG, 2016. Investigation of the applicability of standardised growth curves for OSL dating of quartz from Haua Fteah cave, Libya. Quaternary Geochronology, 35: 1-15.
}

\seealso{
  \link{analyseBINdata}; \link{fitGrowth}; \link{SARdata}; \link{scaleSGCN}; \link{calSGCED}
}

\examples{
  ### Example 1:
  data(SARdata)
  # Use only the first five aliquots of SARdata.
  Data <- SARdata[1:40,]
  res_lsNORM <- lsNORM(Data, model="gok")
  res_lsNORM$norm.SARdata

  ### Example 2 (not run):
  # data(BIN)
  # obj_pickBIN <- pickBINdata(BIN, Position=1:48, Grain=0,
  #                            LType="OSL", view=FALSE)
  # obj_analyseBIN <- analyseBINdata(obj_pickBIN, nfchn=3, nlchn=20)
  # lsNORM(obj_analyseBIN$SARdata, norm.dose=300)
  
  
}
\keyword{growth curve}
\keyword{LS-normalisation}
