\name{sgcED}
\alias{sgcED}
\alias{sgcED.default}
\title{
Fitting a dose response curve and calculating a number of equivalent dose values using the standardised growth curves (SGC) method
}
\description{
It is a function for fitting a dose response curve and calculating an equivalent dose values using the standardised growth curves (SGC) method after Roberts HM and Duller GAT (2004). It is based on the fact that using information collected as part of standard SAR measurements, it is possible to construct reproducible patterns of growth from multiple-grain aliquots of a large number of samples (Roberts HM and Duller GAT, 2004).
}
\usage{
sgcED(Curvedata, Ltx, 
      model = c("line","exp","line+exp"), origin = FALSE,
      nstart = 100, upb = 1, ErrorMethod = c("mc","sp"),
      nsim = 1000, plot = TRUE, samplename = NULL, outfile = NULL)
}
\arguments{
  \item{Curvedata}{\link{data.frame}(\bold{required}): three columns, see function \link{calED}}
  \item{Ltx}{\link{data.frame}(\bold{required}): two columns (standardized signals and its standard errors)}
  \item{model}{\link{character}(with default): a fitting model("linear", "exponential" or "linear+exponential"), see function \link{calED}}
  \item{origin}{\link{logical}(optional): force the fitting to pass the origin (x=0,y=0) or not}
  \item{nstart}{\link{numeric}(with default): the maximum number of attempts used to initialize parameters when fitting a dose response curve, see function \link{calED}}
  \item{upb}{\link{numeric}(with default): the upper limit on b value, see function \link{calED}}
  \item{ErrorMethod}{\link{character}(with default): a method ("sp" or "mc") for estimating the standard errors of equivalent doses, see function \link{calED}}
  \item{nsim}{\link{numeric}(with default): the allowed maximum number of simulations if a Monte Carlo method is used for error assessing}
  \item{plot}{\link{logical}(with default): plot the results or not}
  \item{samplename}{\link{character}(optional): the name of the sample}
  \item{outfile}{\link{character}(optional): if specified, calculated equivalent doses will be written to a file named "outfile" in ".csv" format will and be saved to the current work directory}
}
\details{
see function \link{calED}() for details.
}
\value{
  Return a visible list that contains following elements:
  \item{LMpars}{the characteristic parameters and their standard errors of the fitted dose response curve}
  \item{residual}{the sum of square of residuals}
  \item{fit.value}{the fitted values}
  \item{ED}{the estimated equivalent doses and their standard errors}
}
\references{
Roberts, H.M. and Duller, G.A.T., 2004. Standardised growth curves for optical dating of sediment using multiple-grain aliquots. Radiation Measurements 38(2), pp. 241-252.

Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived from single aliquot regenerative dose measurements. Ancient TL 25(1), pp. 15-24.

Roberts, H.M., Durcan, J.A., Duller, G.A.T., 2009. Exploring procedures for the rapid assessment of optically stimulated luminescence range-finder ages. Radiation Measurements, 44(5-6), pp. 582-587.

Jorge More, Burton Garbow, Kenneth Hillstrom, User Guide for MINPACK-1, Technical Report ANL-80-74, Argonne National Laboratory, 1980.

\bold{Further reading} 

Long, H., Lai, Z.P., Fan, Q.S., Sun, Y.J., Liu, X.J., 2010. Applicability of a quartz OSL standardised growth curve for De determination up to
400 Gy for lacustrine sediments from the Qaidam Basin of the Qinghai-Tibetan Plateau. Quaternary Geochronology 5(2-3), pp. 212-217.

Yang, L.H., Lai, Z.P., Long, H., Zhang, J.R., 2011. Construction of a quartz OSL standardised growth curve (SGC) for aeolian samples from the Horqin dunefield in northeastern China. Geochronometria 38(4), pp. 391-396.
}
\seealso{
  \link{calED}
}
\examples{
 ltx1<-c(0.032,1.61,2.55,3.21,3.87,0.031,1.55) # Lx/Tx for the first aliquot
 ltx2<-c(0.025,1.44,2.47,3.35,4.17,0.033,1.47) # Lx/Tx for the second aliquot
 ltx3<-c(0.027,1.51,2.68,3.52,4.41,0.021,1.39) # Lx/Tx for the third aliquot
 ltx4<-c(0.018,1.71,2.28,3.81,4.03,0.017,1.62) # Lx/Tx for the four aliquot
 ltx5<-c(0.026,1.49,1.99,3.43,4.17,0.015,2.01) # Lx/Tx for the five aliquot
 ltx<-cbind(ltx1,ltx2,ltx3,ltx4,ltx5)
 ltx<-cbind(apply(ltx,MARGIN=1,mean),
            apply(ltx,MARGIN=1,sd))  # means an standard deviations
 redose<-c(0,18,36,54,72,0,18)       # the same ReDose for the five aliquots
 Curvedata<-data.frame(redose,ltx)
 Ltx<-data.frame(c(0.5,1.0,1.8,2.3,2.8,3.1,3.6,4.0),
                 rep(0.1,8))         # Lx/Txs from which EDs are expected
 sgcED(Curvedata,Ltx,model="line+exp",origin=TRUE) # fitting y=a*(1-exp(-b*x))+c*x
}
\keyword{equivalent dose}
