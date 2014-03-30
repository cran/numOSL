\name{decompc}
\alias{decompc}
\alias{decompc.default}
\title{
Decomposing an OSL decay curve (allows a constant component being subtracted)
}
\description{
Decomposing an OSL decay curve (type "CW" or "LM") to a specified number of components (plus a constant component).
}
\usage{
decompc(Sigdata, ncomp = 2, typ = c("cw", "lm"), 
        control.args = list(), LEDpower = 60, LEDwavelength = 470, 
        plot = TRUE, xlog = TRUE, lwd = 3, samplename = NULL, outfile = NULL)
}
\arguments{
  \item{Sigdata}{\link{data.frame}(\bold{required}): two columns, i.e. time values (t1,t2,etc) and signal values (s1,s2,etc)}
  \item{ncomp}{\link{numeric}(with default): the number of components to be decomposed, which should be an integer from the space [1,7]}
  \item{typ}{\link{character}(with default): the type ("cw" or "lm") of an OSL decay curve, default typ="cw"
}
  \item{control.args}{\link{list}(with default): some further arguments used in the differential evolution algorithm, see function \link{decomp}}
  \item{LEDpower}{\link{numeric}(with default): the LED power (mW/cm2) used in the measurement}
  \item{LEDwavelength}{\link{numeric}(with default): the wavelength (nm) used in the measurement}
  \item{plot}{\link{logical}(with default): plot the fitted results or not}
  \item{xlog}{\link{logical}(with default): tansform x-axis to log-scale or not}
  \item{lwd}{\link{numeric}(with default): the width of curves in the plot}
  \item{samplename}{\link{character}(optional): the name of the sample}
  \item{outfile}{\link{character}(optional): if specified, decomposed signal values will be written to a file named "outfile" in ".csv" format and will be saved to the current work directory}
}
\details{
Function decompc() is used to decompose a decay curves with a subtracted constant component (a constant component may be regarded as the sum of signals of some slower components plus the background). For a CW-OSL decay curves:\cr
I(t)=a1*b1*exp(-b1*t)+a2*b2*exp(-b2*t)+...+ak*bk*exp(-bk*t)+c , k=1, 2, ..., 7, \cr\cr
And for LM-OSL decay curves:\cr
I(t)=a1*b1*(t/P)*exp(-b1*t^2/(2*P))+a2*b2*(t/P)*exp(-b2*t^2/(2*P))+...+\cr
ak*bk*(t/P)*exp(-bk*t^2/(2*P))+c*(t/P) , k=1, 2, ..., 7\cr
Note that unlike function decomp(), the fitting model is different and a transformation (transf) of parameters is not allowed here. For more details, see function \link{decomp}.
}
\value{
Return an invisible list containing following elements:\cr
  \item{Comp.Signal}{the fitted decaying signals for each component}
  \item{pars}{the estimates and associated standard errors}
  \item{constant}{the estimated constant component}
  \item{value}{the sum of squared residuals}
  \item{errorflag}{an error message generated during the calculation, if parameters' standard errors are available, errorflag=0, or else errorflag=1}
}
\seealso{
  \link{decomp}; \link{fastED}
}
\examples{
  data(Signaldata)
  res<-decompc(Signaldata$lm1,ncomp=5,typ="lm")
  list(res$pars,res$constant)
}
\keyword{fitting OSL decay curve}
