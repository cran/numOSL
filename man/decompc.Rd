\name{decompc}
\alias{decompc}
\alias{decompc.default}
\title{
Decomposition decay curve (with constant component subtracted) in OSL dating
}
\description{
Decomposing the OSL decay curve (type "CW" or "LM") to a specified number of components (plus a constant component).
}
\usage{
decompc(Sigdata, ncomp = 2, typ = c("cw", "lm"), 
        control.args = list(), LEDpower = 60, LEDwavelength = 470, 
        plot = TRUE, xlog = TRUE, lwd = 3, samplename = NULL, outfile = NULL)
}
\arguments{
  \item{Sigdata}{\link{data.frame}(\bold{required}): two columns, time (t1,t2,etc) and correspond signal values (s1,s2,etc) to be decomposed
}
  \item{ncomp}{\link{numeric}(with default): the number of components to be decomposed, which must be an integer in the space [1,7]
}
  \item{typ}{\link{character}(with default): type ("cw" or "lm") of the OSL decay curve, default typ="cw"
}
  \item{control.args}{\link{list}(with default): further arguments used in differential evolution algorithm, see function \link{decomp}
}
  \item{LEDpower}{\link{numeric}(with default): LED power (mW/cm2) used for intensity ramping during the measurement
}
  \item{LEDwavelength}{\link{numeric}(with default): LED wavelength (nm) used for stimulation during the measurement
}
  \item{plot}{\link{logical}(with default): whether plot the fitted results or not
}
  \item{xlog}{\link{logical}(with default): tansform x-axis to log-scale for plotting or not
}
  \item{lwd}{\link{numeric}(with default): the width of curves for plotting
}
  \item{samplename}{\link{character}(optional): name of the sample
}
  \item{outfile}{\link{character}(optional): if specified, decomposed signal values will be written to a file of name outfile in .csv format and saved to the current work directory
}
}
\details{
Function decompc() is used to decompose decay curves with a subtracted constant component (the constant component could be regarded as the summation of some slow components and the background). For CW-OSL decay curves:\cr
I(t)=a1*b1*exp(-b1*t)+a2*b2*exp(-b2*t)+...+ak*bk*exp(-bk*t)+c , k=1, 2, ..., 7, \cr\cr
And for LM-OSL decay curves:\cr
I(t)=a1*b1*(t/P)*exp(-b1*t^2/(2*P))+a2*b2*(t/P)*exp(-b2*t^2/(2*P))+...+\cr
ak*bk*(t/P)*exp(-bk*t^2/(2*P))+c*(t/P) , k=1, 2, ..., 7\cr
Note that unlike function decomp(), transformation (transf) of parameters is not allowed in this function. For more details, see function \link{decomp}.
}
\value{
Return an invisible list containing following elements:\cr
  \item{Comp.Signal}{fitted values for observations and each component}
  \item{pars}{estimated parameters and associated standard errors}
  \item{constant}{estimated constant component (a constant component for CW-OSL and an increased linearly component for LM-OSL)}
  \item{value}{sum of squared residuals}
  \item{errorflag}{error message generated during the calculation, if parameters' standard errors can be estimated, errorflag=0, else 1}
}
\seealso{
  \link{decomp}; \link{fastED}
}
\examples{
  data(Signaldata)
  res<-decompc(Signaldata$lm1,ncomp=5,typ="lm")
  list(res$pars,res$constant)
}
\keyword{CW-OSL decay curve}
\keyword{LM-OSL decay curve}
\keyword{decay curve decomposition}
