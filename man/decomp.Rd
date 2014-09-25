\name{decomp}
\alias{decomp}
\alias{decomp.default}
\title{
Decomposing an OSL decay curve
}
\description{
Decomposing an OSL decay curve (type "CW" or "LM") to a given number of first-order exponential components using a combination of differential evolution and Levenberg-Marquardt algorithm.
}

\usage{
decomp(Sigdata, ncomp = 2, constant = TRUE, 
       typ = c("cw","lm"), control.args = list(), 
       transf = TRUE, LEDpower = 60, LEDwavelength = 470, 
       plot = TRUE, xylog= FALSE, lwd = 3, outfile = NULL)
}

\arguments{
  \item{Sigdata}{\link{matrix}(\bold{required}): two columns, i.e., time values and photon counts}

  \item{ncomp}{\link{numeric}(with default): the number of components to be decomposed}

  \item{constant}{\link{logical}(with default): subtract a constant component or not}

  \item{typ}{\link{character}(with default): type (\code{"cw"} or \code{"lm"}) of a decay curve, default \code{typ="cw"}}

  \item{control.args}{\link{list}(with default): arguments used by the differential evolution algorithm, see details}

  \item{transf}{\link{logical}(with default): do not use}

  \item{LEDpower}{\link{numeric}(with default): LED power (mW/cm2) of the instrument}

  \item{LEDwavelength}{\link{numeric}(with default): LED wavelength (nm) of the instrument}

  \item{plot}{\link{logical}(with default): draw a plot or not}

  \item{xylog}{\link{logical}(with default): using both logged x-axis and y-axis or only logged x-axis}

  \item{lwd}{\link{numeric}(with default): the width of a line}

  \item{outfile}{\link{character}(optional): if specified, decomposed signal values will be written to a file named \code{"outfile"} in \code{".csv"} format and saved to the current work directory}
}

\details{
This function decomposes an OSL decay curve to a specified number of components using a combination of differential evolution and Levenberg-Marquardt algorithm. Both OSL decay curves of type "CW" and "LM" can be decomposed.\cr\cr
For an CW-OSL decay curve, the fitting model (Bluszcz and Adamiec, 2006) is:\cr
\code{I(t)=a1*b1*exp(-b1*t)+...+ak*bk*exp(-bk*t)},\cr
where \code{k=1, 2, ..., 7}, \code{I(t)} is the luminescence intensity as a function of time, \code{a} is the number of trapped electrons, and \code{b} is the detrapping rate. The constant component will be \code{c} if \code{constant=TRUE}.\cr\cr
For an LM-OSL decay curve, the fitting model (Bulur, 2000) is:\cr
\code{I(t)=a1*b1*(t/P)*exp[-b1*t^2/(2*P)]+...+ak*bk*(t/P)*exp[-bk*t^2/(2*P)]},\cr
where \code{k=1, 2, ..., 7}, and \code{I(t)} is the luminescence intensity as a function of time, \code{P} is the total stimulation time, \code{a} is the number of trapped electrons, and \code{b} is the detrapping rate. The constant component will be \code{c*(t/P)} if 
\code{constant=TRUE}.\cr\cr
Parameters are initialized using a differential evolution method suggested by Bluszcz and Adamiec (2006), then the Levenberg-Marquardt algorithm (available at \url{http://ftp.netlib.org/minpack}) will be performed to fit the curve. Arguments (\code{control.args}) that control the performance of the differential evolution algorithm include:\cr\cr
\emph{factor}: control the number of parents, \code{np=factor*ncomp}, default \code{factor=10} \cr\cr
\emph{f}: a weighting factor that lies between \code{0} and \code{1.2}, default \code{f=0.5} \cr\cr
\emph{cr}: a crossover constant that lies between \code{0} and \code{1}, default \code{cr=0.99} \cr\cr
\emph{maxiter}: the allowed maximum number of iterations, default \code{maxiter=500} \cr\cr
\emph{tol}: a tolerance for stopping the iteration, the iteration will be terminated if the relative standard deviations of parameters are smaller than \code{tol}, defalut \code{tol=0.1}\cr
}

\value{
Return an list containing following elements:
  \item{pars}{the decay rates, the number of trapped electrons and their standard errors}
  \item{constant}{the constant component (it will be \code{0} if \code{constant=FALSE)}}
  \item{value}{the minimized sum of squared residuals}
}

\seealso{
  \link{fastED}
}

\references{
Bluszcz A, 1996. Exponential function fitting to TL growth data and similar applications. Geochronometria, 13: 135-141.

Bluszcz A, Adamiec G, 2006. Application of differential evolution to fitting OSL decay curves. Radiation Measurements, 41: 886-891.

Bulur E, 2000. A simple transformation for converting CW-OSL curves to LM-OSL curves. Radiation Measurements, 32: 141-145.

Differential evolution algorithm, \url{http://en.wikipedia.org/wiki/Differential_evolution}
  
Jain M, Murray AS, Boetter-Jensen L, 2003. Characterisation of blue-light stimulated luminescence components in different quartz samples: implications for dose measurement. Radiation Measurements, 37: 441-449.

More JJ, 1978. "The Levenberg-Marquardt algorithm: implementation and theory," in Lecture Notes in Mathematics: Numerical Analysis, 
Springer-Verlag: Berlin. 105-116.

\bold{Further reading} 

Adamiec G, 2005. OSL decay curves-relationship between single- and multiple-grain aliquots. Radiation Measurements, 39: 63-75.

Choi JH, Duller GAT, Wintle AG, 2006. Analysis of quartz LM-OSL curves. Ancient TL, 24: 9-20.

Li SH, Li B, 2006. Dose measurement using the fast component of LM-OSL signals from quartz. Radiation Measurements, 41: 534-541.
}


\note{
Arugments "LEDpower" and "LEDwavelength" are used for determining the absolute photoionisation cross section. However, the estimated absolute photoionisation cross section (\code{pcs}) may not be so reliable as it seems, for both the stimulation power and the wavelength vary with the age of the instrument.\cr\cr
For a given model, this routine will return an error if any standard errors of parameters of the decay curve cannot be estimated by numerical approximation.\cr\cr
Function \emph{decompc} in previous version was bundled to function \emph{decomp}.\cr\cr
I thank Andrzej Bluszcz, who helps me a lot during the programming of this function. Jeong-Heon Choi and Bo Li are thanked for providing published data to test this routine.
}

\examples{
 data(Signaldata)
 decomp(Signaldata$lm,ncomp=3,typ="lm",control.args=list(maxiter=10))
}
\keyword{OSL decay curve}
