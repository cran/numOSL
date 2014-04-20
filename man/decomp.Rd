\name{decomp}
\alias{decomp}
\alias{decomp.default}
\title{
Decomposing an OSL decay curve
}
\description{
Function decomp() decomposes an OSL decay curve (type "CW" or "LM") to a specified number of components.
}
\usage{
decomp(Sigdata, ncomp = 3, typ = c("cw","lm"),
       control.args = list(), transf = FALSE, 
       LEDpower = 60, LEDwavelength = 470, plot = TRUE,
       xlog= TRUE, lwd = 3, samplename = NULL, outfile = NULL)
}

\arguments{
  \item{Sigdata}{\link{data.frame}(\bold{required}): two columns, i.e. time values (t1,t2,etc) and corresponded signal values (s1,s2,etc)}
  \item{ncomp}{\link{numeric}(with default): the number of components to be decomposed, which should be an integer from the space [1,7]}
  \item{typ}{\link{character}(with default): type ("cw" or "lm") of a decay curve, default typ="cw"}
  \item{control.args}{\link{list}(with default): some further arguments used in the differential evolution algorithm, see details}
  \item{transf}{\link{logical}(with default): transform estimated parameters or not, see details}
  \item{LEDpower}{\link{numeric}(with default): LED power (mW/cm2) used in the measurement}
  \item{LEDwavelength}{\link{numeric}(with default): LED wavelength (nm) used in the measurement}
  \item{plot}{\link{logical}(with default): plot the fitted results or not}
  \item{xlog}{\link{logical}(with default): use a log-scale x-axis for plotting or not}
  \item{lwd}{\link{numeric}(with default): the width of curves in the plot}
  \item{samplename}{\link{character}(optional): the name of the sample}
  \item{outfile}{\link{character}(optional): if specified, decomposed signal values will be written to a file named "outfile" in ".csv" format and will be saved to the current work directory}
}
\details{
This function decomposes an OSL decay curve to a specified number of components using a combination of differential evolution and Levenberg-Marquardt algorithm. Both OSL decay curves of type "CW" and "LM" can be decomposed.\cr
For an OSL decay curve of type "CW", the fitting model (accroding to Bluszcz and Adamiec, 2006) is:\cr\cr
I(t)=a1*exp(-b1*t)+a2*exp(-b2*t)+...+ak*exp(-bk*t) , k=1, 2, ..., 7, which is identical to:\cr\cr
I(t)=a1*b1*exp(-b1*t)+a2*b2*exp(-b2*t)+...+ak*bk*exp(-bk*t) , k=1, 2, ..., 7 if "transf" is set to be TRUE.\cr\cr
where I(t) is the luminescence intensity as a function of time(t), a is the number of trapped electrons, b is the detrapping rates. Parameters are initialized using a differential evolution method suggested by Bluszcz A and Adamiec G (2006), then the Levenberg-Marquardt algorithm (minpack: original FORTRAN 77 version by Jorge More, Burton Garbow, Kenneth Hillstrom. FORTRAN 90 version by John Burkardt) will be performed to obtain the final estimate. In some radical cases, the differential evolution or the Levenberg-Marquardt algorithm may fails, then a serise of simple attempts will be carried out with the Levenberg-Marquardt algorithm. Arguments (control.args) that employed to control the performance of a differential evolution algorithm include:\cr\cr
\bold{factor}: control the number of parents (NP), np=factor*ncomp, default factor=20; \cr
\bold{f}: a weighting factor that lies between 0 and 1.2, default f=0.5; \cr
\bold{cr}: a crossover constant that lies between 0 and 1, default cr=0.9; \cr
\bold{maxiter}: the allowed maximum number of iterations, default maxiter=1e3; \cr
\bold{tol}: a tolerance for stopping the iteration, the iteration will be terminated if magnitudes of standard deviations of parameters are smaller than tol. default tol=0.1 \cr

For an OSL decay curve of type "LM", the fitting model (accroding to Bulur E, 2000) is:\cr\cr
I(t)=a1*(t/P)*exp(-b1*t^2/(2*P))+a2*(t/P)*exp(-b2*t^2/(2*P))+...+\cr
ak*(t/P)*exp(-bk*t^2/(2*P)) , k=1, 2, ..., 7, which is identical to:\cr\cr
I(t)=a1*b1*(t/P)*exp(-b1*t^2/(2*P))+a2*b2*(t/P)*exp(-b2*t^2/(2*P))+...+\cr
ak*bk*(t/P)*exp(-bk*t^2/(2*P)) , k=1, 2, ..., 7 if "transf" is set to be TRUE\cr\cr
where I(t) is the luminescence intensity as a function of time(t), P is the total stimulation time, a is the number of trapped electrons, b is the detrapping rates. Initial detrapping rates are set equal to that listed in Jain et al. (2003). Possible combinations of detrapping rates are calculated with permutation and combinations. For example, if ncomp=3, then the possible number of combinations is C(3,7)=35. Then for each combination that contains three detrapping rates, the number of trapped electrons can be coarsely obtained using a Linear Algebra method (Bluszcz, 1996). The fitting is performed with the Levenberg-Marquardt algorithm.\cr\cr
}
\value{
Return an invisible list containing following elements:
  \item{Comp.Signal}{the fitted decaying signals for each components}
  \item{pars}{the estimates and associated standard errors}
  \item{value}{the sum of squared residuals}
  \item{errorflag}{an error message generated during the calculation, if parameters' standard errors are available, errorflag=0, or else errorflag=1}
}
\seealso{
  \link{decompc}; \link{fastED}
}
\references{
Bluszcz, A., 1996. Exponential function fitting to TL growth data and similar applications. Geochronometria, 13, pp.135-141.

Bluszcz, A., Adamiec, G., 2006. Application of differential evolution to fitting OSL decay curves. Radiation Measurements, 41 (7-8), pp.886-891.

Bulur, E., 2000. A simple transformation for converting CW-OSL curves to LM-OSL curves. Radiation Measurements, 32, pp.141-145.

Differential evolution, http://en.wikipedia.org/wiki/Differential_evolution
  
Jain, M., Murray, A.S., Boetter-Jensen, L., 2003. Characterisation of blue-light stimulated luminescence components in different quartz samples: implications for dose measurement. Radiation Measurements, 37 (4-5), pp.441-449.

Jorge, M., Burton, G., Kenneth, H., 1980. User Guide for MINPACK-1, Technical Report ANL-80-74, Argonne National Laboratory. 

\bold{Further reading} 

Adamiec, G., 2005. OSL decay curves-relationship between single- and multiple-grain aliquots. Radiation Measurements, 39, pp.63-75.

Choi, J.H., Duller G.A.T., Wintle, A.G., 2006. Analysis of quartz LM-OSL curves. Ancient TL, 24, pp. 9-20.

Li, S.H., Li, B., 2006. Dose measurement using the fast component of LM-OSL signals from quartz. Radiation Measurements, 41, pp.534-541.
}

\note{Arugments "LEDpower" and "LEDwavelength" are used for determining the absolute photoionisation cross section. But the calculated absolute photoionisation cross section (pcs) may not be so reliable as it seems, for both the stimulation power and the wavelength vary with the age of the instrument. For fitting an CW-OSL decay curve, in most cases, parameters can be initialized by the differential evolution method, but sometimes the Levenberg-Marquardt algorithm may fails(that's exactly the truth as the increasing of the number of components, which may indicate overcharging of the model). If this happened, parameters initialized with a differential evolution method may be regarded as the final estimation. \cr\cr
I thank Andrzej Bluszcz, who helps me a lot during the programming of this function. Jeong-Heon Choi (Choi et al, 2006) and Bo Li (Li and Li, 2006) are thanked for providing published data to test this routine.
}

\examples{
 data(Signaldata)
 res<-decomp(Signaldata$cw[,c(1,2)],ncomp=3,typ="cw",samplename="GL1")
 print(res)
}
\keyword{fitting OSL decay curve}
