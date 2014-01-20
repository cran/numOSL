\name{decomp}
\alias{decomp}
\alias{decomp.default}
\title{
Decomposition decay curve in OSL dating
}
\description{
Function decomp() decomposes the OSL decay curve (type "CW" or "LM") to a specified number of components.
}
\usage{
decomp(Sigdata, ncomp = 3, typ = c("cw","lm"),
       control.args = list(), transf = FALSE, 
       LEDpower = 60, LEDwavelength = 470, plot = TRUE,
       xlog= TRUE, lwd = 3, samplename = NULL, outfile = NULL)
}

\arguments{
  \item{Sigdata}{\link{data.frame}(\bold{required}): two columns, time (t1,t2,etc) and correspond signal values (s1,s2,etc) to be decomposed
}
  \item{ncomp}{\link{numeric}(with default): the number of components to be decomposed, which must be an integer in the space [1,7]
}
  \item{typ}{\link{character}(with default): type ("cw" or "lm") of the OSL decay curve, default typ="cw"
}
  \item{control.args}{\link{list}(with default): further arguments used in differential evolution algorithm, see details
}
  \item{transf}{\link{logical}(with default): whether transform estimated parameters or not, see details
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
This function decomposes OSL decay curve to a specified number of components using a combination of differential evolution and
Levenberg-Marquardt algorithm. Both OSL decay curve of type "CW" and "LM" can be decomposed.\cr
For OSL decay curve of type "CW", the fitting model (accroding to Bluszcz and Adamiec, 2006) is:\cr\cr
I(t)=a1*exp(-b1*t)+a2*exp(-b2*t)+...+ak*exp(-bk*t) , k=1, 2, ..., 7, which is identical to:\cr\cr
I(t)=a1*b1*exp(-b1*t)+a2*b2*exp(-b2*t)+...+ak*bk*exp(-bk*t) , k=1, 2, ..., 7 if "transf" is set to be TRUE.\cr\cr
where I(t) is the luminescence intensity as a function of time(t), a is the number of trapped electrons, b is the detrapping rates. Initial parameters are estimated using the differential evolution method suggested by Bluszcz A and Adamiec G (2006), then Levenberg-Marquardt algorithm (minpack: original FORTRAN 77 version by Jorge More, Burton Garbow, Kenneth Hillstrom. FORTRAN 90 version by John Burkardt) will be performed to obtain a final optimization. In some radical cases, differential evolution will fail in parameter initializing, or Levenberg-Marquardt algorithm may fails, then a serise of simple attempts will be carryed out using Levenberg-Marquardt algorithm. Arguments (control.args) that employed to control the performance of differential evolution algorithm are:\cr\cr
\bold{factor}: used for control the number of parents (NP), np=factor*ncomp, default factor=20; \cr
\bold{f}: weighting factor, ranging from 0 to 1.2, default f=0.5; \cr
\bold{cr}: crossover constant, ranging from 0 to 1, default cr=0.9; \cr
\bold{maxiter}: the maximum iterative number allowed, default maxiter=1e3; \cr
\bold{tol}: a tolerance for stopping the iteration, if magnitude of all standard deviations of parameters within the whole population are smaller than tol, the iteration will be terminated. default tol=0.1 \cr\cr

For OSL decay curve of type "LM", the fitting model (accroding to Bulur E, 2000) is:\cr\cr
I(t)=a1*(t/P)*exp(-b1*t^2/(2*P))+a2*(t/P)*exp(-b2*t^2/(2*P))+...+\cr
ak*(t/P)*exp(-bk*t^2/(2*P)) , k=1, 2, ..., 7, which is identical to:\cr\cr
I(t)=a1*b1*(t/P)*exp(-b1*t^2/(2*P))+a2*b2*(t/P)*exp(-b2*t^2/(2*P))+...+\cr
ak*bk*(t/P)*exp(-bk*t^2/(2*P)) , k=1, 2, ..., 7 if "transf" is set to be TRUE\cr\cr
where I(t) is the luminescence intensity as a function of time(t), P is the total stimulation time, a is the number of trapped electrons, b is the detrapping rates. Initial detrapping rates are estimated using the detrapping values described after Jain et al. (2003). Possible combinations of detrapping rates for a given component are calculated using permutation and combinations. For example, if ncomp=3, then the possible number of combinations is C(3,7)=35. Then for each combination that contains three detrapping rates, the correspond number of trapped electrons can be coarsely obtained using Linear Algebra method outlined by Bluszcz A (1996). Then the Levenberg-Marquardt algorithm will be used to optimize the finial result.\cr\cr
}
\value{
Return an invisible list containing following elements:\cr
  \item{Comp.Signal}{fitted values for observations and each component}
  \item{pars}{estimated parameters and associated standard errors}
  \item{value}{sum of square of residual errors}
  \item{errorflag}{error message generated during the calculation, if parameters' standard errors can be estimated, errorflag=0, else 1}
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

\note{Both LEDpower and LEDwavelength are used for determining the absolute photoionisation cross section. However, the estimated absolute photoionisation cross section (pcs) may not so reliable as it seems to be, because both LEDpower and LEDwavelength will vary with the age of the instrument. For fitting OSL decay curve of type "CW", in most cases, parameters can be initialized by differential evolution method, but sometimes Levenberg-Marquardt algorithm fails to obtain a final optimization (that's exactly the truth as the increasing of the number of components, which might indicate overcharging of the model), if this happened, parameters initialized with differential evolution method will be regarded as the final estimation. For better performance in CW-OSL curve fitting, trying remodifying default parameters in argument "control.args".\cr\cr
I thank Andrzej Bluszcz, who helps me a lot during the writting of this function. Jeong-Heon Choi (Choi et al, 2006) and Bo Li (Li and Li, 2006) are thanked for providing published data to test this routine.
}

\examples{
 # loading signal data
 data(Signaldata)
 # decompose natural OSL signal of type "CW" to 3 components
 res<-decomp(Signaldata$cw[,c(1,2)],ncomp=3,typ="cw",samplename="GL1")
 res 
 # decompose OSL signal of type "LM" to 4 components
 res<-decomp(Signaldata$lm,ncomp=4,typ="lm",transf=TRUE)
 print(res)
}
\keyword{CW-OSL decay curve}
\keyword{LM-OSL decay curve}
\keyword{decay curve decomposition}
