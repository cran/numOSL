\name{decomp}
\alias{decomp}
\alias{decomp.default}
\title{
OSL decay curve decomposition
}
\description{
Decomposing a CW-OSL or LM-OSL decay curve to a given number of first-order exponential components using a combination of differential evolution and Levenberg-Marquardt algorithm.
}

\usage{
decomp(Sigdata, ncomp = 2, constant = TRUE, typ = "cw", 
       control.args = list(), transf = TRUE, weight = FALSE, 
       plot = TRUE, xylog = FALSE, lwd = 3, outfile = NULL)
}

\arguments{
  \item{Sigdata}{\link{matrix}(\bold{required}): a two-column matrix (i.e., stimulation time values and photon counts)}

  \item{ncomp}{\link{integer}(with default): number of decomposed components}

  \item{constant}{\link{logical}(with default): subtract a constant component or not}

  \item{typ}{\link{character}(with default): type (\code{"cw"} or \code{"lm"}) of a decay curve}

  \item{control.args}{\link{list}(with default): arguments used by the differential evolution algorithm, see details}

  \item{transf}{\link{logical}(with default): do not use (for backward compatibility purpose)}

  \item{weight}{\link{logical}(with default): perform the optimization using a weighted or an un-weighted procedure, see details}

  \item{plot}{\link{logical}(with default): draw a plot or not}

  \item{xylog}{\link{logical}(with default): use a logarithmic scale for both x-axis and y-axis or only for x-axis}

  \item{lwd}{\link{double}(with default): line width}

  \item{outfile}{\link{character}(optional): if specified, decomposed signal values will be written to a CSV file named \code{"outfile"} and saved to the current work directory}
}

\details{
Function \link{decomp} decomposes an OSL decay curve to a specified number of components using a combination of differential evolution and Levenberg-Marquardt algorithm. Both CW-OSL and LM-OSL decay curves can be decomposed.\cr\cr
For a CW-OSL decay curve, the fitting model (Bluszcz and Adamiec, 2006) is:\cr
\code{I(t)=a1*b1*exp(-b1*t)+...+ak*bk*exp(-bk*t)},\cr
where \code{k=1, 2, ..., 7}, \code{I(t)} is the luminescence intensity as a function of time, \code{a} is the number of trapped electrons, and \code{b} is the detrapping rate. The constant component is \code{c} if \code{constant=TRUE}.\cr\cr
For a LM-OSL decay curve, the fitting model (Bulur, 2000) is:\cr
\code{I(t)=a1*b1*(t/P)*exp[-b1*t^2/(2*P)]+...+ak*bk*(t/P)*exp[-bk*t^2/(2*P)]},\cr
where \code{k=1, 2, ..., 7}, and \code{I(t)} is the luminescence intensity as a function of time, \code{P} is the total stimulation time, \code{a} is the number of trapped electrons, and \code{b} is the detrapping rate. The constant component is \code{c*(t/P)} if \code{constant=TRUE}.\cr\cr
Parameters are initialized using a differential evolution method suggested by Bluszcz and Adamiec (2006), then the Levenberg-Marquardt algorithm (minpack: Fortran 90 version by John Burkardt, freely available at \url{http://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.html}) will be performed to optimize the parameters. If \code{weight=TRUE}, the photon counts will be assumed to follow a Possion distribution with a standard error equal to the square root of the number of photon counts, and the decay curve will be fitted using a weighted procedure. Setting \code{weight=TRUE} gives more weight to photon counts from slower decaying components.\cr\cr
Arguments in \code{control.args} that control the differential evolution algorithm include:\cr
\bold{(1)} \emph{factor}: the number of population members, \code{np=factor*ncomp}, default \code{factor=10}; \cr
\bold{(2)} \emph{f}: a weighting factor that lies between \code{0} and \code{1.2}, default \code{f=0.5}; \cr
\bold{(3)} \emph{cr}: a crossover constant that lies between \code{0} and \code{1}, default \code{cr=0.99}; \cr
\bold{(4)} \emph{maxiter}: maximum number of iterations, default \code{maxiter=500}; \cr
\bold{(5)} \emph{tol}: tolerance for stopping the iteration, the procedure will be terminated if \cr
all relative standard deviations of parameters are smaller than \code{tol}, defalut \code{tol=0.1}.
}

\value{
Return a list containing the following elements:
  \item{LMpars}{optimized parameters for the fitted decay curve}
  \item{constant}{estimated constant component (return \code{0} if \code{constant=FALSE)}}
  \item{value}{minimized objective for the fitted decay curve}
}

\seealso{
  \link{fastED}
}

\references{
Bluszcz A, 1996. Exponential function fitting to TL growth data and similar applications. Geochronometria, 13: 135-141.

Bluszcz A, Adamiec G, 2006. Application of differential evolution to fitting OSL decay curves. Radiation Measurements, 41(7-8): 886-891.

Bulur E, 2000. A simple transformation for converting CW-OSL curves to LM-OSL curves. Radiation Measurements, 32(2): 141-145.

Differential evolution algorithm, \url{http://en.wikipedia.org/wiki/Differential_evolution}
  
Jain M, Murray AS, Boetter-Jensen L, 2003. Characterisation of blue-light stimulated luminescence components in different quartz samples: implications for dose measurement. Radiation Measurements, 37(4-5): 441-449.

More JJ, 1978. "The Levenberg-Marquardt algorithm: implementation and theory," in Lecture Notes in Mathematics: Numerical Analysis, 
Springer-Verlag: Berlin. 105-116.

\bold{Further reading} 

Adamiec G, 2005. OSL decay curves-relationship between single- and multiple-grain aliquots. Radiation Measurements, 39(1): 63-75.

Choi JH, Duller GAT, Wintle AG, 2006. Analysis of quartz LM-OSL curves. Ancient TL, 24(1): 9-20.

Li SH, Li B, 2006. Dose measurement using the fast component of LM-OSL signals from quartz. Radiation Measurements, 41(5): 534-541.

Peng J, Dong ZB, Han FQ, Han YH, Dai XL, 2014. Estimating the number of components in an OSL decay curve using the Bayesian Information Criterion. Geochronometria, 41(4): 334-341.
}


\note{
The model to be optimized should not be underdetermined. This means that the number of data points should exceed (or at least be equal to) the number of parameters. For a given model, this routine will return an error if any standard errors of parameters cannot be estimated by numerical difference-approximation. Function \code{decompc} in previous versions was bundled to function \link{decomp}.\cr\cr
We would like to thank Professor Andrzej Bluszcz who helps us a lot during the programming of this function. Dr Jeong-Heon Choi is thanked for providing published data sets to test this routine.
}

\examples{
 data(Signaldata)
 decomp(Signaldata$lm,ncomp=3,typ="lm",
        control.args=list(maxiter=10))
}
\keyword{decay curve}
