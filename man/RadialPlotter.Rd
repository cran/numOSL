\name{RadialPlotter}
\alias{RadialPlotter}
\alias{RadialPlotter.default}
\title{
 Estimating parameters of Rex Galbraith's statistical age model (with classic numeric methods) and drawing radial plot  
}
\description{
Depending on the specified number of components, this function performs statistical age models analysis reviewed in Galbraith and Roberts (2012) dynamically. Age models
that can be performed includes: central age model (CAM), minimum age model (MAM), finite mixture age model (FMM). 
}
\usage{
RadialPlotter(EDdata, ncomp = 0, addsigma = 0, 
              maxiter = 500, maxcomp = 6, 
              algorithm = c("lbfgsb","port"),
              eps = .Machine$double.eps^0.5, 
              plot = TRUE, pcolor = "blue", psize = 1.5, 
              kratio = 0.3, zscale = NULL, samplename = NULL)
}
\arguments{
  \item{EDdata}{\link{data.frame}(\bold{required}): equivalent dose values and their standard errors (two columns). The routine do the fitting in a log-scale, so all equivalent dose values should larger than zero}
  \item{ncomp}{\link{numeric}(with default): the number of components, 0 for analyzing finite mxiture age models automatically, 1 for fitting a central age model. As to minimum age models, -1 means fitting a 3-parameter minimum age model while -2 means fitting a 4-parameter minimum age model}
  \item{addsigma}{\link{numeric}(with default): the added spread to the relative errors of equivalent dose values, i.e., an additional uncertainty other than those result from counting statistics, instrument reproducibility and curve-fitting}
  \item{maxiter}{\link{numeric}(with default): the allowed maximum numer of iterations}
  \item{maxcomp}{\link{numeric}(with default): the allowed maximum number of components (up to 9) for fitting finite mixture age models}
  \item{algorithm}{\link{character}(with default): an algorithm used for optimizing minimum age models, available algorithms are "lbfgsb" and "port", default algorithms="lbfgsb"}
  \item{eps}{\link{numeric}(with default): the allowed maximum tolerance for stopping the iterative process in fitting of a central age model or a finite mixture age model}
  \item{plot}{\link{logical}(with default): draw a radial plot with the estimates or not}
  \item{pcolor}{\link{character}(with default): the color of points in the plot, input \link{colors}() see available colors}
  \item{psize}{\link{numeric}(with default): the size of points in the plot}
  \item{kratio}{\link{numeric}(with default): an option used for modifying the shape of zscale}
  \item{zscale}{\link{numeric}(optional): an option used for modifying the z-scale of a radial plot. Example: \code{zscale = seq(min(EDdata),max(EDdata), by = 3L)}}
  \item{samplename}{\link{character}(optional): the name of the sample}
}
\details{
Both the central age model and the finite mixture age model are fitted using the maximum likelihood method outlined by Galbraith (1988), while minimum age model can be estimated using either the "L-BFGS-B" subroutine written by Zhu et al (1994) or the "port" rountines (R function \link{nlminb} in package "stats"). Notes that all these calculations are performed in a log-scale, so any minus equivalent dose is not allowed to be analyzed.\cr
\cr
Because both the maximum likelihood method and the L-BFGS-B algorithm (so is the "port" routine) result in a local optimization instead of a globle one, to obtain a feasible estimation, all models (except the central age model) try various initial values during the fitting. For the fitting of minimum age models (of type 3 or type 4), the lower and upper limits are set as following:\cr
\cr
lower boundary for type 3: (p=1e-4, gama=min(ED), sigma=1e-3);\cr
upper boundary for type 3: (p=0.99, gama=max(ED), sigma=5).\cr
lower boundary for type 4: (p=1e-4, gama=min(ED), mu=min(ED), sigma=1e-3);\cr
upper boundary for type 4: (p=0.99, gama=max(ED), mu=max(ED), sigma=5).\cr
\cr 
Also note that, if parameters' standard errors in a finite mixture age model or a minimum age model cannot be estimated, the estimates will be warned (no characteristic dose values will be drawn). For a 4-parameter minimum age model, if the estimate gives a gamma value that is larger than the mu value, the result will be suppressed.
}
\value{
Return a radial plot, and an invisible list of class "RadialPlotter" that contains following elements:

  \item{errorflag}{error message generated during the fitting, 0 indicates a successful estimate. For a finite mixture age model, 1 implicates that at least one parameter's standard error can not be approximated. For a minimum age model, 1 indicates that at least one parameter is near the specified boundary, or at least one parameter's standard error cannot be approximated, or estimated gama value is larger than mu value. As to a central age model, a successful work is not a problem 
}
  \item{commonED}{an equivalent dose calculated using the common age model}
  \item{ncomp}{the used number of components in a finite mixture age model}
  \item{maxcomp}{the allowed maximum number of components for fitting a finite mixture age model}
  \item{loopcomp}{fit a finite mixture age model automatically or not, 1 for "yes" and 0 for "no"}
  \item{pars}{the estimated parameters and standard errors for a given age model}
  \item{BIC}{the calculated Bayesian Information Criterion (BIC) value}
  \item{maxlik}{the estimated logged maximum likelihood value}
  \item{BILI}{BIC values and logged maximum likelihood values for a varying number of components in finite mixture age models (if ncomp=0)}
}
\references{
Galbraith, R.F., 1988. Graphical Display of Estimates Having Differing Standard Errors. Technometrics, 30 (3), pp. 271-281.

Galbraith, R.F., 1990. The radial plot: Graphical assessment of spread in ages. International Journal of Radiation Applications and Instrumentation. Part D. Nuclear Tracks and Radiation Measurements, 17 (3), pp. 207-214.

Galbraith, R.F., Green, P., 1990. Estimating the component ages in a finite mixture. International Journal of Radiation Applications and Instrumentation. Part D. Nuclear Tracks and Radiation Measurements, 17 (3), pp. 197-206.

Galbraith, R.F., Laslett, G.M., 1993. Statistical models for mixed fission track ages. Nuclear Tracks And Radiation Measurements, 21 (4), pp. 459-470.

Galbraith, R.F., 1994. Some Applications of Radial Plots. Journal of the American Statistical Association, 89 (428), pp. 1232-1242.

Galbraith, R.F., Roberts, R.G., Laslett, G.M., Yoshida, H. & Olley, J.M., 1999. Optical dating of single grains of quartz from Jinmium rock shelter, northern Australia. Part I: experimental design and statistical models. Archaeometry, 41 (2), pp. 339-364.

Galbraith, R.F., 2005. Statistics for Fission Track Analysis, Chapman & Hall/CRC, Boca Raton.

Galbraith, R.F., 2010. On plotting OSL equivalent doses. Ancient TL, 28 (1), pp. 1-10.

Galbraith, R.F., Roberts, R.G., 2012. Statistical aspects of equivalent dose and error calculation and display in OSL dating: An overview and some recommendations. Quaternary Geochronology, 11, pp. 1-27.

Zhu, C., Byrd, R.H., Lu, P., Nocedal, J., 1994. 'L-BFGS-B: FORTRAN Subroutines for Large Scale Bound Constrained Optimization' Tech. Report, NAM-11, EECS Department, Northwestern University.

\bold{Further reading}

Bailey, R.M., Arnold, L.J., 2006. Statistical modelling of single grain quartz De distributions and an assessment of procedures for estimating burial dose. Quaternary Science Reviews, 25 (19-20), pp. 2475-2502.

Duller, G.A.T., 2008. Single-grain optical dating of Quaternary sediments: why aliquot size matters in luminescence dating. Boreas, 37 (4), pp. 589-612.

Rodnight, H., 2008. How many equivalent dose values are needed to obtain a reproducible distribution? Ancient TL, 26 (1), pp. 3-10.

Rodnight, H., Duller, G.A.T., Wintle, A.G., Tooth, S., 2006. Assessing the reproducibility and accuracy of optical dating of fluvial deposits. Quaternary Geochronology, 1, pp. 109-120.

Schmidt, S., Tsukamoto, S., Salomon, E., Frechen, M., Hetzel, R., 2012. Optical dating of alluvial deposits at the orogenic front of the andean precordillera (Mendoza, Argentina). Geochronometria, 39 (1), pp. 62-75.

Sebastian, K., Christoph, S., Margret, C., F., Michael, D., Manfred, F., Markus, F., 2012. Introducing an R package for luminescence dating analysis. Ancient TL, 30 (1), pp. 1-8.

Vermeesch, P., 2009. RadialPlotter: a Java application for fission track, luminescence and other radial plots, Radiation Measurements, 
44 (4), pp. 409-410.
}

\note{The original S code for drawing radial plot is written by Rex Galbraith and is transformed to R code by Sebastian Kreutzer. The code for drawing radial plot in this package is modified from package "Luminescence" written by Sebastian Kreutzer et al (2012). I thank Rex Galbraith for his permission to modify and bundle the code to this function. I also thank Silke Schmidt (Silke Schmidt et al, 2012) and Helena Rodnight (Rodnight 2008; Rodnight et al, 2006) for providing published OSL data  for checking and testing this function.
}

\seealso{
  \link{dbED}; \link{print.RadialPlotter}; \link{mcMAM}; \link{mcFMM}; \link{optim}; \link{nlminb}
}

\examples{
# Loading equivalent dose data
  data(EDdata)
# Finding the appropriate number of components of a finite mixture age model
  obj<-RadialPlotter(EDdata$al3,zscale=seq(24,93,7),samplename="AL3")
  print(obj)
# Fitting a 3-parameter minimum age model 
  obj<-RadialPlotter(EDdata$gl11,ncomp=-1,maxiter=100,zscale=seq(20,37,3))
  unclass(obj)
}
\keyword{Radial Plot}
\keyword{equivalent dose}
\keyword{Central Age Model}
\keyword{Finite Mixture Age Model}
\keyword{Minimum Age Model}
