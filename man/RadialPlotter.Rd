\name{RadialPlotter}
\alias{RadialPlotter}
\alias{RadialPlotter.default}
\title{
 Estimating parameters of Rex Galbraith's statistical age model in OSL dating and drawing radial plot  
}
\description{
Depending on the specified number of components, this function performs statistical age models analysis reviewed in Galbraith and Roberts (2012) dynamically. Statistical OSL age models
that can be performed includes: central age model (CAM), minimum age model (MAM), finite mixture age model (FMM). 
}
\usage{
RadialPlotter(EDdata, 
              ncomp = 0, addsigma = 0, 
              maxiter = 500, maxcomp = 9, 
              algorithm = c("lbfgsb","port"),
              eps = .Machine$double.eps^0.5, 
              plot = TRUE, pcolor = "blue", psize = 1.5, 
              kratio = 0.3, zscale = NULL, samplename = NULL)
}
\arguments{
  \item{EDdata}{\link{data.frame}(\bold{required}): equivalent doses and associated errors (a total of two columns), the routine do analyzes using log-scale, hence all equivalent doses used for fitting must be larger than 0
}
  \item{ncomp}{\link{numeric}(with default): the number of components used for optimizing, 0 means analyzing automatically with finite mxiture age model, 1 for fitting central age model, any other positive integer will fitting finite mixture age model with the specified number of components. As to minimum age models, -1 means performing minimum age model analysis with three parameters, -2 means fitting minimum age model with four parameters
}
  \item{addsigma}{\link{numeric}(with default): the added spread to the relative error of equivalent dose (log-scale), usually a number ranging from 0 to 0.2 is feasible
}
  \item{maxiter}{\link{numeric}(with default): the maximum iterative number allowed in optimization. For central age model and finite mixture age model, often a maximum iterative number of 300 to 500 will be enough in most cases, while in some special situations, 1000 or more will be needed. As to minimum age model, a maximum iterative number of 50 to 100 will be totally enough
}
  \item{maxcomp}{\link{numeric}(with default): the maximum number of components allowed for fitting finite mixture age model, the up limitation
is set to be equal to the length of EDdata, but it will be both time consuming and make no sense if it is set to be too large. It does not exert any effect if minimum age model or central age model is analyzed
}
  \item{algorithm}{\link{character}(with default): the algorithm for optimizing minimum age models, available algorithms are "lbfgsb" and "port", default algorithms="lbfgsb"
}
  \item{eps}{\link{numeric}(with default): the maximum tolerance allowed for stopping the iterative process when performing central age model or finite mixture age model analysis, iterative process will stop if updated parameters satisfying ||pars-newpars||<=eps. It does not exert any effect if minimum age model is analyzed
}
  \item{plot}{\link{logical}(with default): whether drawing the radial plot using estimated parameters or not
}
  \item{pcolor}{\link{character}(with default): color of points for plotting, input \link{colors}() see available colors
}
  \item{psize}{\link{numeric}(with default): size of points for plotting
}
  \item{kratio}{\link{numeric}(with default): an option used for modifying zscale, a bit like a curvature for z-axis
}
  \item{zscale}{\link{numeric}(optional): option for modifying the z-scale of radial plot manually. Example: \code{zscale = seq(min(EDdata),max(EDdata), by = 3L)}
}
  \item{samplename}{\link{character}(optional): the name of the sample, such as "Tengger"
}
}
\details{
Specifying a number of components that not smaller than 2, this function will fit finite mixture age model with the specified number of components. If ncomp is set to be 0, the routine will do the optimization automatically by trying various initial values to pick out a number of components that given the minimum BIC value. Central age model will be analyzed if the number of components is specified to be 1. Minimum age model can also be fitted by specifying ncomp to be -1 (type 3) or -2 (type 4). Both central age model and finite mixture age model are fitted using the maximum likelihood method outlined by Galbraith (1988), while minimum age model can be estimated using the L-BFGS-B subroutine written by Zhu et al (1994) or the "port" rountines (R function \link{nlminb} in package "stats"). Notes that all these calculations are performed using log-scale, which means that any minus equivalent dose will not be allowed to be analyzed by this routine.\cr
\cr
Because both maximum likelihood method and L-BFGS-B algorithm (so is the "port" routine) find out optimized parameters that correspond to a local minimum value instead of a globle minimum value, hence to achieve an acceptable estimation, all routines (except central age model) try various initial values then picking out optimized parameters that given the minimum minus logged likelihood value. When fitting minimum age models (of type 3 or type 4), the lower and upper boundaries are set as following:\cr
\cr
lower boundary for type 3: (p=1e-4, gama=min(ED), sigma=1e-3);\cr
upper boundary for type 3: (p=0.99, gama=max(ED), sigma=5).\cr
lower boundary for type 4: (p=1e-4, gama=min(ED), mu=min(ED), sigma=1e-3);\cr
upper boundary for type 4: (p=0.99, gama=max(ED), mu=max(ED), sigma=5).\cr
\cr 
Also note that, for both finite mixture age model and minimum age model, if parameters' standard error can not be estimated, the result will be suppressed. For minimum age
model with four parameters, if the optimization gives a gama value which is larger than the mu value, the result will be suppressed too.
}
\value{
Return a radial plot, and an invisible list of S3 class "RadialPlotter" that contains following elements:

  \item{errorflag}{error message generated during the fitting, 0 indicates a successful work for analyzed statistical age models. For finite mixture age model, 1 implicates that at least one parameter's standard error can not be approximated. For minimum age model, 1 indicates at least one estimated parameter is near the boundary, or at least one parameter's standard error can not be approximated, or estimated gama value is larger than mu value (for type 4), which might attribute to using of an inappropriate addsigma value or an inappropriate number of parameters (type 3 or type 4) for the sample of interest. As to central age model, a successful work will not be a problem 
}
  \item{commonED}{the equivalent dose and corresponded standard error calculated using common age model}
  \item{ncomp}{the number of components used for optimizing}
  \item{maxcomp}{the maximum number of components allowed for fitting finite mixture age model}
  \item{loopcomp}{indicator of minimum BIC search for finite mixture age model, if ncomp=0 then loopcomp will be 1, else it will be 0}
  \item{pars}{the estimated parameters and standard errors for related OSL age models}
  \item{BIC}{the calculated Bayesian Information Criterion (BIC) value}
  \item{maxlik}{the estimated maximum logged likelihood value}
  \item{BILI}{BIC and maximum logged likelihood values for different number of components in a finite mixture age model, only exist if ncomp=0, else it will be NULL}
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

Schmidt, S., Tsukamoto, S., Salomon, E., Frechen, M., Hetzel, R., 2012. Optical dating of alluvial deposits at the orogenic front of the andean precordillera (Mendoza, Argentina). Geochronometria, 39 (1), pp. 62-75.

Sebastian, K., Christoph, S., Margret, C., F., Michael, D., Manfred, F., Markus, F., 2012. Introducing an R package for luminescence dating analysis. Ancient TL, 30 (1), pp. 1-8.

Vermeesch, P., 2009. RadialPlotter: a Java application for fission track, luminescence and other radial plots, Radiation Measurements, 
44 (4), pp. 409-410.
}

\note{The original S code for drawing radial plot is written by Rex Galbraith and is transformed to R code by Sebastian Kreutzer. The code for drawing radial plot in this package is modified from package "Luminescence" written by Sebastian Kreutzer et al (2012). I thank Rex Galbraith for his permission to modify and bundle the code to this function. I also thank Silke Schmidt for providing published OSL data (Silke Schmidt et al, 2012) for checking and testing this function.
}

\seealso{
  \link{dbED}; \link{print.RadialPlotter}; \link{optim}; \link{nlminb}
}

\examples{
# Loading equivalent dose data
  data(EDdata)
# Specifying the number of component to be 1 to do central age model analysis
  obj<-RadialPlotter(EDdata$gl11,ncomp=1,zscale=seq(20,37,3))
  obj
# Fitting finite mxiture age model with  a number of components of 7 (k=7)
  obj<-RadialPlotter(EDdata$al3,ncomp=7,zscale=seq(25,90,5),pcolor="gray")
  print(obj)
# Picking out the appropriate number of components then fitting automatically
  obj<-RadialPlotter(EDdata$al3,zscale=seq(24,93,7),samplename="AL3",pcolor="brown")
  obj
# Setting zscale manually
  zscale<-seq(min(EDdata$al3[,1])-1,max(EDdata$al3[,1])+1,by = 15L)
  obj<-RadialPlotter(EDdata$al3,maxcomp=13,zscale=zscale,kratio=0.5)
  print(obj)
# Performing minimum age model analysis with three parameters
  obj<-RadialPlotter(EDdata$gl11,ncomp=-1,maxiter=100,zscale=seq(20,37,3))
  class(obj)
# Performing minimum age model analysis with three parameters (addsigma=0.1)
  print(RadialPlotter(EDdata$al3,ncomp=-1,addsigma=0.1,kratio=0.5,
        maxiter=100,algorithm="port",zscale=seq(24,93,7)))
# Performing minimum age model analysis with four parameters
  obj<-RadialPlotter(EDdata$al3,ncomp=-2,maxiter=100,zscale=seq(24,93,7))
  unclass(obj)
}
\keyword{statistical age models}
\keyword{Central Age Model}
\keyword{Finite Mixture Age Model}
\keyword{Minimum Age Model}
