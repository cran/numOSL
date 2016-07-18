\name{RadialPlotter}
\alias{RadialPlotter}
\alias{RadialPlotter.default}
\title{
Statistical age model optimization (with a Maximum Likelihood Estimation method) 
}
\description{
Depending on the specified number of components, this function performs statistical age models analysis reviewed in Galbraith and Roberts (2012) dynamically using a Maximum Likelihood Estimation method. Age models that can be applied include: central age model (CAM), minimum age model (MAM), and finite mixture age model (FMM). 
}
\usage{
RadialPlotter(EDdata, ncomp = 0, addsigma = 0, 
              maxcomp = 6, algorithm = c("port","lbfgsb"),
              plot = TRUE, pcolor = "blue", psize = 1.5, 
              kratio = 0.3, zscale = NULL)
}
\arguments{
  \item{EDdata}{\link{matrix}(\bold{required}): a two-column matrix (i.e., equivalent dose values and associated standard errors)}

  \item{ncomp}{\link{integer}(with default): number of components, \code{-1=MAM3}, \code{-2=MAM4}, \code{0} means fitting FMM automatically, \code{1=CAM}, and \code{>=1} means fitting FMM with a given number of components}

  \item{addsigma}{\link{double}(with default): additional uncertainty}

  \item{maxcomp}{\link{integer}(with default): maximum number of components in FMM}

  \item{algorithm}{\link{character}(with default): algorithm used for optimizing MAM,\cr
default \code{algorithm="port"}}

  \item{plot}{\link{logical}(with default): draw a radial plot or not}

  \item{pcolor}{\link{character}(with default): color of a data point, input \link{colors}() to see more available colors}

  \item{psize}{\link{double}(with default): size of a data point}

  \item{kratio}{\link{double}(with default): argument controlling the shape of zscale}

  \item{zscale}{\link{vector}(optional): argument controlling the scale of z-axis.\cr
Example: \code{zscale=seq(min(EDdata),max(EDdata),by=3L)}}
}

\details{ 
Both CAM and FMM are fitted using a iterative Maximum Likelihood Estimation procedure outlined by Galbraith (1988), while MAM can be estimated using either the "L-BFGS-B" algorithm (\bold{\emph{R}} function \link{optim} in package \bold{\emph{stats}}) or the "port" algorithm (\bold{\emph{R}} function \link{nlminb} in package \bold{\emph{stats}}).
}

\value{
Return a list of S3 class \code{"RadialPlotter"} that contains the following elements:
  \item{pars}{optimizaed parameters}
  \item{bic}{calculated Bayesian Information Criterion (BIC) value}
  \item{maxlik}{optimized maximum logged likelihood value}
}

\references{
Galbraith RF, 1988. Graphical display of estimates having differing standard errors. Technometrics, 30(3): 271-281.

Galbraith RF, 1990. The radial plot: Graphical assessment of spread in ages. International Journal of Radiation Applications and Instrumentation. Part D. Nuclear Tracks and Radiation Measurements, 17(3): 207-214.

Galbraith RF, Green P, 1990. Estimating the component ages in a finite mixture. International Journal of Radiation Applications and Instrumentation. Part D. Nuclear Tracks and Radiation Measurements, 17: 197-206.

Galbraith RF, Laslett GM, 1993. Statistical models for mixed fission track ages. Nuclear Tracks And Radiation Measurements, 21(4): 459-470.

Galbraith RF, 1994. Some applications of radial plots. Journal of the American Statistical Association, 89(428): 1232-1242.

Galbraith RF, Roberts RG, Laslett GM, Yoshida H, Olley JM, 1999. Optical dating of single grains of quartz from Jinmium rock shelter, northern Australia. Part I: experimental design and statistical models. Archaeometry, 41(2): 339-364.

Galbraith RF, 2005. Statistics for fission track analysis. Chapman & Hall/CRC Press.

Galbraith RF, 2010. On plotting OSL equivalent doses. Ancient TL, 28(1): 1-10.

Galbraith RF, Roberts RG, 2012. Statistical aspects of equivalent dose and error calculation and display in OSL dating: an overview and some recommendations. Quaternary Geochronology, 11: 1-27.

\bold{Further reading}

Duller GAT, 2008. Single-grain optical dating of Quaternary sediments: why aliquot size matters in luminescence dating. Boreas, 37(4): 589-612.

Kreutzer S, Schmidt C, Fuchs MC, Dietze M, Fischer M, Fuchs M, 2012. Introducing an R package for luminescence dating analysis. Ancient TL, 30(1): 1-8. 
Software is freely available at \url{http://CRAN.R-project.org/package=Luminescence}.

Rodnight H, 2008. How many equivalent dose values are needed to obtain a reproducible distribution? Ancient TL, 26(1): 3-10.

Rodnight H, Duller GAT, Wintle AG, Tooth S, 2006. Assessing the reproducibility and accuracy of optical dating of fluvial deposits. Quaternary Geochronology, 1(2): 109-120.

Schmidt S, Tsukamoto S, Salomon E, Frechen M, Hetzel R, 2012. Optical dating of alluvial deposits at the orogenic front of the andean precordillera (Mendoza, Argentina). Geochronometria, 39(1): 62-75.

Vermeesch P, 2009. RadialPlotter: a Java application for fission track, luminescence and other radial plots. Radiation Measurements, 44: 409-410. 
Software is freely available at \url{http://www.ucl.ac.uk/~ucfbpve/radialplotter/}.
}

\note{
Function \link{RadialPlotter} was given the same name as the \bold{\emph{Java}} package \bold{\emph{RadialPlotter}} written by Pieter Vermeesch (Vermeesch, 2009). Note that this function fits a model in log-scale, hence any minus equivalent dose value is not allowed, and that the procedure will return an error if any standard error of a parameter cannot be estimated by numerical difference-approximation.\cr\cr
The original \bold{\emph{S}} code for drawing a radial plot was written by Rex Galbraith and was transformed to \bold{\emph{R}} by Sebastian Kreutzer. The code for drawing radial plot in this function was modified from package \bold{\emph{Luminescence}} written by Kreutzer et al. (2012). We thank Dr Rex Galbraith for his permission to modify and bundle the code to this function. We also thank Dr Silke Schmidt, Dr Helena Rodnight, Dr Xian-Jiao Ou, and Dr Amanda Keen-Zebert for providing published OSL data sets to test this routine.
}

\seealso{
  \link{dbED}; \link{mcMAM}; \link{mcFMM}; \link{psRadialPlot}
}

\examples{
  data(EDdata)
  # Find out the appropriate number of components 
  # in FMM and fit automatically.
  RadialPlotter(EDdata$al3,zscale=seq(24,93,7))
  # Fit MAM3. 
  RadialPlotter(EDdata$gl11,ncomp=-1,zscale=seq(20,37,3))
}
\keyword{statistical age models}
\keyword{radial plot}
\keyword{equivalent dose}
\keyword{Central Age Model}
\keyword{Finite Mixture Age Model}
\keyword{Minimum Age Model}
