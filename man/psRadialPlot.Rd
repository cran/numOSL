\name{psRadialPlot}
\alias{psRadialPlot}
\alias{psRadialPlot.default}
\title{
Pseudo radial plot drawing
}
\description{
Drawing a pseudo (simplified) radial plot.
}
\usage{
psRadialPlot(EDdata, addsigma = 0, dose = NULL, zmin = NULL, 
             zmax = NULL, ntick = 6, digits = 2, pcolor = "blue", 
             psize = 1, rg = 2, zlabel = "De (Gy)")
}

\arguments{
  \item{EDdata}{\link{matrix}(\bold{required}): a two-column matrix, equivalent dose values and associated errors are stored in the first and second column, respectively}

  \item{addsigma}{\link{double}(with default): additional uncertainty}

  \item{dose}{\link{vector}(optional): dose population(s) to be drawn}

  \item{zmin}{\link{double}(with default): lower limit on z-axis}

  \item{zmax}{\link{double}(with default): upper limit on z-axis}

  \item{ntick}{\link{integer}(with default): desired number of ticks in z-axis}

  \item{digits}{\link{integer}(with default): number of decimal places or significant digits for values shown in z-axis}

  \item{pcolor}{\link{character}(with default): color of a data point, input \link{colors}() to see more available colors}

  \item{psize}{\link{double}(with default): size of a data point}
  
  \item{rg}{\link{integer}(with default): range of a dose population to be drawn, \code{0=dose},\cr 
\code{1=dose+/-sigma}, \code{2=dose+/-2sigma}}

  \item{zlabel}{\link{character}(with default): title for the z-axis}
}

\details{
Function \emph{psRadialPlot} is used for drawing a simplified radial plot in which the z-axis is a straight line. The pseudo radial plot is easier to construct compared to the regular radial plot. This function can be adopted to display estimates that have different error estimates in any field of the analytical sciences. Note that the function handles a data set in log-scale, hence any minus observation is not allowed.
}

\value{
Return a pseudo radial plot
}
\references{
Galbraith RF, 1988. Graphical display of estimates having differing standard errors. Technometrics, 30(3): 271-281.

Galbraith RF, 1994. Some applications of radial plots. Journal of the American Statistical Association, 89(428): 1232-1242.

Galbraith RF, 2010. On plotting OSL equivalent doses. Ancient TL, 28(1): 1-10.

Galbraith RF, Roberts RG, Laslett GM, Yoshida H, Olley JM, 1999. Optical dating of single grains of quartz from Jinmium rock shelter, northern Australia. Part I: experimental design and statistical models. Archaeometry, 41(2): 339-364.
}

\seealso{
  \link{dbED}; \link{RadialPlotter}
}

\examples{
   data(EDdata)
   psRadialPlot(EDdata$al3, addsigma=0.10, 
                dose=c(39.14, 51.27, 79.14), digits=1,
                zmin=30, zmax=100, ntick=10, rg=1)
}
\keyword{equivalent dose}
\keyword{radial plot}
