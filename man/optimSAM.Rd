\name{optimSAM}
\alias{optimSAM}
\alias{optimSAM.default}
\title{
Optimization of statistical age models
}
\description{
Estimating the parameters of statistical age models, including the common age model (COM), the central age model (CAM), the minimum age model (MAM), the maximum age model (MXAM), and the finite mixture age model (FMM), using the Maximum Likelihood Estimation method.
}
\usage{
optimSAM(EDdata, model = "cam", addsigma = 0, iflog = TRUE, maxcomp = 6)
}
\arguments{
  \item{EDdata}{\link{matrix}(\bold{required}): a two-column matrix (i.e., equivalent dose values and\cr
associated standard errors)}

  \item{model}{\link{character}(with default): the fitting model, one of \code{"com"}, \code{"cam"}, \code{"mam3"}, \code{"mam4"},
\code{"mxam3"}, \code{"mxam4"}, \code{"fmm0"}, \code{"fmm1"}, \code{"fmm2"}, ..., \code{"fmm9"}}

  \item{addsigma}{\link{numeric}(with default): additional uncertainty, i.e., the sigmab value}

  \item{iflog}{\link{logical}(with default): transform equivalent dose values to log-scale or not}

  \item{maxcomp}{\link{integer}(with default): the maximum number of components in FMM}
}

\value{
Return a \link{list} that contains the following elements:
  \item{pars}{optimized parameters, the name of the parameter of COM is \code{"COM.De"}, that of CAM are \code{c("CAM.OD","CAM.De")}, that of MAM3 are \code{c("Prop","MAM3.De","Sigma")}, that of MXAM3 are \code{c("Prop","MXAM3.De","Sigma")}, that of MAM4 are\cr \code{c("Prop","MAM4.De","Mu","Sigma")}, that of MXAM4 are\cr \code{c("Prop","MXAM4.De","Mu","Sigma")}, and that of FMM are 
\code{c("Prop","FMM.De")}}
  \item{maxlik}{optimized maximum logged likelihood value}
  \item{bic}{calculated Bayesian Information Criterion (BIC) value}
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
}

\seealso{
 \link{mcMAM}; \link{mcFMM}; \link{dbED}; \link{psRadialPlot}; \link{RadialPlotter}; \link{EDdata}
}
\examples{
  data(EDdata)
  
  ### Fitting a 3-component FMM.
  optimSAM(EDdata$al3, model="fmm3", addsigma=0)

  ### Fitting a 4-parameter MXAM.
  optimSAM(EDdata$al3, model="mxam4", addsigma=0.1)
}
\keyword{statistical age models}
\keyword{equivalent dose}
\keyword{Central Age Model}
\keyword{Finite Mixture Age Model}
\keyword{Minimum Age Model}
\keyword{Maximum Age Model}
