\name{numOSL-package}
\alias{numOSL-package}
\docType{package}
\title{Package for tackling basic numeric problems in optically stimulated luminescence dating}
\description{
Available functions in package \bold{\emph{numOSL}} include:\cr\cr
\bold{Function} \link{analyst} calculates an equivalent dose value using data from a SAR sequence. \cr\cr
\bold{Function} \link{calED} calculates equivalent dose values and assesses their standard errors (Duller, 2007).\cr\cr
\bold{Function} \link{dbED} summarizes the distribution of equivalent dose values. \cr\cr
\bold{Function} \link{decomp} decomposes an OSL decay curve (type "CW" or "LM") to a specified number of components, which perform the fitting using a combination of differential evolution and Levenberg-Marquardt algorithm suggested by Bluszcz and Adamiec (2006).\cr\cr
\bold{Function} \link{fastED} estimates a fast-component equivalent dose using data from a SAR sequence.\cr\cr
\bold{Function} \link{fitGrowth} fits a growth curve.\cr\cr
\bold{Function} \link{mcFMM} and \link{mcMAM} estimate parameters of statistical age models with a Markov chain Monte Carlo method.\cr\cr
\bold{Function} \link{RadialPlotter} performs statistical age models (log-scale) analysis (Galbraith and Roberts, 2012) using the Maximum Likelihood Estimation method.
}
\details{
\tabular{ll}{
Package: \tab numOSL\cr
Type: \tab Package\cr
Version: \tab 1.7\cr
Date: \tab 2014-10-04\cr
License: \tab GPL-3\cr
}
}
\author{
\tabular{ll}{
Peng Jun \tab University of Chinese Academy of Sciences, Beijing, China
}
\bold{Package maintainer}

Peng Jun
\email{pengjun10@mails.ucas.ac.cn}
}
\references{
Peng J, Dong ZB., Han FQ, Long H, Liu XJ, 2013. R package numOSL: numeric routines for optically stimulated luminescence dating. Ancient TL, 31: 41-48.
}
\keyword{OSL dating}
\keyword{statistical age models}
\keyword{OSL decay curve}
\keyword{fast-component equivalent dose}
\keyword{equivalent dose}
