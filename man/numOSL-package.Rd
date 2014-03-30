\name{numOSL-package}
\alias{numOSL-package}
\docType{package}
\title{Package for tackling basic numeric problems in optically stimulated luminescence dating}
\description{
Available functions in package "numOSL" include:\cr\cr
\bold{1)}Function calED() fits a dose-response curve and calculates an equivalent dose, the standard error of which can be assessed using either simple transformation or monte Carlo method (Duller, 2007).\cr
\bold{2)}Function dbED() summarizes the distribution of equivalent dose values and draws a histogram and a probability density plot. \cr
\bold{3)}Function decomp() and decompc() decompose an OSL decay curve (type "CW" or "LM") to a specified number of components, which perform the fitting using a combination of differential evolution and
Levenberg-Marquardt algorithm suggested by Bluszcz and Adamiec (2006).\cr
\bold{4)}Function fastED() estimates a fast-component equivalent dose with data abtained by the single aliquot regenerative (SAR) protocol (Murray and Wintle, 2000).\cr
\bold{5)}Function RadialPlotter() performs statistical age models (log-scale) analysis (Galbraith and Roberts, 2012) using classic numeric methods and draw a radial plot(Galbraith, 1988).\cr
\bold{6)}Function sgcED() fits a dose-response curve and calculates a number of equivalent doses using standardised growth curves (SGC) method (Roberts and Duller, 2004).\cr
\bold{7)}Function mcFMM() and mcMAM() estimate parameters of statistical age models with a Markov chain Monte Carlo method (i.e. Slice Sampling, Neal, 2003).
}
\details{
\tabular{ll}{
Package: \tab numOSL\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2014-03-31\cr
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
Peng, J., Dong, Z.B., Han, F.Q., Long, H., Liu, X.J., 2013. R package numOSL: numeric routines for optically stimulated luminescence dating. Ancient TL, 31(2), pp. 41-48.
}
\keyword{OSL dating}
\keyword{statistical age models}
\keyword{fitting OSL decay curve}
\keyword{fast-component equivalent dose}
\keyword{equivalent dose}
