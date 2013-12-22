\name{numOSL-package}
\alias{numOSL-package}
\docType{package}
\title{Package for tackling basic numeric problems in optically stimulated luminescence dating}
\description{
Seven functions are available in package "numOSL" currently.\cr\cr
Function calED() fits dose-response curve and calculates equivalent dose in optically stimulated luminescence dating, assiociated standard error of equivlalent dose can also be estimated using either simple transformation or monte Carlo method (Duller, 2007).\cr\cr
Function dbED() summarizes the distribution of equivalent dose values, drawing histogram and probability density plot. \cr\cr
Function decomp() and decompc() decomposes the OSL decay curve (type "CW" or "LM") to a specified number of components, which supports the combination of differential evolution and
Levenberg-Marquardt algorithm suggested by Bluszcz and Adamiec (2006).\cr\cr
Function fastED() approximates a fast-component equivalent dose with data abtained by the single aliquot regenerative (SAR) protocol (Murray and Wintle, 2000).\cr\cr
Function RadialPlotter() performs statistical age models (log-scale) analysis (Galbraith and Roberts, 2012) and radial plot drawing (Galbraith, 1988). Routines for statistical age models analysis include: \cr
Central age model (CAM) \cr
Finite mixture age model (FMM) \cr
Minimum age model with three parameters (MAM3) \cr
Minimum age model with four parameters (MAM4) \cr\cr
Function sgcED() is a wrapped version of function calED(), fitting OSL dose-response curve and calculating a series of equivalent doses using standardised growth curves (SGC) method (Roberts and Duller, 2004).
}
\details{
\tabular{ll}{
Package: \tab numOSL\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-22\cr
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
\keyword{OSL dating}
\keyword{statistical age models}
\keyword{decay curve decomposition}
\keyword{fast-component De}
\keyword{SGC method}
