\name{numOSL-package}
\alias{numOSL-package}
\docType{package}
\title{Package for tackling basic numeric problems in optically stimulated luminescence dating}
\description{
Eight functions are available in package "numOSL" currently.\cr\cr
\bold{1)}Function calED() fits dose-response curve and calculates equivalent dose in optically stimulated luminescence dating, assiociated standard error of equivlalent dose can also be estimated using either simple transformation or monte Carlo method (Duller, 2007).\cr
\bold{2)}Function dbED() summarizes the distribution of equivalent dose values, drawing histogram and probability density plot. \cr
\bold{3)}Function decomp() and decompc() decomposes the OSL decay curve (type "CW" or "LM") to a specified number of components, which supports the combination of differential evolution and
Levenberg-Marquardt algorithm suggested by Bluszcz and Adamiec (2006).\cr
\bold{4)}Function fastED() approximates a fast-component equivalent dose with data abtained by the single aliquot regenerative (SAR) protocol (Murray and Wintle, 2000).\cr
\bold{5)}Function RadialPlotter() performs statistical age models (log-scale) analysis (Galbraith and Roberts, 2012) and radial plot drawing (Galbraith, 1988). Routines for statistical age models analysis include: \cr
a)Central age model (CAM) \cr
b)Finite mixture age model (FMM) \cr
c)Minimum age model with three parameters (MAM3) \cr
d)Minimum age model with four parameters (MAM4) \cr
\bold{6)}Function sgcED() is a wrapped version of function calED(), fitting OSL dose-response curve and calculating a series of equivalent doses using standardised growth curves (SGC) method (Roberts and Duller, 2004).
}
\details{
\tabular{ll}{
Package: \tab numOSL\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2014-01-20\cr
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
\keyword{decay curve decomposition}
\keyword{fast-component De}
\keyword{SGC method}
