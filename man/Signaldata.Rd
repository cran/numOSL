\name{Signaldata}
\alias{Signaldata}
\docType{data}
\title{
Signaldata is used for exemplifying the usage of function decomp(), decompc() and fastED()
}
\description{
Decay curves of type "CW" and "LM".
}
\usage{data(Signaldata)}
\format{
  A list that contains both decay curves of type "CW" and "LM"
  \describe{
    \item{cw}{a number of decay curves of type "CW" of sand samples from Tengger Desert in northern china (Peng and Han, 2013)}
    \item{lm}{LM-OSL decay curve from Li and Li (2006)}
    \item{lm1}{LM-OSL decay curve from Choi et al (2006)}
}
}
\references{
Li, S.H., Li, B., 2006. Dose measurement using the fast component of LM-OSL signals from quartz. Radiation Measurements, 41, pp.534-541.

Choi, J.H., Duller G.A.T., Wintle, A.G., 2006. Analysis of quartz LM-OSL curves. Ancient TL, 24, pp. 9-20.

Peng, J., Han, F.Q., 2013. Selections of fast-component OSL signal using sediments from the south edge of Tengger Desert. Acta Geoscientica Sinica, 34 (6), pp. 757-762.
}
\examples{
data(Signaldata)
names(Signaldata)
}
\keyword{OSL dating}
\keyword{fitting OSL decay curve}
