\name{Signaldata}
\alias{Signaldata}
\docType{data}
\title{
Decay curves datasets
}
\description{
CW-OSL and LM-OSL decay curves.
}
\usage{data(Signaldata)}
\format{
  A list that contains CW-OSL and LM-OSL decay curves:
  \describe{
    \item{cw}{a number of CW-OSL decay curves of a sand sample from the Tengger Desert in northern china (Peng and Han, 2013)}
    \item{lm}{a LM-OSL decay curve from Li and Li (2006)}
}
}
\references{
Li SH, Li B, 2006. Dose measurement using the fast component of LM-OSL signals from quartz. Radiation Measurements, 41(5): 534-541.

Peng J, Han FQ, 2013. Selections of fast-component OSL signal using sediments from the south edge of Tengger Desert. Acta Geoscientica Sinica, 34(6): 757-762.
}

\seealso{
  \link{decomp}; \link{fastED}
}

\examples{
# Not run.
# data(Signaldata)
# names(Signaldata)
}
\keyword{OSL dating}
\keyword{decay curve}
