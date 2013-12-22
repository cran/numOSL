\name{EDdata}
\alias{EDdata}
\docType{data}
\title{
EDdata is used for exemplifying the usage of function RadialPlotter()
}
\description{
Two sets equivalent dose values.
}
\usage{data(EDdata)}
\format{
A list contains two sets of equivalent dose values
\describe{
  \item{gl11}{a data.frame contains 35 equivalent dose values from aeolain at Tengger Desert (Peng and Han, 2013)}
  \item{al3}{a data.frame contains 84 equivalent dose values from alluvial deposits at andean precordillera (Schmidt et al, 2012)}
}
}
\references{
Schmidt, S., Tsukamoto, S., Salomon, E., Frechen, M., Hetzel, R., 2012. Optical dating of alluvial deposits at the orogenic front of the andean precordillera (Mendoza, Argentina). Geochronometria, 39 (1), pp. 62-75.

Peng, J., Han, F.Q., 2013. Selections of fast-component OSL signal using sediments from the south edge of Tengger Desert. Acta Geoscientica Sinica, 34 (6), pp. 757-762.
}
\examples{
# loading the equivalent dose data
  data(EDdata)
  names(EDdata)
}
\keyword{equivalent dose}
