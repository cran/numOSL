\name{BIN}
\alias{BIN}
\docType{data}

\title{
BIN data
}

\description{
BIN data for aeolian sample GL2-1 from the south margin of the Tengger Desert (Peng et al., 2013).
}

\usage{data(BIN)}

\format{
A S3 class object "loadBIN" produced by function \link{loadBINdata} that contains the following two elements:
\describe{
  \item{records}{a \link{list} consists of loaded data records for each aliquot (grain)}
  \item{tab}{a \link{data.frame} used for summarizing loaded data records}
}
}

\references{
Peng J, Han FQ, 2013. Selections of fast-component OSL signal using sediments from the south edge of Tengger Desert. 
Acta Geoscientica Sinica, 34(6): 757-762.
}

\seealso{
  \link{loadBINdata}; \link{pickBINdata}; \link{analyseBINdata}
}

\examples{
 # Not run.
 # data(BIN)
 # class(BIN)
}
\keyword{OSL dating}
\keyword{BIN file}
