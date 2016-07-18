\name{select}
\alias{select}
\alias{select.default}
\title{
BIN file records selection
}
\description{
Extracting data records according to certain criteria.
}
\usage{
select(obj, position = NULL, runNumber = NULL, nPoints = NULL, 
       low = NULL, high = NULL, irrTime = NULL, lumType = NULL, 
       dataType = NULL, lightSource = NULL)
}

\arguments{
  \item{obj}{\link{list}(\bold{required}): an object of S3 class \code{"binfile"}, which is produced by function \link{viewer}}

  \item{position}{\link{vector}(optional): carousel position(s)}

  \item{runNumber}{\link{vector}(optional): run number(s)}

  \item{nPoints}{\link{vector}(optional): number(s) of points of records}

  \item{low}{\link{double}(optional): lower limit on temperature, time, or wavelength}

  \item{high}{\link{double}(optional): upper limit on temperature, time, or wavelength}

  \item{irrTime}{\link{vector}(optional): irradiation time(s)}

  \item{lumType}{\link{character}(optional): luminescence type(s), see \link{viewer} for available luminescence types}

  \item{dataType}{\link{character}(optional): data type(s), see \link{viewer} for available data types}

  \item{lightSource}{\link{character}(optional): light source(s), see \link{viewer} for available light sources}
}

\details{
Function \link{select} is used for extracting records that satisfy a number of specified criteria (attributes) from an object of S3 class \code{"binfile"}. Extracted equal-length records are stored in a matrix column by column and can be used in further analysis.
}

\value{
Return a list containing the following elements:
\item{records}{a matrix (if records are of equal-length) or a list (if records have different lengths) that contains selected data records}
\item{tab}{a table that summaries the data structure and attributes of selected records}
}

\seealso{
  \link{pickSARdata}; \link{pickSARED}; \link{viewer}
}

\examples{
   # Do not run.
   # obj <- viewer("foo.bin")
   # select(obj, position=2, lumType="osl")
}
\keyword{BIN file}
\keyword{data extraction}
