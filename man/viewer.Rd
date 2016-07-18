\name{viewer}
\alias{viewer}
\alias{viewer.default}
\title{
BIN file importing
}
\description{
Importing a BIN file into R platform to visualize and extract measurement records.
}
\usage{
viewer(filename)
}

\arguments{
  \item{filename}{\link{character}(\bold{required}): name of a file (in ".bin" or "binx" format), the user need ensure that the file is available from the current work directory}
}

\details{
Function \link{viewer} is used for importing data records stored in a BIN file produced by a Risoe reader into the R platform. It visualizes datasets in a manner similar to software \bold{\emph{Viewer}}. Two kinds of binary data file are supported: (1) file in ".bin" format with data structure described in Appendix 1 in user mannual (version 3.22b) of the \bold{\emph{Analyst}} software; (2) file in ".binx" format with data structure outlined in user manual (February 2013) of the \bold{\emph{Viewer}} software. \cr\cr
Each data record stored in a BIN file is associated with a number of items (attributes). Only items that facilate record identification and equivalent dose calculation are reserved. Items of this kind may include:\cr
\bold{(1)} \emph{position}: carousel position of an aliquot; \cr
\bold{(2)} \emph{runNumber}: run number of measurement; \cr
\bold{(3)} \emph{nPoints}: number of data points of the current record; \cr
\bold{(4)} \emph{low}: lower limit on temperature (TL stimulation), time (OSL stimulation), or wavelength; \cr
\bold{(5)} \emph{high}: upper limit on temperature (TL stimulation), time (OSL stimulation), or wavelength; \cr
\bold{(6)} \emph{irrTime}: irradiation time; \cr
\bold{(7)} \emph{lumtype}: luminescence type (includes: \code{tl}, \code{osl}, \code{irsl}, \code{m-ir}, \code{m-vis}, \code{tol},\cr 
\code{pulsed}, \code{trposl}, \code{rir}, \code{rbr}, \code{user}, \code{posl}, \code{sgosl}, \code{rl}, \code{xrf}); \cr
\bold{(8)} \emph{dataType}: data type (includes: \code{natural}, \code{n+dose}, \code{bleach}, \code{bleach+dose},\cr
\code{natural(bleach)}, \code{n+dose(bleach)}, \code{dose}, \code{background});\cr
\bold{(9)} \emph{lightSource}: light source (includes: \code{none}, \code{lamp}, \code{ir-diodes/ir-laster},\cr
\code{calibraition-led}, \code{blue-diodes}, \code{white-light}, \code{green-laser(single-grain)},\cr
\code{ir-laser(single-grain)}). Note that long strings \code{ir-diodes/ir-laser},\cr
\code{green-laser(single-grain)}, \code{ir-laser(single-grain)} are abbreviated as \cr
\code{ir-diodes}, \code{green-laser(sg)}, \code{ir-laser(sg)}, respectively.
}

\value{
Return an invisible list of S3 class \code{"binfile"} containing the following elements:
\item{records}{a list that contains imported data records, each record is a numeric vector of \code{nPoints} long with 9 attributes (see above)}
\item{tab}{a table that summaries the data structure and attributes of imported data records}
}

\references{
Both documents (user manuals of \bold{\emph{Analyst}} and \bold{\emph{Viewer}}) are available in \url{http://www.nutech.dtu.dk/english/Products-and-Services/Dosimetry/Radiation-Measurement-Instruments/TL_OSL_reader/Manuals}.
}

\note{
We would like to appreciate Dr Lei Gao who prompts me to write this function and provides measured data sets to test this procedure.
}
\seealso{
  \link{select}
}

\examples{
   ### please ensure that file "foo.bin" is 
   ### available in the current work directory.
   # obj <- viewer("foo.bin")
   # class(obj)
   # obj$tab
}

\keyword{data importing}
\keyword{BIN file}
