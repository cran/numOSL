\name{loadBINdata}
\alias{loadBINdata}
\alias{loadBINdata.default}
\title{
BIN file loading (importing)
}
\description{
Loading (importing) a BIN file into the R platform.
}

\usage{
loadBINdata(filename, view = TRUE)
}

\arguments{
  \item{filename}{\link{character}(\bold{required}): name(s) of file(s) (with file extension ".BIN", ".bin", "BINX", or "binx"), the file(s) must be available from the current working directory. Example: \code{filename=c("foo1.bin","foo2.binx")}}

  \item{view}{\link{logical}(optional): logical value indicating if the loaded data should be visualized in a Summary Table}
}

\details{
Function \link{loadBINdata} is used for loading BIN (BINX) files into the R platform. Five versions of binary files (V3, V4, V6, V7, and V8) are loadable. It can load a single BIN (BINX) file or a number of files into R simultaneously.\cr
Items reserved during the loading process include:\cr
\bold{(1)} \emph{Position}: Carousel position;

\bold{(2)} \emph{Grain}: Grain number;

\bold{(3)} \emph{Run}: Run number;

\bold{(4)} \emph{Set}: Set number;

\bold{(5)} \emph{DType}: Data type, includes: \code{Natural}, \code{N+dose}, \code{bleach}, \code{Bleach+dose},\cr
  \code{Natural(Bleach)}, \code{N+dose(Bleach)}, \code{Dose}, \code{Background};

\bold{(6)} \emph{IRRTime}: Irradiation time;

\bold{(7)} \emph{NPoints}: number of data points;

\bold{(8)} \emph{LType}: Luminescence type, includes: \code{TL}, \code{OSL}, \code{IRSL}, \code{M-IR}, \code{M-VIS},\cr 
  \code{TOL}, \code{TRPOSL}, \code{RIR}, \code{RBR}, \code{USER}, \code{POSL}, \code{SGOSL}, \code{RL}, \code{XRF};

\bold{(9)} \emph{Low}: Low (temperature, time, wavelength);

\bold{(10)} \emph{High}: High (temperature, time, wavelength);

\bold{(11)} \emph{Rate}: Rate (temperature, time, wavelength);

\bold{(12)} \emph{Temperature}: Sample temperature;

\bold{(13)} \emph{Delay}: TOL "delay" channels;

\bold{(14)} \emph{On}: TOL "on" channels;

\bold{(15)} \emph{Off}: TOL "off" channels;

\bold{(16)} \emph{LightSource}: Light source, includes: \code{None}, \code{Lamp}, \code{IRDiodes},\cr 
  \code{CalibraitionLED}, \code{BlueDiodes}, \code{WhiteLight}, \code{GreenLaser}, \code{IRLaser};

\bold{(17)} \emph{AnTemp}: Annealing temperature;

\bold{(18)} \emph{TimeSinceIrr}: Time since irradiation;

\bold{(19)} \emph{Time}: Data collection time;

\bold{(20)} \emph{Date}: Data collection date.

}

\value{
Return an invisible \link{list} of S3 class object \code{"loadBIN"} containing the following elements:
\item{records}{a \link{list} containing loaded data records}
\item{tab}{a table (\link{data.frame}) summarizing items of loaded data records}
}

\references{
Duller GAT, 2016. Analyst (v4.31.9), User Mannual (\url{http://www.nutech.dtu.dk/english/Products-and-Services/Dosimetry/Radiation-Measurement-Instruments/TL_OSL_reader/Manuals}).
}

\note{
We would like to appreciate Dr Lei Gao who prompts me to write this function and provides measured data sets to test this procedure.
}

\seealso{
  \link{pickBINdata}; \link{analyseBINdata}; \link{BIN}
}

\examples{
   ### Not run.
   ### Ensure that file "foo.bin" is available 
   ### from the current working directory.
   # obj_loadBIN <- loadBINdata("foo.bin", view=TRUE)
   # class(obj_loadBIN)
   # obj_loadBIN$records
}

\keyword{data importing}
\keyword{BIN file}
