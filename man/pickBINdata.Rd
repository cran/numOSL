\name{pickBINdata}
\alias{pickBINdata}
\alias{pickBINdata.default}

\title{
BIN data set selection
}
\description{
Extracting data sets from a loaded BIN (BINX) file. 
}
\usage{
pickBINdata(obj_loadBIN, Position = NULL, Grain = NULL, Run = NULL, 
            Set = NULL, DType = NULL, IRRTime = NULL, NPoints = NULL,
            LType = NULL, Low = NULL, High = NULL, Rate = NULL, 
            Temperature = NULL, Delay = NULL, On = NULL, Off = NULL,
            LightSource = NULL, AnTemp = NULL, TimeSinceIrr = NULL, 
            view = TRUE, manual.select = FALSE, force.matrix = FALSE)
}


\arguments{
  \item{obj_loadBIN}{\link{list}(\bold{required}): an object of S3 class \code{"loadBIN"} produced by\cr 
function \link{loadBINdata}}

  \item{Position}{\link{vector}(optional): carousel position, Example: \code{Position=1:48}}

  \item{Grain}{\link{vector}(optional): grain number}

  \item{Run}{\link{vector}(optional): run number}

  \item{Set}{\link{vector}(optional): set number}

  \item{DType}{\link{character}(optional): a character vector indicating data type,\cr
Example: \code{DType=c("Natural", "N+dose")}}

  \item{IRRTime}{\link{vector}(optional): irradiation time}

  \item{NPoints}{\link{vector}(optional): number of data points}

  \item{LType}{\link{character}(optional): a character vector indicating luminescence types,\cr 
Example: \code{LType="OSL"}}

  \item{Low}{\link{vector}(optional): lower limit on temperature, time, or wavelength}

  \item{High}{\link{vector}(optional): upper limit on temperature, time, or wavelength}

  \item{Rate}{\link{vector}(optional): increasing rate of temperature, time, or wavelength}

  \item{Temperature}{\link{vector}(optional): a vector indicating the sample temperatures}

  \item{Delay}{\link{vector}(optional): TOL "delay" channels}

  \item{On}{\link{vector}(optional): TOL "on" channels}

  \item{Off}{\link{vector}(optional): TOL "off" channels}

  \item{LightSource}{\link{character}(optional): a character vector indicating light source,\cr
Example: \code{LightSource="BlueDiodes"}}

  \item{AnTemp}{\link{vector}(optional): annealing temperature}

  \item{TimeSinceIrr}{\link{vector}(optional): time since irradiation}

  \item{view}{\link{logical}(with default): logical value indicating if the loaded data should be\cr
 visualized in a Summary Table}

  \item{manual.select}{\link{logical}(with default): logical value indicating if the loaded data should be\cr
 selected manually using a Summary Table}

  \item{force.matrix}{\link{logical}(with default): logical value indicating if the picked data sets of an aliquot (grain) should be transformed into a matrix}

}


\details{
Function \link{pickBINdata} is used for pick up data sets from an object of S3 class \code{"loadBIN"} generated using function \link{loadBINdata}.
Set \code{force.matrix=TRUE} will transform data sets of an aliquot (grain) into a matrix, the transformation fails if data sets have different x (temperature, time, or wavelength) coordinates.
Transformed data sets stored in a matrix can be visualize via \link{matplot} (see examples).  
}

\value{
Return an invisible \link{list} of S3 class object \code{"pickBIN"} containing two elements:
  \item{BINdata}{selected BIN data}
  \item{agID}{Aliquot or grain ID (i.e., \code{c("NO","Position","Grain")}) of selected data sets, it returns \code{NULL} if \code{force.matrix=TRUE}}
}

\references{
Duller GAT, 2016. Analyst (v4.31.9), User Mannual.
}

\seealso{
  \link{loadBINdata}; \link{analyseBINdata}; \link{BIN};\cr
  \link{decomp}; \link{fastED}
}

\examples{
   ### Example 1 (visualize decay curves for Position=2):
   data(BIN)
   obj_pickBIN <- pickBINdata(BIN, Position=2, view=FALSE,
                              LType="OSL", force.matrix=TRUE)
   matplot(obj_pickBIN$BINdata[[1]][,1], 
           obj_pickBIN$BINdata[[1]][,-1], 
           main="Decay curves",
           xlab="Stimulation time (s)",
           ylab="Photon counts",
           type="l", log="xy")

  ### Example 2 (visualize test-dose decay curves for Position=2):
  obj_pickBIN <- pickBINdata(BIN, Position=2, Run=c(5,11,17,23,29,34,40), 
                             view=FALSE, LType="OSL", force.matrix=TRUE)
  matplot(obj_pickBIN$BINdata[[1]][,1], 
          obj_pickBIN$BINdata[[1]][,-1], 
          main="Test-dose decay curves",
          xlab="Stimulation time (s)",
          ylab="Photon counts",
          type="l", log="xy")

  ### Example 3 (visualize regenerative-dose decay curves for Position=2):
  obj_pickBIN <- pickBINdata(BIN, Position=2, Run=c(8,14,20,26,31,37), 
                             view=FALSE, LType="OSL", force.matrix=TRUE)
  matplot(obj_pickBIN$BINdata[[1]][,1], 
          obj_pickBIN$BINdata[[1]][,-1], 
          main="Regenerative-dose decay curves",
          xlab="Stimulation time (s)",
          ylab="Photon counts",
          type="l", log="xy")

}

\keyword{BIN file}
\keyword{data extraction}
